#include <math.h>
#include "ptpriv.h"
#include "ksort.h"
#include "gfa-priv.h"

#define mz_key(z) ((z).x)
KRADIX_SORT_INIT(mz, pt_mz1_t, mz_key, 8)

void pt_opt_init(pt_pdopt_t *opt)
{
	memset(opt, 0, sizeof(pt_pdopt_t));
	opt->k = 51, opt->w = 51, opt->is_hpc = 1;
	opt->max_occ = 20;
	opt->min_cnt = 5;
	opt->min_sim = 0.9;
	opt->diff_thres = 0.5;
}

static pt_mz1_t *pt_collect_minimizers(const pt_pdopt_t *opt, const gfa_t *g, uint32_t *n_mz_)
{
	uint32_t i;
	pt_mz1_v mz = {0,0,0};
	for (i = 0; i < g->n_seg; ++i)
		pt_sketch(g->seg[i].seq, g->seg[i].len, opt->w, opt->k, i, opt->is_hpc, &mz);
	if (pt_verbose >= 3)
		fprintf(stderr, "[%s::%.3f] collected %d minimizers\n", __func__, pt_realtime(), mz.n);
	radix_sort_mz(mz.a, mz.a + mz.n);
	if (pt_verbose >= 3)
		fprintf(stderr, "[%s::%.3f] sorted %d minimizers\n", __func__, pt_realtime(), mz.n);
	*n_mz_ = mz.n;
	return mz.a;
}

static pt128_t *pt_collect_anchors(const gfa_t *g, uint32_t n_mz, const pt_mz1_t *mz, int32_t max_occ, int64_t *n_a_, uint32_t *cnt, uint32_t *ucnt)
{
	uint32_t st, j;
	int64_t m_a = 0, n_a = 0;
	pt128_t *a = 0;
	for (j = 1, st = 0; j <= n_mz; ++j) {
		if (j == n_mz || mz[j].x != mz[st].x) {
			uint32_t k, l;
			if (j - st == 1) ++ucnt[mz[st].rid];
			if (j - st == 1 || j - st > max_occ) goto end_anchor;
			for (k = st; k < j; ++k) {
				++cnt[mz[k].rid];
				for (l = k + 1; l < j; ++l) {
					uint32_t span, rev = (mz[k].rev != mz[l].rev);
					int32_t lk = g->seg[mz[k].rid].len, ll = g->seg[mz[l].rid].len;
					if (n_a == m_a) PT_EXPAND(a, m_a);
					span = mz[l].span;
					a[n_a].x = (uint64_t)mz[k].rid << 33 | mz[l].rid << 1 | rev;
					a[n_a++].y = (uint64_t)mz[k].pos << 32 | (rev? ll - (mz[l].pos + 1 - span) - 1 : mz[l].pos);
					if (n_a == m_a) PT_EXPAND(a, m_a);
					span = mz[k].span;
					a[n_a].x = (uint64_t)mz[l].rid << 33 | mz[k].rid << 1 | rev;
					a[n_a++].y = (uint64_t)mz[l].pos << 32 | (rev? lk - (mz[k].pos + 1 - span) - 1 : mz[k].pos);
				}
			}
end_anchor:	st = j;
		}
	}
	if (pt_verbose >= 3)
		fprintf(stderr, "[%s::%.3f] collected %ld anchors\n", __func__, pt_realtime(), (long)n_a);
	radix_sort_pt128x(a, a + n_a);
	if (pt_verbose >= 3)
		fprintf(stderr, "[%s::%.3f] sorted %ld anchors\n", __func__, pt_realtime(), (long)n_a);
	*n_a_ = n_a;
	return a;
}

int32_t pt_lis_64(int32_t n, const uint64_t *a, int32_t *b)
{
	int32_t i, k, L = 0, *M, *P = b;
	PT_MALLOC(M, n+1);
	for (i = 0; i < n; ++i) {
		int32_t lo = 1, hi = L, newL;
		while (lo <= hi) {
			int32_t mid = (lo + hi + 1) >> 1;
			if (a[M[mid]] < a[i]) lo = mid + 1;
			else hi = mid - 1;
		}
		newL = lo, P[i] = M[newL - 1], M[newL] = i;
		if (newL > L) L = newL;
	}
	k = M[L];
	memcpy(M, P, n * sizeof(int32_t));
	for (i = L - 1; i >= 0; --i) b[i] = k, k = M[k];
	free(M);
	return L;
}

static pt_match1_t *pt_cal_sim(int64_t n_an, const pt128_t *an, const uint32_t *cnt, const uint32_t *ucnt, int32_t k, int32_t min_cnt, double min_sim, uint32_t *n_ma_)
{
	int64_t st, i, j, max = 0;
	uint64_t *a;
	int32_t *b;
	uint32_t n_ma = 0, m_ma = 0;
	pt_match1_t *ma = 0;
	for (st = 0, i = 1; i <= n_an; ++i) // pre-calculate the max size
		if (i == n_an || an[i].x != an[st].x) {
			max = max > i - st? max : i - st;
			st = i;
		}
	PT_MALLOC(a, max);
	PT_MALLOC(b, max);
	for (st = 0, i = 1; i <= n_an; ++i) {
		if (i == n_an || an[i].x != an[st].x) {
			int32_t l, n = 0;
			pt_match1_t m;
			if (i - st < min_cnt) goto end_chain;
			for (j = st; j < i; ++j)
				a[n++] = an[j].y;
			radix_sort_gfa64(a, a + n);
			for (l = 0; l < n; ++l) a[l] = (uint32_t)a[l];
			m.m = pt_lis_64(n, a, b);
			if (m.m < min_cnt) goto end_chain;
			m.sid[0] = an[st].x >> 33;
			m.sid[1] = ((uint32_t)an[st].x) >> 1;
			if (m.sid[0] == m.sid[1]) goto end_chain;
			m.n[0] = cnt[m.sid[0]];
			m.n[1] = cnt[m.sid[1]];
			m.uni[0] = ucnt[0];
			m.uni[1] = ucnt[1];
			m.rev = (an[st].x>>32&1) ^ (an[st].x&1);
			m.sim = pow(2.0 * m.m / (m.n[0] + m.n[1]), 1.0 / k);
			if (m.sim >= min_sim) {
				if (n_ma == m_ma) PT_EXPAND(ma, m_ma);
				ma[n_ma++] = m;
			}
end_chain:	st = i;
		}
	}
	free(b);
	free(a);
	*n_ma_ = n_ma;
	if (pt_verbose >= 3)
		fprintf(stderr, "[%s::%.3f] calculated similarity for %d pairs of unitigs\n", __func__, pt_realtime(), n_ma);
	return ma;
}

static void pt_pdist_idx(pt_match_t *ma)
{
	uint32_t st, i;
	PT_CALLOC(ma->idx, ma->n_seg);
	for (st = 0, i = 1; i <= ma->n_ma; ++i)
		if (i == ma->n_ma || ma->ma[i].sid[0] != ma->ma[st].sid[0])
			ma->idx[ma->ma[st].sid[0]] = (uint64_t)st << 32 | (i - st), st = i;
}

static uint32_t pt_pdist_mark_pair(const pt_match_t *ma, uint8_t *mark, uint32_t sid1, uint32_t sid2)
{
	uint32_t o = ma->idx[sid1] >> 32, c = 0;
	uint32_t n = (uint32_t)ma->idx[sid1], j;
	for (j = o; j < o + n; ++j)
		if (ma->ma[j].sid[1] == sid2)
			mark[j] = 1, ++c;
	return c;
}

static void pt_pdist_flt(pt_match_t *ma, int32_t min_cnt, double diff_thres)
{
	uint32_t i, k;
	uint8_t *mark;
	PT_CALLOC(mark, ma->n_ma);
	for (i = 0; i < ma->n_seg; ++i) {
		uint32_t o = ma->idx[i] >> 32;
		uint32_t n = (uint32_t)ma->idx[i];
		uint32_t j, max = 0;
		if (n == 0) continue;
		for (j = o; j < o + n; ++j)
			max = max > ma->ma[j].m? max : ma->ma[j].m;
		for (j = o; j < o + n; ++j)
			if (ma->ma[j].m >= max * diff_thres || ma->ma[j].m + min_cnt >= max)
				mark[j] = 1;
	}
	for (i = 0; i < ma->n_ma; ++i) {
		uint32_t c;
		if (mark[i] == 0) continue;
		c = pt_pdist_mark_pair(ma, mark, ma->ma[i].sid[1], ma->ma[i].sid[0]);
		if (c == 0) mark[i] = 0;
	}
	for (i = k = 0; i < ma->n_ma; ++i)
		if (mark[i]) ma->ma[k++] = ma->ma[i];
	ma->n_ma = k;
	free(mark);
	pt_pdist_idx(ma);
	if (pt_verbose >= 3)
		fprintf(stderr, "[%s::%.3f] %d pairs of unitigs remain after filtering\n", __func__, pt_realtime(), ma->n_ma);
}

pt_match_t *pt_pdist(const pt_pdopt_t *opt, const gfa_t *g)
{
	pt_match_t *ma;
	int64_t n_an;
	uint32_t n_mz;
	pt_mz1_t *mz;
	pt128_t *an = 0;

	PT_CALLOC(ma, 1);
	ma->n_seg = g->n_seg;
	mz = pt_collect_minimizers(opt, g, &n_mz);
	PT_CALLOC(ma->cnt, g->n_seg);
	PT_CALLOC(ma->ucnt, g->n_seg);
	an = pt_collect_anchors(g, n_mz, mz, opt->max_occ, &n_an, ma->cnt, ma->ucnt);
	free(mz);
	ma->ma = pt_cal_sim(n_an, an, ma->cnt, ma->ucnt, opt->k, opt->min_cnt, opt->min_sim, &ma->n_ma);
	free(an);
	pt_pdist_idx(ma);
	pt_pdist_flt(ma, opt->min_cnt, opt->diff_thres);
	return ma;
}

void pt_match_print(FILE *fp, const gfa_t *g, const pt_match_t *ma)
{
	uint32_t i;
	for (i = 0; i < g->n_seg; ++i) {
		const gfa_seg_t *s = &g->seg[i];
		fprintf(fp, "C\t%s\t%d\t%d\t%d\n", s->name, s->len, ma->cnt[i], ma->ucnt[i]);
	}
	for (i = 0; i < ma->n_ma; ++i) {
		const pt_match1_t *m = &ma->ma[i];
		fprintf(fp, "S\t%s\t%s\t%c\t%d\t%d\t%d\t%.6f\n", g->seg[m->sid[0]].name, g->seg[m->sid[1]].name,
				"+-"[!!m->rev], m->n[0], m->n[1], m->m, m->sim);
	}
}

void pt_match_free(pt_match_t *ma)
{
	free(ma->idx); free(ma->cnt); free(ma->ucnt); free(ma->ma);
	free(ma);
}
