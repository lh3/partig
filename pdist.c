#include "ptpriv.h"
#include "ksort.h"

#define mz_key(z) ((z).x)
KRADIX_SORT_INIT(mz, pt_mz1_t, mz_key, 8)

void pt_opt_init(pt_pdopt_t *opt)
{
	memset(opt, 0, sizeof(pt_pdopt_t));
	opt->k = 51, opt->w = 51, opt->is_hpc = 1;
	opt->max_occ = 20;
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

static pt128_t *pt_collect_anchors(uint32_t n_mz, const pt_mz1_t *mz, int32_t max_occ, int64_t *n_a_, uint32_t *cnt)
{
	uint32_t st, j;
	int64_t m_a = 0, n_a = 0;
	pt128_t *a = 0;
	for (j = 1, st = 0; j <= n_mz; ++j) {
		if (j == n_mz || mz[j].x != mz[st].x) {
			uint32_t k, l;
			if (j - st == 1 || j - st > max_occ) goto end_anchor;
			for (k = st; k < j; ++k) {
				++cnt[mz[k].rid];
				for (l = k + 1; l < j; ++l) {
					uint32_t span, rev = (mz[k].rev != mz[l].rev);
					if (n_a == m_a) PT_EXPAND(a, m_a);
					span = mz[l].span;
					a[n_a].x = (uint64_t)mz[k].rid << 33 | mz[l].rid << 1 | rev;
					a[n_a++].y = (uint64_t)mz[k].pos << 32 | (mz[l].pos + 1 - span);
					if (n_a == m_a) PT_EXPAND(a, m_a);
					span = mz[k].span;
					a[n_a].x = (uint64_t)mz[l].rid << 33 | mz[k].rid << 1 | rev;
					a[n_a++].y = (uint64_t)mz[l].rid << 33 | (mz[k].pos + 1 - span);
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

void pt_pdist(const pt_pdopt_t *opt, const gfa_t *g)
{
	int64_t n_an;
	uint32_t n_mz, *cnt;
	pt_mz1_t *mz;
	pt128_t *an = 0;

	mz = pt_collect_minimizers(opt, g, &n_mz);
	PT_CALLOC(cnt, g->n_seg);
	an = pt_collect_anchors(n_mz, mz, opt->max_occ, &n_an, cnt);
	free(mz);

	free(cnt);
	free(an);
}
