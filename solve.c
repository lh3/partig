#include <assert.h>
#include <string.h>
#include "ptpriv.h"
#include "gfa-priv.h"

void pt_svopt_init(pt_svopt_t *opt)
{
	memset(opt, 0, sizeof(pt_svopt_t));
	opt->seed = 11;
	opt->topn = 1<<30;
	opt->n_perturb = 1000;
	opt->f_perturb = 0.1;
}

uint64_t *pt_cc_core(const pt_match_t *ma)
{
	uint32_t i, x, y, *flag, *stack = 0, ns = 0, ms = 0;
	uint64_t *group;

	PT_MALLOC(flag, ma->n_seg);
	for (i = 0; i < ma->n_seg; ++i)
		flag[i] = (uint32_t)-1;

	// connected componets
	for (i = 0; i < ma->n_seg; ++i) {
		if (flag[i] != (uint32_t)-1) continue;
		if (ns == ms) PT_EXPAND(stack, ms);
		stack[ns++] = i;
		while (ns > 0) {
			uint32_t k, j, n, s;
			k = stack[--ns];
			flag[k] = i;
			n = (uint32_t)ma->idx[k];
			s = ma->idx[k] >> 32;
			for (j = 0; j < n; ++j) {
				uint32_t t = ma->ma[s + j].sid[1];
				if (flag[t] != (uint32_t)-1) continue;
				if (ns == ms) PT_EXPAND(stack, ms);
				stack[ns++] = t;
			}
		}
	}
	free(stack);

	// precalculate the size of each group
	PT_CALLOC(group, ma->n_seg);
	for (i = 0; i < ma->n_seg; ++i)
		group[i] = (uint64_t)flag[i] << 32 | i;
	radix_sort_gfa64(group, group + ma->n_seg);
	for (i = 1, x = y = 0; i <= ma->n_seg; ++i) {
		if (i == ma->n_seg || group[i]>>32 != group[x]>>32) {
			uint32_t j;
			for (j = x; j < i; ++j)
				group[j] = (uint64_t)y << 32 | (uint32_t)group[j];
			++y, x = i;
		}
	}
	free(flag);
	if (pt_verbose >= 3)
		fprintf(stderr, "[%s::%.3f] identified %d connected components\n", __func__, pt_realtime(), y);
	return group;
}

void pt_cc(pt_match_t *ma)
{
	ma->cc = pt_cc_core(ma);
}

typedef struct {
	uint32_t m, n, *shuffled;
	uint32_t off, size; // offset in pt_match_t::cc; size of the component
	uint64_t *a, *buf;
	int8_t *s, *s_tmp;
} solve_aux_t;

static int64_t pt_score(const pt_match_t *ma, int32_t topn, solve_aux_t *aux)
{
	uint32_t i;
	int64_t z = 0;
	for (i = 0; i < aux->size; ++i) {
		uint32_t k = (uint32_t)ma->cc[aux->off + i];
		uint32_t o = ma->idx[k] >> 32;
		uint32_t n = (uint32_t)ma->idx[k], j;
		for (j = 0; j < n; ++j)
			aux->buf[j] = (uint64_t)((uint32_t)-1 - ma->ma[o + j].m) << 32 | (o + j);
		radix_sort_gfa64(aux->buf, aux->buf + n);
		for (j = 0; j < n && j < topn; ++j) {
			const pt_match1_t *m = &ma->ma[(uint32_t)aux->buf[j]];
			z += -(int64_t)m->m * aux->s[m->sid[0]] * aux->s[m->sid[1]];
		}
	}
	return z;
}

static int64_t pt_solve1_init_phase(const pt_match_t *ma, int32_t topn, uint64_t *x, solve_aux_t *aux)
{
	uint32_t i;
	aux->n = 0;
	for (i = 0; i < aux->size; ++i) {
		uint32_t k = (uint32_t)ma->cc[aux->off + i];
		uint32_t o = ma->idx[k] >> 32;
		uint32_t n = (uint32_t)ma->idx[k], j;
		aux->shuffled[i] = k;
		for (j = 0; j < n; ++j) {
			const pt_match1_t *m = &ma->ma[o + j];
			if (aux->n == aux->m) PT_EXPAND(aux->a, aux->m);
			aux->a[aux->n++] = (uint64_t)((uint32_t)-1 - m->m) << 32 | (o + j);
		}
	}
	radix_sort_gfa64(aux->a, aux->a + aux->n);
	for (i = 0; i < aux->n; ++i) { // from the strongest edge to the weakest
		const pt_match1_t *m = &ma->ma[(uint32_t)aux->a[i]];
		if (aux->s[m->sid[0]] == 0 && aux->s[m->sid[1]] == 0) {
			*x = kr_splitmix64(*x);
			aux->s[m->sid[0]] = *x&1? 1 : -1;
			aux->s[m->sid[1]] = -aux->s[m->sid[0]];
		} else if (aux->s[m->sid[0]] == 0) {
			aux->s[m->sid[0]] = -aux->s[m->sid[1]];
		} else if (aux->s[m->sid[1]] == 0) {
			aux->s[m->sid[1]] = -aux->s[m->sid[0]];
		}
	}
	return pt_score(ma, topn, aux);
}

static void pt_solve1_perturb(const pt_svopt_t *opt, const pt_match_t *ma, uint64_t *x, solve_aux_t *aux)
{
	uint32_t i;
	double y;
	for (i = 0; i < aux->size; ++i) {
		uint32_t k = (uint32_t)ma->cc[aux->off + i];
		y = kr_drand_r(x);
		if (y < opt->f_perturb)
			aux->s[k] = -aux->s[k];
	}
}

static int64_t pt_solve1_optimize(const pt_match_t *ma, uint32_t topn, uint64_t *x, solve_aux_t *aux, uint32_t *n_iter)
{
	uint32_t i;
	while (1) {
		uint32_t n_flip = 0;
		++(*n_iter);
		ks_shuffle_uint32_t(aux->size, aux->shuffled, x);
		for (i = 0; i < aux->size; ++i) {
			uint32_t k = aux->shuffled[i];
			uint32_t o = ma->idx[k] >> 32;
			uint32_t n = (uint32_t)ma->idx[k], j;
			uint32_t z[2];
			int8_t s;
			for (j = 0; j < n; ++j) {
				const pt_match1_t *m = &ma->ma[o + j];
				assert(m->sid[0] == k);
				aux->buf[j] = (uint64_t)((uint32_t)-1 - m->m) << 32 | (o + j);
			}
			radix_sort_gfa64(aux->buf, aux->buf + n);
			for (j = 0, z[0] = z[1] = 0; j < n && j < topn; ++j) {
				const pt_match1_t *m = &ma->ma[(uint32_t)aux->buf[j]];
				if (aux->s[m->sid[1]] > 0) z[0] += m->m;
				else if (aux->s[m->sid[1]] < 0) z[1] += m->m;
			}
			if (z[0] == z[1]) continue;
			s = z[0] > z[1]? -1 : 1;
			if (aux->s[k] != s)
				aux->s[k] = s, ++n_flip;
		}
		if (n_flip == 0) break;
	}
	return pt_score(ma, topn, aux);
}

uint32_t pt_solve1(const pt_svopt_t *opt, const pt_match_t *ma, uint64_t *x, solve_aux_t *aux)
{
	uint32_t j, k, n_iter = 0;
	int64_t sc_ori, sc_opt = -(1<<30), sc;
	if (aux->size < 2) return 0;

	// first guess
	sc_ori = pt_solve1_init_phase(ma, opt->topn, x, aux);
	if (aux->size == 2) return 0;

	// optimize
	sc_opt = pt_solve1_optimize(ma, opt->topn, x, aux, &n_iter);
	for (j = 0; j < aux->size; ++j)
		aux->s_tmp[aux->shuffled[j]] = aux->s[aux->shuffled[j]];
	for (k = 0; k < opt->n_perturb; ++k) {
		pt_solve1_perturb(opt, ma, x, aux);
		sc = pt_solve1_optimize(ma, opt->topn, x, aux, &n_iter);
		if (sc > sc_opt) {
			for (j = 0; j < aux->size; ++j)
				aux->s_tmp[aux->shuffled[j]] = aux->s[aux->shuffled[j]];
			sc_opt = sc;
		} else {
			for (j = 0; j < aux->size; ++j)
				aux->s[aux->shuffled[j]] = aux->s_tmp[aux->shuffled[j]];
		}
	}
	for (j = 0; j < aux->size; ++j)
		aux->s[aux->shuffled[j]] = aux->s_tmp[aux->shuffled[j]];
	fprintf(stderr, "[%s::%.3f] group:%d, size:%d, #edges:%d, #iter:%d, sc_ori:%ld, sc_opt:%ld\n", __func__, pt_realtime(),
			(uint32_t)(ma->cc[aux->off]>>32), aux->size, aux->n, n_iter, (long)sc_ori, (long)sc_opt);
	return n_iter;
}

int8_t *pt_partition_core(const pt_svopt_t *opt, const pt_match_t *ma)
{
	int8_t *s;
	uint32_t st, i, max = 0;
	uint64_t x = opt->seed;
	solve_aux_t *aux;
	PT_CALLOC(aux, 1);
	PT_CALLOC(aux->s, ma->n_seg);
	PT_CALLOC(aux->s_tmp, ma->n_seg);
	for (i = 0; i < ma->n_seg; ++i) {
		uint32_t n = (uint32_t)ma->idx[i];
		max = max > n? max : n;
	}
	PT_MALLOC(aux->buf, max);
	PT_MALLOC(aux->shuffled, ma->n_seg); // FIXME: this is over-allocation for convenience
	for (st = 0, i = 1; i <= ma->n_seg; ++i) {
		if (i == ma->n_seg || ma->cc[st]>>32 != ma->cc[i]>>32) {
			if (i - st >= 2) {
				aux->off = st, aux->size = i - st;
				pt_solve1(opt, ma, &x, aux);
			}
			st = i;
		}
	}
	s = aux->s;
	free(aux->a); free(aux->buf); free(aux->shuffled); free(aux->s_tmp);
	free(aux);
	if (pt_verbose >= 3)
		fprintf(stderr, "[%s::%.3f] finished partitioning\n", __func__, pt_realtime());
	return s;
}

void pt_partition(const pt_svopt_t *opt, pt_match_t *ma)
{
	uint32_t i;
	int8_t *s;
	uint64_t *buf;
	pt_cc(ma);
	s = pt_partition_core(opt, ma);
	PT_MALLOC(buf, ma->n_ma); // FIXME: this is over-allocation for convenience
	for (i = 0; i < ma->n_seg; ++i) {
		uint32_t z[2];
		uint32_t o = ma->idx[i] >> 32;
		uint32_t n = (uint32_t)ma->idx[i], j;
		ma->info[i].s = s[i];
		for (j = 0; j < n; ++j) {
			const pt_match1_t *m = &ma->ma[o + j];
			buf[j] = (uint64_t)((uint32_t)-1 - m->m) << 32 | (o + j);
		}
		radix_sort_gfa64(buf, buf + n);
		for (j = 0, z[0] = z[1] = 0; j < n; ++j) {
			const pt_match1_t *m = &ma->ma[(uint32_t)buf[j]];
			if (s[m->sid[1]] > 0) z[0] += m->m;
			else if (s[m->sid[1]] < 0) z[1] += m->m;
		}
		ma->info[i].m[0] = z[0], ma->info[i].m[1] = z[1];
	}
	free(buf);
	free(s);
}
