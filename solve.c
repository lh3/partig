#include <assert.h>
#include "ptpriv.h"
#include "gfa-priv.h"

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
	return group;
}

void pt_cc(pt_match_t *ma)
{
	ma->cc = pt_cc_core(ma);
}

static inline uint64_t kr_splitmix64(uint64_t x)
{
	uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
	z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
	z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
	return z ^ (z >> 31);
}

static inline double kr_drand_r(uint64_t *x)
{
    union { uint64_t i; double d; } u;
	*x = kr_splitmix64(*x);
    u.i = 0x3FFULL << 52 | (*x) >> 12;
    return u.d - 1.0;
}

typedef struct {
	uint32_t m, n;
	uint64_t *a;
	int8_t *s;
} solve_aux_t;

static int64_t pt_score(const pt_match_t *ma, solve_aux_t *aux)
{
	uint32_t i;
	int64_t z = 0;
	for (i = 0; i < aux->n; ++i) {
		const pt_match1_t *m = &ma->ma[(uint32_t)aux->a[i]];
		z += -(int64_t)m->m * aux->s[m->sid[0]] * aux->s[m->sid[1]];
	}
	return z;
}

static void ks_shuffle_uint32_t(size_t n, uint32_t a[], uint64_t *x)
{
	size_t i, j;
	for (i = n; i > 1; --i) {
		uint32_t tmp;
		j = (size_t)(kr_drand_r(x) * i);
		tmp = a[j]; a[j] = a[i-1]; a[i-1] = tmp;
	}
}

uint32_t pt_solve1(const pt_match_t *ma, uint32_t off, uint32_t cnt, uint64_t *x, solve_aux_t *aux)
{
	uint32_t i, n_iter = 0, *b;
	int64_t sc_ori, sc_opt;
	if (cnt < 2) return 0;

	// find the initial states
	aux->n = 0;
	for (i = 0; i < cnt; ++i) {
		uint32_t k = (uint32_t)ma->cc[off + i];
		uint32_t o = ma->idx[k] >> 32;
		uint32_t n = (uint32_t)ma->idx[k], j;
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
	sc_ori = pt_score(ma, aux);
	//fprintf(stderr, "X\t%d\t%d\t%d\t%ld\n", off, cnt, aux->n, (long)sc_ori);
	if (cnt == 2) return 0;
	PT_MALLOC(b, cnt);
	for (i = 0; i < cnt; ++i)
		b[i] = (uint32_t)ma->cc[off + i];
	while (1) {
		uint32_t n_flip = 0;
		++n_iter;
		ks_shuffle_uint32_t(cnt, b, x);
		for (i = 0; i < cnt; ++i) {
			uint32_t k = b[i];
			uint32_t o = ma->idx[k] >> 32;
			uint32_t n = (uint32_t)ma->idx[k], j;
			uint32_t z[2];
			int8_t s;
			z[0] = z[1] = 0;
			for (j = 0; j < n; ++j) {
				const pt_match1_t *m = &ma->ma[o + j];
				assert(m->sid[0] == k);
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
	free(b);
	sc_opt = pt_score(ma, aux);
	fprintf(stderr, "[%s] group:%d, size:%d, #edges:%d, #iter:%d, sc_ori:%ld, sc_opt:%ld\n", __func__,
			(uint32_t)(ma->cc[off]>>32), cnt, aux->n, n_iter, (long)sc_ori, (long)sc_opt);
	return n_iter;
}

int8_t *pt_solve_core(const pt_match_t *ma, uint64_t x)
{
	int8_t *s;
	uint32_t st, i;
	solve_aux_t *aux;
	PT_CALLOC(aux, 1);
	PT_CALLOC(aux->s, ma->n_seg);
	for (st = 0, i = 1; i <= ma->n_seg; ++i) {
		if (i == ma->n_seg || ma->cc[st]>>32 != ma->cc[i]>>32) {
			if (i - st >= 2)
				pt_solve1(ma, st, i - st, &x, aux);
			st = i;
		}
	}
	s = aux->s;
	free(aux->a);
	free(aux);
	return s;
}

void pt_solve(pt_match_t *ma, uint64_t x)
{
	uint32_t i;
	int8_t *s;
	pt_cc(ma);
	s = pt_solve_core(ma, x);
	for (i = 0; i < ma->n_seg; ++i)
		ma->info[i].s = s[i];
	free(s);
}
