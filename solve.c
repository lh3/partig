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
	ma->group = pt_cc_core(ma);
}

static inline uint64_t kr_splitmix64(uint64_t x)
{
	uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
	z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
	z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
	return z ^ (z >> 31);
}

typedef struct {
	uint32_t m, n;
	uint64_t *a;
	int8_t *s;
} solve_aux_t;

void pt_solve1(const pt_match_t *ma, uint32_t off, uint32_t cnt, uint64_t *x, solve_aux_t *aux)
{
	uint32_t i;
	if (cnt < 2) return;

	// find the initial states
	aux->n = 0;
	for (i = 0; i < cnt; ++i) {
		uint32_t o = ma->idx[off + i] >> 32;
		uint32_t n = (uint32_t)ma->idx[off + i], j;
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
			aux->s[m->sid[0]] = *x&1? 1 : -1;
			*x = kr_splitmix64(*x);
			aux->s[m->sid[1]] = -aux->s[m->sid[0]];
		} else if (aux->s[m->sid[0]] == 0) {
			aux->s[m->sid[0]] = -aux->s[m->sid[1]];
		} else if (aux->s[m->sid[1]] == 0) {
			aux->s[m->sid[1]] = -aux->s[m->sid[0]];
		}
	}
}

int8_t *pt_solve_core(const pt_match_t *ma, uint64_t x)
{
	int8_t *s;
	uint32_t st, i;
	solve_aux_t *aux;
	PT_CALLOC(aux, 1);
	PT_CALLOC(aux->s, ma->n_seg);
	x = kr_splitmix64(x);
	for (st = 0, i = 1; i <= ma->n_seg; ++i) {
		if (i == ma->n_seg || ma->group[st]>>32 != ma->group[i]>>32) {
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
	pt_cc(ma);
	ma->s = pt_solve_core(ma, x);
}
