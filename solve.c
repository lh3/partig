#include "ptpriv.h"
#include "gfa-priv.h"

uint64_t *pt_cc(pt_match_t *ma)
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
