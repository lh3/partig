#include "ptpriv.h"

void pt_mc_index(pt_mcgraph_t *mcg)
{
	uint32_t j, st;
	free(mcg->idx);
	PT_CALLOC(mcg->idx, mcg->n_node);
	for (st = 0, j = 1; j <= mcg->n_edge; ++j) {
		if (j == mcg->n_edge || mcg->edge[j].x>>32 != mcg->edge[st].x>>32) {
			mcg->idx[mcg->edge[st].x>>32] = (uint64_t)st << 32 | (j - st);
			st = j;
		}
	}
}

void pt_mc_destroy(pt_mcgraph_t *mcg)
{
	free(mcg->edge); free(mcg->idx); free(mcg->s); free(mcg);
}
