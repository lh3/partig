#include <stdlib.h>
#include <stdio.h>
#include "partig.h"
#include "ketopt.h"
#include "gfa.h"

static void print_usage(FILE *fp, const pt_pdopt_t *po)
{
	fprintf(stderr, "Usage: partig [options] <in.gfa>\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "  -k INT     k-mer size [%d]\n", po->k);
	fprintf(stderr, "  -w INT     minimizer window size [%d]\n", po->w);
	fprintf(stderr, "  -c INT     max occurrance [%d]\n", po->max_occ);
	fprintf(stderr, "  -m FLOAT   mini k-mer similarity [%.3g]\n", po->min_sim);
}

int main(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int i, c;
	pt_match_t *ma;
	gfa_t *g;
	pt_pdopt_t po;

	pt_realtime();
	pt_pdopt_init(&po);
	while ((c = ketopt(&o, argc, argv, 1, "k:w:c:m:", 0)) >= 0) {
		if (c == 'k') po.k = atoi(o.arg);
		else if (c == 'w') po.w = atoi(o.arg);
		else if (c == 'c') po.max_occ = atoi(o.arg);
		else if (c == 'm') po.min_sim = atof(o.arg);
	}
	if (o.ind == argc) {
		print_usage(stderr, &po);
		return 0;
	}

	g = gfa_read(argv[o.ind]);
	if (pt_verbose >= 3)
		fprintf(stderr, "[%s::%.3f] read the graph\n", __func__, pt_realtime());

	ma = pt_pdist(&po, g);
	/*
	pt_partition(&so, ma);
	if (fn_link) {
		pt_mcgraph_t *mcg;
		mcg = pt_phase_gen_graph(g, ma, fn_link, weight);
	}
	*/
	pt_match_print(stdout, g, ma);

	pt_match_free(ma);
	gfa_destroy(g);

	if (pt_verbose >= 3) {
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, PT_VERSION);
		fprintf(stderr, "[M::%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, pt_realtime(), pt_cputime(), pt_peakrss() / 1024.0 / 1024.0 / 1024.0);
	}
	return 0;
}
