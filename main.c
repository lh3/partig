#include <stdlib.h>
#include <stdio.h>
#include "partig.h"
#include "ketopt.h"
#include "gfa.h"

static void print_usage(FILE *fp, const pt_pdopt_t *po, const pt_svopt_t *so)
{
	fprintf(stderr, "Usage: partig [options] <in.gfa>\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "  -k INT     k-mer size [%d]\n", po->k);
	fprintf(stderr, "  -w INT     minimizer window size [%d]\n", po->w);
	fprintf(stderr, "  -c INT     max occurrance [%d]\n", po->max_occ);
	fprintf(stderr, "  -m FLOAT   mini k-mer similarity [%.3g]\n", po->min_sim);
	fprintf(stderr, "  -n INT     inspect top INT edges [%d]\n", so->topn);
	fprintf(stderr, "  -s INT     RNG seed [%ld]\n", (long)so->seed);
	fprintf(stderr, "  -p INT     rounds of perturbations [%d]\n", so->n_perturb);
	fprintf(stderr, "  -f FLOAT   fraction to flip for perturbation [%.3g]\n", so->f_perturb);
	fprintf(stderr, "  -l FILE    unitig links []\n");
}

int main(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int i, c;
	pt_match_t *ma;
	gfa_t *g;
	pt_pdopt_t po;
	pt_svopt_t so;
	char *fn_link = 0;
	pt128_t *link = 0;
	uint32_t n_link = 0;

	pt_realtime();
	pt_pdopt_init(&po);
	pt_svopt_init(&so);
	while ((c = ketopt(&o, argc, argv, 1, "k:w:c:n:s:p:f:m:l:", 0)) >= 0) {
		if (c == 'k') po.k = atoi(o.arg);
		else if (c == 'w') po.w = atoi(o.arg);
		else if (c == 'c') po.max_occ = atoi(o.arg);
		else if (c == 'm') po.min_sim = atof(o.arg);
		else if (c == 'n') so.topn = atoi(o.arg);
		else if (c == 's') so.seed = atol(o.arg);
		else if (c == 'p') so.n_perturb = atoi(o.arg);
		else if (c == 'f') so.f_perturb = atof(o.arg);
		else if (c == 'l') fn_link = o.arg;
	}
	if (o.ind == argc) {
		print_usage(stderr, &po, &so);
		return 0;
	}

	g = gfa_read(argv[o.ind]);
	if (pt_verbose >= 3)
		fprintf(stderr, "[%s::%.3f] read the graph\n", __func__, pt_realtime());

	if (fn_link) {
		link = pt_read_links(g, fn_link, &n_link);
		if (pt_verbose >= 3)
			fprintf(stderr, "[%s::%.3f] read %d links\n", __func__, pt_realtime(), n_link);
	}

	ma = pt_pdist(&po, g);
	pt_solve(&so, ma);
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
