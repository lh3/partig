#include <stdlib.h>
#include <stdio.h>
#include "partig.h"
#include "ketopt.h"
#include "gfa.h"

static void print_usage(FILE *fp)
{
	fprintf(stderr, "Usage: partig [options] <in.gfa>\n");
}

int main(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int c;
	gfa_t *g;
	pt_pdopt_t po;

	pt_realtime();
	pt_opt_init(&po);
	while ((c = ketopt(&o, argc, argv, 1, "k:w:c:", 0)) >= 0) {
		if (c == 'k') po.k = atoi(o.arg);
		else if (c == 'w') po.w = atoi(o.arg);
		else if (c == 'c') po.max_occ = atoi(o.arg);
	}
	if (o.ind == argc) {
		print_usage(stderr);
		return 0;
	}

	g = gfa_read(argv[o.ind]);
	if (pt_verbose >= 3)
		fprintf(stderr, "[%s::%.3f] read the graph\n", __func__, pt_realtime());
	pt_pdist(&po, g);
	gfa_destroy(g);
	return 0;
}
