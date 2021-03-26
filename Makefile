CC=			gcc
CFLAGS=		-g -Wall -Wc++-compat -std=c99 -O2
CPPFLAGS=
INCLUDES=
OBJS=		kalloc.o gfa-base.o gfa-io.o misc.o sketch.o pdist.o solve.o
PROG=		partig
LIBS=		-lz -lpthread -lm

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

.SUFFIXES:.c .o
.PHONY:all clean depend

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

partig:$(OBJS) main.o
		$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

clean:
		rm -fr gmon.out *.o a.out $(PROG) *~ *.a *.dSYM

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE

gfa-base.o: gfa-priv.h gfa.h kstring.h khash.h kalloc.h ksort.h
gfa-io.o: kstring.h gfa-priv.h gfa.h kseq.h
kalloc.o: kalloc.h
main.o: partig.h gfa.h ketopt.h
misc.o: ptpriv.h partig.h gfa.h ksort.h
pdist.o: ptpriv.h partig.h gfa.h ksort.h gfa-priv.h
sketch.o: partig.h gfa.h kvec.h
solve.o: ptpriv.h partig.h gfa.h
