This is an *experimental* tool to estimate the approximate distances between
all pairs of unitigs. It takes a GFA or FASTA file as input and outputs a
TAB-delimited format as follows:
```txt
C  seqName   seqLen    minimizerConsidered    mzUnique
S  seqName1  seqName2  strand  mzConsidered1  mzConsidered2  mzShared  kmerSimilarity
```
The algorithm behind is fairly simple: collect all minimizers, sort them,
identify the longest colinear chain based on longest increasing sequence, and
finally estimate k-mer similarity. It can find the similarity between 6Gb
of unitigs in four minutes on a single thread.

This mini-project was initially intended for unitig partition. That is how this
repo is named. It has demonstrated the feasibility of similarity based
partition. The basic idea will be folded into hifiasm. The experimental code is
left here. It probably won't receive update but it can be helpful if you want
to quickly get an idea of redundancy or to learn how to do that efficiently.
