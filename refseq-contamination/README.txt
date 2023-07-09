
GENERATING DATA
--------------------

Scripts run-kmers-chm13.sh and run-kmers-grch38.sh will run 31 jobs each to
compute 31 sets of outputs in chm13v2/ and grch38/.  Slurm parameters are in
the script itself.

  sbatch run-kmers-chm13.sh
  sbatch run-kmers-grch38.sh

Outputs:

  hpq - hits per query
      - col 1 - the number of kmers in the reference that match the query
      - col 2 - the number of kmers in the query that match the reference
      - col 3 - the length of the query
      - col 4 - the name of the query

  mpb - the number of kmers in the queries that match each position of the
        reference
      - col 1 position in chrY
      - col 2 number of kmers (counting duplicates) in the queries that match

  qpb - the number of queries that have at least one matching kmer for each
        position of the reference
      - col 1 position in chrY
      - col 2 number of queries that have a kmer match to this position

SUMMARIZING DATA
--------------------

The hpq files can be simply concatenated.

  cat chm13v2/hpq.?? > chm13v2-hits-per-query
  cat grch38/hpq.?? > grch38-hits-per-query


Script with-without.pl computes the two histograms of query length with and
without a hit to chrY.  Both are on stdout.  The first histogram is
logarithmic in query size, the second is for 1 Kbp windows.

  perl with-without.pl > with-without.out


The per-base results need a helper script to sum the values per base.

  perl combine.pl chm13v2/mpb.?? > chm13v2-mers-per-base
  perl combine.pl chm13v2/qpb.?? > chm13v2-queries-per-base

  perl combine.pl grch38/mpb.?? > grch38-mers-per-base
  perl combine.pl grch38/qpb.?? > grch38-queries-per-base


The per-base results can also be summed to find min/ave/max per 10k window.

  perl coverage.pl chm13v2/mpb.?? > chm13v2-mers-per-10k-window
  perl coverage.pl chm13v2/qpb.?? > chm13v2-queries-per-10k-window

  perl coverage.pl grch38/mpb.?? > grch38-mers-per-10k-window
  perl coverage.pl grch38/qpb.?? > grch38-queries-per-10k-window


EXTRACTING READS
--------------------

  awk '{ if ($1 > 0) print }' < chm13v2-hits-per-query > chm13v2-hits-per-query.with-hit

  ../seqrequester/build/bin/seqrequester extract \
    -sequences chm13v2-hits-per-query.with-hit release213/bacteria/*gz\
> chm13v2-hits-per-query.with-hit.fasta
