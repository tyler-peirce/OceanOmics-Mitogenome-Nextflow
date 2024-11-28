#!/bin/bash

USER=tpeirce
RUNDIR=/scratch/pawsey0812/$USER/MITOGENOMES/ilmn

for DIR in $RUNDIR/OG*/OG*; do
  echo "$DIR"
  rm $DIR/emma_prefix.txt
  rm -r $DIR/mtdna/*extended_spades
  rm -r $DIR/mtdna/*seed
  rm -r $DIR/emma/proteins
  rm $DIR/emma/cds/*ND*
  rm $DIR/emma/cds/*CYTB*
  rm $DIR/emma/cds/*CO2*
  rm $DIR/emma/cds/*CO3*
  rm $DIR/emma/cds/*ATP*
  rm $DIR/mtdna/*extended*.fq
  rm $DIR/mtdna/prefix.txt
  rm $DIR/lca/blast.*.emma102.tsv
done