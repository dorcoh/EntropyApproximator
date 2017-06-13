#!/bin/bash
#PBS -l walltime=04:00:00,nodes=1
#PBS -j oe
#PBS -N batchtest
#PBS -q default


set -x

cd $PBS_O_WORKDIR

cp $basedir/../srcILOG/WH_cplex $TMPDIR/iloglue

cd $TMPDIR

timeout $timeout ./iloglue -paritylevel 1 -number $level $basedir/$file > out

# Copy output files to your output folder
cp -f $TMPDIR/out $outdir/`basename $file`.xor$level.loglen$len.$sample.ILOGLUE.uai.LOG