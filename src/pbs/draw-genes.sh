#!/usr/bin/env bash

cd $HOME/workspace/gbm
for idx in raw int; do
    qsub <<EOF
        #PBS -N $idx
        #PBS -l walltime=240:0:0
        #PBS -l nodes=comput3:ppn=1
        #PBS -e ./log/${idx}.draw.stderr.log
        #PBS -o ./log/${idx}.draw.stdout.log

        cd $HOME/workspace/gbm
        Rscript src/R/seurat.draw.genes.R $idx
EOF
done
