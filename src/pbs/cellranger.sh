cd $HOME/workspace/gbm

for idx in S20220222-1 S20220222-2; do

    qsub << \
EOF
    #PBS -N $idx
    #PBS -l walltime=240:0:0
    #PBS -l nodes=1:ppn=10
    #PBS -l mem=70GB
    #PBS -e ./log/cellranger-${idx}.stderr.log
    #PBS -o ./log/cellranger-${idx}.stdout.log

    source cellranger-7.1.0.sh
    cd $HOME/workspace/gbm/cellranger

    cellranger count \
        --id=$idx \
        --transcriptome=$HOME/refgenome/spaceranger/refdata-gex-GRCh38-2020-A \
        --fastqs=$HOME/workspace/gbm/Data/SC/${idx} \
        --sample=$idx \
        --localcores=10 \
        --localmem=64

    cellranger mat2csv \
        ${idx}/outs/filtered_feature_bc_matrix \
        ${idx}/outs/filtered_feature_bc_matrix/${idx}.filtered.csv

    cellranger mat2csv \
        ${idx}/outs/raw_feature_bc_matrix \
        ${idx}/outs/raw_feature_bc_matrix/${idx}.raw.csv
EOF

done
