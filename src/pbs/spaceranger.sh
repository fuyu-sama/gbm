cd $HOME/workspace/gbm
declare -A slides=(
    ["21B-603-5"]="V43J31-122"
    ["22F-10823-3"]="V43J31-122"
    ["22F-21576-1"]="V43J31-095"
    ["22F-23738-2"]="V43J31-095"
)

declare -A areas=(
    ["21B-603-5"]="D1"
    ["22F-10823-3"]="A1"
    ["22F-21576-1"]="A1"
    ["22F-23738-2"]="D1"
)

probe_set=$HOME/refgenome/spaceranger/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv

for idx in "${!slides[@]}"; do
    image_file=${slides[$idx]}_${areas[$idx]}_${idx}.tif

    qsub << \
EOF
    #PBS -N $idx
    #PBS -l walltime=240:0:0
    #PBS -l nodes=1:ppn=10
    #PBS -l mem=70GB
    #PBS -e ./log/spaceranger-${idx}.stderr.log
    #PBS -o ./log/spaceranger-${idx}.stdout.log

    source spaceranger-2.1.0.sh
    cd $HOME/workspace/gbm/spaceranger

    spaceranger count \
        --id=$idx \
        --transcriptome=$HOME/refgenome/spaceranger/refdata-gex-GRCh38-2020-A \
        --probe-set=$probe_set \
        --fastqs=$HOME/workspace/gbm/Data/ST/${idx} \
        --sample=$idx \
        --cytaimage=$HOME/workspace/gbm/Data/ST/cytaimage/${image_file} \
        --image=$HOME/workspace/gbm/Data/ST/brightfield/${image_file} \
        --reorient-images=true \
        --slide=${slides[$idx]} \
        --slidefile=$HOME/workspace/gbm/Data/slides/${slides[$idx]}.gpr \
        --area=${areas[$idx]} \
        --localcores=10 \
        --localmem=64

    spaceranger mat2csv \
        ${idx}/outs/filtered_feature_bc_matrix \
        ${idx}/outs/filtered_feature_bc_matrix/${idx}.filtered.csv

    spaceranger mat2csv \
        ${idx}/outs/raw_feature_bc_matrix \
        ${idx}/outs/raw_feature_bc_matrix/${idx}.raw.csv
EOF

done
