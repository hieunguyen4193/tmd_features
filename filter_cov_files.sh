outdir="/media/hieunguyen/GSHD_HN01/raw_data/bismark_cov/TMD_cov";
mkdir -p ${outdir}/filtered_5reads_cov;
mkdir -p ${outdir}/filtered_3reads_cov;
mkdir -p ${outdir}/filtered_10reads_cov;
mkdir -p ${outdir}/filtered_15reads_cov;
mkdir -p ${outdir}/filtered_20reads_cov;
files=$(ls ${outdir}/raw_TMD_cov);

for file in $files;do filename=${file%.cov*} && \
echo -e "working on sample " $filename "\n" && \
# awk -F'\t' '$5 + $6 >= 5' ${outdir}/raw_TMD_cov/${filename}.cov > ${outdir}/filtered_5reads_cov/${filename}.filtered.cov && \
# awk -F'\t' '$5 + $6 >= 3' ${outdir}/raw_TMD_cov/${filename}.cov > ${outdir}/filtered_3reads_cov/${filename}.filtered.cov && \
# awk -F'\t' '$5 + $6 >= 10' ${outdir}/raw_TMD_cov/${filename}.cov > ${outdir}/filtered_10reads_cov/${filename}.filtered.cov;
awk -F'\t' '$5 + $6 >= 15' ${outdir}/raw_TMD_cov/${filename}.cov > ${outdir}/filtered_15reads_cov/${filename}.filtered.cov;
awk -F'\t' '$5 + $6 >= 20' ${outdir}/raw_TMD_cov/${filename}.cov > ${outdir}/filtered_20reads_cov/${filename}.filtered.cov;
done