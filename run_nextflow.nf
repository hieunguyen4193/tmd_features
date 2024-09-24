params.SampleSheet=""
params.OUTDIR=""
params.bed=""
params.src=""

bedfile=file(params.bed)
src=file(params.src)
Channel
    .fromPath( params.SampleSheet )
    .splitCsv( header:true )
    .map { row -> tuple(row.SampleID,  file(row.TM_BAM))}  
    .view()
    .set { Input_ch }


process Preprocess_bam_file {
    cache "deep"; tag "${sample_id}"
    maxForks 30
    errorStrategy "ignore"
    publishDir "$params.OUTDIR/reads_from_450_regions", mode: 'copy'

    input: 
        tuple sample_id, file(bam_file) from Input_ch
        file bedfile
        file src
    output: 
        tuple sample_id, file("*.csv") into output_ch
    script:
    """
    samtools sort ${bam_file} -o ${sample_id}.sorted.bam
    samtools index ${sample_id}.sorted.bam
    python $src --bam ${sample_id}.sorted.bam --bed ${bedfile} --output .
    """
}

// to run this nextflow pipeline in HPC-R
// nextflow run run_nextflow.nf --SampleSheet /datassd/hieunguyen/tmp/ecd_features/metadata/metadata_cfDNA_lowpdepth_TMD_bam_cov.filtered.csv 
// --OUTDIR output/ --bed /datassd/hieunguyen/tmp/ecd_features/resources/panel.hg19_liftover_to_hg38.bed 
// --src /datassd/hieunguyen/tmp/ecd_features/fetch_reads_from_450_regions.py -resume