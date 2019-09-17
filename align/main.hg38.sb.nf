#!/usr/bin/env nextflow

// Get sample info from sample sheet
// Minimum information that are needed in the sample sheet are SampleID, forward and reverse read file location
Channel.fromPath( file(params.sample_sheet) )
        .splitCsv(header: true, sep: '\t')
        .map{row ->
            def sample_id = row['SampleID']
            def fastq_r1_file = file(row['FastqR1'])
            def fastq_r2_file = file(row['FastqR2'])
            return [ sample_id, fastq_r1_file, fastq_r2_file ]
        }.into{samples_1; samples_2; samples_3; samples_4; samples_5}

process print_sample_info {
    tag { "${sample_id}" }
    echo true
    input:
    set val(sample_id), file(fastq_r1_file), file(fastq_r2_file) from samples_1
    script:
    """
    printf "[sample_info] sample: ${sample_id}\tFastqR1: ${fastq_r1_file}\tFastqR2: ${fastq_r2_file}\n"
    """
}

process run_bwa {
    tag { "${params.project_name}.${sample_id}.rBwa" }
    memory { 48.GB * task.attempt }
    cpus { "${params.bwa_threads}" }
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false
    label 'bwa_samtools'

    input:
    set val(sample_id), val(fastq_r1_file), val(fastq_r2_file) from samples_2

    output:
    set val(sample_id), file("${sample_id}.sam")  into raw_sam

    script:
    nr_threads = task.cpus - 1
    readgroup_info="@RG\\tID:$sample_id.0\\tLB:LIBA\\tSM:$sample_id\\tPL:Illumina"
    """
    ${params.bwa_base}/bwa mem \
    -R \"${readgroup_info}\" \
    -t ${nr_threads}  \
    -M \
    ${params.ref_seq} \
    ${fastq_r1_file} \
    ${fastq_r2_file} > ${sample_id}.sam
    """
}

process run_mark_duplicates {
    tag { "${params.project_name}.${sample_id}.rMD" }
    memory { 32.GB * task.attempt }
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false
    label 'samblaster'

    input:
    set val(sample_id), file(sam_file) from raw_sam

    output:
    set val(sample_id), file("${sample_id}.md.bam") into md_bam
   
    script:
    """
    samblaster -M \
    -i ${sam_file} \
    -o ${sample_id}.md.bam
    """
}

process run_samtools_sort {
    tag { "${params.project_name}.${sample_id}.sS" }
    memory { 32.GB * task.attempt }
    cpus { "${params.bwa_threads}" }
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false
    label 'bwa_samtools'

    input:
    set val(sample_id), file(bam_file) from md_bam

    output:
    set val(sample_id), file("${sample_id}.md.sorted.bam") into md_sorted_bam

    script:
      mem = task.memory.toGiga() - 2
      nr_threads = task.cpus - 1 
    """
    samtools \
    sort \
    --threads ${params.bwa_threads} \
    ${bam_file} > ${sample_id}.md.sorted.bam \
    """
}

process run_samtools_index {
    tag { "${params.project_name}.${sample_id}.sI" }
    memory { 4.GB * task.attempt }
    cpus { "${params.bwa_threads}" }
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false
    label 'bwa_samtools'

    input:
    set val(sample_id), file(bam_file) from md_sorted_bam

    output:
    set val(sample_id), file(bam_file), file("${sample_id}.md.sorted.bai") into md_sorted_indexed_bam

    script:
      mem = task.memory.toGiga() - 2
      nr_threads = task.cpus - 1 
    """
    samtools \
    index \
    -m memory_per_thread \
    -@ ${params.bwa_threads} \
    ${bam_file} ${sample_id}.md.sorted.bai \
    """
}


process run_create_recalibration_table {
    tag { "${params.project_name}.${sample_id}.rCRT" }
    memory { 96.GB * task.attempt }
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false

    input:
    set val(sample_id), file(bam_file), file(bam_file_index) from md_sorted_indexed_bam

    output:
    set val(sample_id), file(bam_file), file("${sample_id}.md.sorted.recal.table")  into md_sorted_indexed_recal_table
    
    script:
      mem = task.memory.toGiga() - 2
    """
    ${params.gatk_base}/gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
    BaseRecalibrator \
    --input ${bam_file} \
    --output ${sample_id}.md.sorted.recal.table \
    --TMP_DIR ${params.gatk_tmp_dir} \
    -R ${params.ref_seq} \
    --known-sites ${params.dbsnp} \
    --known-sites ${params.known_indels_1} \
    --known-sites ${params.known_indels_2}
    """
}

process run_recalibrate_bam {
    tag { "${params.project_name}.${sample_id}.rRB" }
    memory { 96.GB * task.attempt }
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false

    input:
    set val(sample_id), file(bam_file), file(recal_table_file) from md_sorted_indexed_recal_table

    output:
    set val(sample_id), file("${sample_id}.md.sorted.recal.bam"), file("${sample_id}.md.sorted.recal.bai") into md_sorted_indexed_recal_bam
    
    script:
      mem = task.memory.toGiga() - 2
    """
    ${params.gatk_base}/gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
     ApplyBQSR \
    --input ${bam_file} \
    --output ${sample_id}.md.sorted.recal.bam \
    --TMP_DIR ${params.gatk_tmp_dir} \
    -R ${params.ref_seq} \
    --bqsr-recal-file ${recal_table_file}
    """
}


process run_samtools_flagstat {
    tag { "${params.project_name}.${sample_id}.rSF" }
    memory { 4.GB * task.attempt } 
    cpus { "${params.bwa_threads}" }
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false
    label 'bwa_samtools'

    input:
    set val(sample_id), file(bam_file), file(bam_file_index) from md_sorted_indexed_recal_bam

    output:
    set val(sample_id), file("${sample_id}.md.sorted.recal.flagstat")  into recal_stats

    """
    samtools flagstat  \
    --threads ${params.bwa_threads} \
    ${bam_file} > ${sample_id}.md.sorted.recal.flagstat  \
    """
}


workflow.onComplete {

    println ( workflow.success ? """
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """ : """
        Failed: ${workflow.errorReport}
        exit status : ${workflow.exitStatus}
        """
    )
}
