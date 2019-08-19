#!/usr/bin/env nextflow

// Get sample info from sample sheet
// Minimum information that are needed in the sample sheet are SampleID, Gender and BAM file location
Channel.fromPath( file(params.sample_sheet) )
        .splitCsv(header: true, sep: '\t')
        .map{row ->
            def sample_id = row['SampleID']
            def gender = row['Gender']
            def bam_file = file(row['BAM'])
            return [ sample_id, gender, bam_file ]
        }.into{samples_1; samples_2; samples_3; samples_4; samples_5}

autosomes = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22".split(',')

process print_sample_info {
    tag { "${sample_id}" }
    echo true
    input:
    set val(sample_id), val(gender), file(bam_file) from samples_1
    script:
    """
    printf "[sample_info] sample: ${sample_id}\tGender: ${gender}\tBAM: ${bam_file}\n"
    """
}

process run_haplotype_caller_on_autosomes {
    tag { "${params.project_name}.${sample_id}.${chr}.rHCoA" }
    memory { 8.GB * task.attempt }
    cpus { 4 }
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false

    input:
    set val(sample_id), val(gender), val(bam_file) from samples_2
	  each chr from autosomes

    output:
	  set val(sample_id), file("${sample_id}.${chr}.vcf.gz")  into autosome_calls
	  set val(sample_id), file("${sample_id}.${chr}.vcf.gz.tbi") into autosome_calls_indexes

    script:
    call_conf = 30 // set default
    if ( params.sample_coverage == "high" )
      call_conf = 30
    else if ( params.sample_coverage == "low" )
      call_conf = 10
    """
    ${params.gatk_base}/gatk --java-options "-Xmx${params.gatk_hc_mem}"  \
    HaplotypeCaller \
    -R ${params.ref_seq} \
    -I $bam_file \
    --dbsnp ${params.dbsnp_sites} \
    --L chr$chr \
    --genotyping-mode DISCOVERY \
    -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
    -stand-call-conf ${call_conf} \
    --sample-ploidy 2 \
    -O ${sample_id}.${chr}.vcf.gz
    """
}

// Now do X and Y calling
samples_3.filter{it[1] == 'M'}.into{samples_male_1; samples_male_2; samples_male_3; samples_male_4; samples_male_5; samples_male_6}
samples_4.filter{it[1] == 'F'}.set{samples_female}

// Males
process run_haplotype_caller_on_x_par1_male {
     tag { "${params.project_name}.${sample_id}.rHCoXP1M" }
     memory { 8.GB * task.attempt }
     cpus { 4 }
     publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false

     input:
     set val(sample_id), val(gender), val(bam_file) from samples_male_1

     output:
	   set val(sample_id), file("${sample_id}.X_PAR1.vcf.gz") into x_par1_calls
	   set val(sample_id), file("${sample_id}.X_PAR1.vcf.gz.tbi") into x_par1_calls_indexes

     script:
     call_conf = 30 // set default
     if ( params.sample_coverage == "high" )
       call_conf = 30
     else if ( params.sample_coverage == "low" )
       call_conf = 10
     """
     ${params.gatk_base}/gatk --java-options "-Xmx${params.gatk_hc_mem}"  \
     HaplotypeCaller \
     -R ${params.ref_seq} \
     -I $bam_file \
     --dbsnp ${params.dbsnp_sites} \
     --L chrX:10001-2781479 \
     --genotyping-mode DISCOVERY \
     -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
     -stand-call-conf ${call_conf} \
     --sample-ploidy 2 \
     -O ${sample_id}.X_PAR1.vcf.gz
     """
}

process run_haplotype_caller_on_x_par2_male {
     tag { "${params.project_name}.${sample_id}.rHCoXP2M" }
     memory { 4.GB * task.attempt }
     cpus { 4 }
     publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false

     input:
     set val(sample_id), val(gender), val(bam_file) from samples_male_2

     output:
	   set val(sample_id), file("${sample_id}.X_PAR2.vcf.gz") into x_par2_calls
	   set val(sample_id), file("${sample_id}.X_PAR2.vcf.gz.tbi") into x_par2_calls_indexes

     script:
     call_conf = 30 // set default
     if ( params.sample_coverage == "high" )
       call_conf = 30
     else if ( params.sample_coverage == "low" )
       call_conf = 10
     """
     ${params.gatk_base}/gatk --java-options "-Xmx${params.gatk_hc_mem}"  \
     HaplotypeCaller \
     -R ${params.ref_seq} \
     -I $bam_file \
     --dbsnp ${params.dbsnp_sites} \
     --L chrX:155701383-156030895 \
     --genotyping-mode DISCOVERY \
     -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
     -stand-call-conf ${call_conf} \
     --sample-ploidy 2 \
     -O ${sample_id}.X_PAR2.vcf.gz
     """
}

process run_haplotype_caller_on_x_nonpar_male {
     tag { "${params.project_name}.${sample_id}.rHCoXNPM" }
     memory { 8.GB * task.attempt }
     cpus { 4 }
     publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false

     input:
     set val(sample_id), val(gender), val(bam_file) from samples_male_3

     output:
	   set val(sample_id), file("${sample_id}.X_nonPAR.vcf.gz") into x_nonpar_calls
	   set val(sample_id), file("${sample_id}.X_nonPAR.vcf.gz.tbi") into x_nonpar_calls_indexes

     script:
     call_conf = 30 // set default
     if ( params.sample_coverage == "high" )
       call_conf = 30
     else if ( params.sample_coverage == "low" )
       call_conf = 10
     """
     ${params.gatk_base}/gatk --java-options "-Xmx${params.gatk_hc_mem}"  \
     HaplotypeCaller \
     -R ${params.ref_seq} \
     -I $bam_file \
     --dbsnp ${params.dbsnp_sites} \
     --L chrX -XL chrX:10001-2781479 -XL chrX:155701383-156030895 \
     --genotyping-mode DISCOVERY \
     -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
     -stand-call-conf ${call_conf} \
     --sample-ploidy 1 \
     -O ${sample_id}.X_nonPAR.vcf.gz
     """
}

process run_haplotype_caller_on_y_par1_male {
     tag { "${params.project_name}.${sample_id}.rHCoYP1M" }
     memory { 8.GB * task.attempt }
     cpus { 4 }
     publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false

     input:
     set val(sample_id), val(gender), val(bam_file) from samples_male_4

     output:
	   set val(sample_id), file("${sample_id}.Y_PAR1.vcf.gz") into y_par1_calls
	   set val(sample_id), file("${sample_id}.Y_PAR1.vcf.gz.tbi") into y_par1_calls_indexes

     script:
     call_conf = 30 // set default
     if ( params.sample_coverage == "high" )
       call_conf = 30
     else if ( params.sample_coverage == "low" )
       call_conf = 10
     """
     ${params.gatk_base}/gatk --java-options "-Xmx${params.gatk_hc_mem}"  \
     HaplotypeCaller \
     -R ${params.ref_seq} \
     -I $bam_file \
     --dbsnp ${params.dbsnp_sites} \
     --L chrY:10001-2781479 \
     --genotyping-mode DISCOVERY \
     -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
     -stand-call-conf ${call_conf} \
     --sample-ploidy 2 \
     -O ${sample_id}.Y_PAR1.vcf.gz
     """
}

process run_haplotype_caller_on_y_par2_male {
     tag { "${params.project_name}.${sample_id}.rHCoYP2M" }
     memory { 8.GB * task.attempt }
     cpus { 4 }
     publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false

     input:
     set val(sample_id), val(gender), val(bam_file) from samples_male_5

     output:
	   set val(sample_id), file("${sample_id}.Y_PAR2.vcf.gz") into y_par2_calls
	   set val(sample_id), file("${sample_id}.Y_PAR2.vcf.gz.tbi") into y_par2_calls_indexes

     script:
     call_conf = 30 // set default
     if ( params.sample_coverage == "high" )
       call_conf = 30
     else if ( params.sample_coverage == "low" )
       call_conf = 10
     """
     ${params.gatk_base}/gatk --java-options "-Xmx${params.gatk_hc_mem}"  \
     HaplotypeCaller \
     -R ${params.ref_seq} \
     -I $bam_file \
     --dbsnp ${params.dbsnp_sites} \
     --L chrY:56887903-57217415 \
     --genotyping-mode DISCOVERY \
     -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
     -stand-call-conf ${call_conf} \
     --sample-ploidy 2 \
     -O ${sample_id}.Y_PAR2.vcf.gz
     """
}

process run_haplotype_caller_on_y_nonpar_male {
     tag { "${params.project_name}.${sample_id}.rHCoYNPM" }
     memory { 8.GB * task.attempt }
     cpus { 4 }
     publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false

     input:
     set val(sample_id), val(gender), val(bam_file) from samples_male_6

     output:
	   set val(sample_id), file("${sample_id}.Y_nonPAR.vcf.gz") into y_nonpar_calls
	   set val(sample_id), file("${sample_id}.Y_nonPAR.vcf.gz.tbi") into y_nonpar_calls_indexes

     script:
     call_conf = 30 // set default
     if ( params.sample_coverage == "high" )
       call_conf = 30
     else if ( params.sample_coverage == "low" )
       call_conf = 10
     """
     ${params.gatk_base}/gatk --java-options "-Xmx${params.gatk_hc_mem}"  \
     HaplotypeCaller \
     -R ${params.ref_seq} \
     -I $bam_file \
     --dbsnp ${params.dbsnp_sites} \
     --L chrY -XL chrY:10001-2781479 -XL chrY:56887903-57217415 \
     --genotyping-mode DISCOVERY \
     -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
     -stand-call-conf ${call_conf} \
     --sample-ploidy 1 \
     -O ${sample_id}.Y_nonPAR.vcf.gz
    """
}

// Females
process run_haplotype_caller_on_x_female {
     tag { "${params.project_name}.${sample_id}.rHCoXF" }
     memory { 8.GB * task.attempt }
     cpus { 4 }
     publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false

     input:
     set val(sample_id), val(gender), val(bam_file) from samples_female

     output:
	   set val(sample_id), file("${sample_id}.X.vcf.gz") into x_calls
	   set val(sample_id), file("${sample_id}.X.vcf.gz.tbi") into x_calls_indexes

     script:
     call_conf = 30 // set default
     if ( params.sample_coverage == "high" )
       call_conf = 30
     else if ( params.sample_coverage == "low" )
       call_conf = 10
     """
     ${params.gatk_base}/gatk --java-options "-Xmx${params.gatk_hc_mem}"  \
     HaplotypeCaller \
     -R ${params.ref_seq} \
     -I $bam_file \
     --dbsnp ${params.dbsnp_sites} \
     --L chrX \
     --genotyping-mode DISCOVERY \
     -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
     -stand-call-conf ${call_conf} \
     --sample-ploidy 2 \
     -O ${sample_id}.X.vcf.gz
     """

}

process run_haplotype_caller_on_mt {
    tag { "${params.project_name}.${sample_id}.rHCoMT" }
    memory { 8.GB * task.attempt }
    cpus { 4 }
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false

    input:
    set val(sample_id), val(gender), file(bam_file) from samples_5

    output:
	  set val(sample_id), file("${sample_id}.MT.vcf.gz") into mt_calls
	  set val(sample_id), file("${sample_id}.MT.vcf.gz.tbi") into mt_calls_indexes

    script:
    call_conf = 30 // set default
    if ( params.sample_coverage == "high" )
      call_conf = 30
    else if ( params.sample_coverage == "low" )
      call_conf = 10
    """
    ${params.gatk_base}/gatk --java-options "-Xmx${params.gatk_hc_mem}"  \
    HaplotypeCaller \
    -R ${params.ref_seq} \
    -I $bam_file \
    --dbsnp ${params.dbsnp_sites} \
    --L chrM \
    --genotyping-mode DISCOVERY \
    -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
    -stand-call-conf ${call_conf} \
    --sample-ploidy 2 \
    -O ${sample_id}.MT.vcf.gz
    """
}

autosome_calls.mix(mt_calls,x_par1_calls,x_nonpar_calls,x_par2_calls,x_calls,y_par1_calls,y_nonpar_calls,y_par2_calls).groupTuple().set{all_calls}

process combine_VCFs {
     tag { "${params.project_name}.${sample_id}.cCVCF" }
     memory { 4.GB * task.attempt }
     cpus { 20 }
     publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false

     input:
     set val(sample_id), file(vcf) from all_calls

     output:
	   set val(sample_id), file("${sample_id}.vcf.gz") into combine_calls
	   set val(sample_id), file("${sample_id}.vcf.gz.tbi") into combine_calls_indexes

     script:
     if (vcf.size() == 29) // working with a male sample
     """
     ${params.gatk_base}/gatk --java-options "-Xmx${params.gatk_sv_mem}"  \
     SortVcf \
     -I ${sample_id}.X_PAR1.vcf.gz \
     -I ${sample_id}.X_PAR2.vcf.gz \
     -I ${sample_id}.X_nonPAR.vcf.gz \
     -O ${sample_id}.X.vcf.gz

     ${params.gatk_base}/gatk --java-options "-Xmx${params.gatk_sv_mem}"  \
     SortVcf \
     -I ${sample_id}.Y_PAR1.vcf.gz \
     -I ${sample_id}.Y_PAR2.vcf.gz \
     -I ${sample_id}.Y_nonPAR.vcf.gz \
     -O ${sample_id}.Y.vcf.gz

     echo "${vcf.join('\n')}" | grep "\\.1\\.vcf.gz" > ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.2\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.3\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.4\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.5\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.6\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.7\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.8\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.9\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.10\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.11\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.12\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.13\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.14\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.15\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.16\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.17\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.18\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.19\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.20\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.21\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.22\\.vcf.gz" >> ${sample_id}.vcf.list
     echo ${sample_id}.X.vcf.gz >> ${sample_id}.vcf.list
     echo ${sample_id}.Y.vcf.gz >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.MT\\.vcf\\.gz" >> ${sample_id}.vcf.list
    
     ${params.gatk_base}/gatk --java-options "-Xmx${params.gatk_gv_mem}"  \
     GatherVcfs \
     -I ${sample_id}.vcf.list \
     -O ${sample_id}.vcf.gz # GatherVCF does not index the VCF. The VCF will be indexed in the next tabix operation.

     ${params.tabix_base}/tabix -p vcf ${sample_id}.vcf.gz 
     """
     else if (gvcf.size() == 24) // working with a female  sample
     """
     echo "${vcf.join('\n')}" | grep "\\.1\\.vcf.gz" > ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.2\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.3\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.4\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.5\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.6\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.7\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.8\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.9\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.10\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.11\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.12\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.13\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.14\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.15\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.16\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.17\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.18\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.19\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.20\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.21\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.22\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.X\\.vcf.gz" >> ${sample_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.MT\\.vcf\\.gz" >> ${sample_id}.vcf.list

     ${params.gatk_base}/gatk --java-options "-Xmx${params.gatk_gv_mem}"  \
     GatherVcfs \
     -I ${sample_id}.vcf.list \
     -O ${sample_id}.vcf.gz # GatherVCF does not index the VCF. The VCF will be indexed in the next tabix operation.
     
     ${params.tabix_base}/tabix -p vcf ${sample_id}.vcf.gz
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
