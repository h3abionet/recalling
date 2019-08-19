#!/usr/bin/env nextflow

// Get sample info from sample sheet
// Minimum information that are needed in the sample sheet are SampleID, Gender and gVCF file location
Channel.fromPath( file(params.sample_sheet) )
        .splitCsv(header: true, sep: '\t')
        .map{row ->
            def sample_id = row['SampleID']
            def gvcf_file = file(row['gVCF'])
            return [ sample_id, gvcf_file ]
        }
        .tap{samples_1; samples_2}
        .map { sample_id, gvcf_file ->
            return [ gvcf_file ]
        }
        .collect().set { gvcf_files }

chroms = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM".split(',')

process print_sample_info {
    tag { "${sample_id}" }
    echo true
    
    input:
    set val(sample_id), file(gvcf_file) from samples_1
    
    script:
    """
    printf "[sample_info] sample: ${sample_id}\tgVCF: ${gvcf_file}\n"
    """
}

process create_variant_list {
    tag { "${params.project_name}.${params.cohort_id}.cVL" }
    publishDir "${params.out_dir}/${params.cohort_id}/combine-gvcfs", mode: 'copy', overwrite: false

    input:
    val gvcf_file from gvcf_files
    
    output:
    file("gvcf.list") into gvcf_list
    
    script:
    """
    echo "${gvcf_file.join('\n')}" > gvcf.list
    """
}

process run_combine_gvcfs {
    tag { "${params.project_name}.${params.cohort_id}.${chr}.rCG" }
    label 'bigmem'
    publishDir "${params.out_dir}/${params.cohort_id}", mode: 'copy', overwrite: false

    input:
    file(gvcf_list)
    each chr from chroms

    output:
    file("${params.cohort_id}.${chr}.g.vcf.gz")  into cohort_chr_calls
    file("${params.cohort_id}.${chr}.g.vcf.gz.tbi") into cohort_chr_indexes

    script:
    mem = task.memory.toGiga() - 3
    """
    ${params.gatk_base}/gatk --java-options "-XX:+UseSerialGC -Xms4g -Xmx${mem}g"  \
    CombineGVCFs \
    -R ${params.ref_seq} \
    --L ${chr} \
    --variant ${gvcf_list} \
    -O ${params.cohort_id}.${chr}.g.vcf.gz
    """
}

cohort_chr_calls.toList().into{ cohort_calls }

process run_concat_combine_gvcf {
     tag { "${params.project_name}.${params.cohort_id}.rCCG" }
     label 'bigmem'
     publishDir "${params.out_dir}/${params.cohort_id}", mode: 'copy', overwrite: false

     input:
     file(gvcf) from cohort_calls

     output:
	   set val(params.cohort_id), file("${params.cohort_id}.g.vcf.gz") into combine_calls
	   set val(params.cohort_id), file("${params.cohort_id}.g.vcf.gz.tbi") into combine_calls_indexes

     script:
     mem = task.memory.toGiga() - 3
     """
     echo "${gvcf.join('\n')}" | grep "\\.chr1\\.g.vcf.gz" > ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.chr2\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.chr3\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.chr4\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.chr5\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.chr6\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.chr7\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.chr8\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.chr9\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.chr10\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.chr11\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.chr12\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.chr13\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.chr14\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.chr15\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.chr16\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.chr17\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.chr18\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.chr19\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.chr20\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.chr21\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.chr22\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.chrX\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.chrY\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.chrM\\.g\\.vcf\\.gz" >> ${params.cohort_id}.gvcf.list
    
     ${params.gatk_base}/gatk --java-options "-XX:+UseSerialGC -Xms4g -Xmx${mem}g"  \
     GatherVcfs \
     -I ${params.cohort_id}.gvcf.list \
     -O ${params.cohort_id}.g.vcf.gz # GatherVCF does not index the VCF. The VCF will be indexed in the next tabix operation.
     ${params.tabix_base}/tabix -p vcf ${params.cohort_id}.g.vcf.gz 
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
