#!/usr/bin/env nextflow

"""
Author: Mamana M.
Affiliation: University of Cape Town

Latest modification: Gerrit divided gene and genome calling into two separate nextflow workflows.
"""

chromosomes = params.chromosomes.split(',')

adme_samples_ch = file("{params.adme_dir}/${params.adme_samples}")

bed_file = Channel.fromPath("${params.gene_region_dir}/${params.gene_set}-*bed").map {
  file ->
    m = file =~ /.*${params.gene_set}-([X0-9]+).bed/  
  return [m[0][1],file]
}.filter { it[1].readLines().size()> 0 }
.view  { it -> it }


Channel.fromFilePairs("${params.gvcf_dir}/*.[X1-9]*.g.vcf*/") 
  { file -> 
           b = file.baseName
           m = b =~ /.*\.([X0-9]+)\.g.*/  
           return m[0][1]
   }.set{ gvcf_file_cha }

bed_file.join(gvcf_file_cha).map { chr, a, b -> [chr, a, b[0], b[1]]} . set { gt_ch }


process run_genotype_gvcf_on_genes {
    tag { "${params.project_name}.${params.cohort_id}.${chr}.rGGoG" }
    label "bigmem"
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
    input:
      set val(chr), file(bed), file(vcf), file (tbi)  from gt_ch
    output:
      set val(chr), file(vcf_out), file(vcf_index_out) into vcf
    script:
         call_conf = 30 // set default
         if ( params.sample_coverage == "high" )
           call_conf = 30
         else if ( params.sample_coverage == "low" )
           call_conf = 10
        vcf_out = "${params.gene_set}-${chr}-genes.vcf.gz"
        vcf_index_out = "${vcf_out}.tbi"
	
        """
        ${params.gatk_base}/gatk \
            GenotypeGVCFs \
            -R ${params.ref_seq} \
            -L ${bed} \
            -V ${vcf}\
            -stand-call-conf ${call_conf} \
            -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
            -O ${vcf_out}
        """
}

process extract_adme {
  input:
  set val(chr), file(vcf), file(tbi) from vcf
  file(adme_samples) from adme_samples_ch
  output:
    set file(data), file("${data}.tbi")
    publishDir "${params.adme_dir}/vcf", mode: 'copy', overwrite: false
  script:
    data = "adme-${params.gene_set}-${chr}.vcf.gz"
    """
      vcftools --gzvcf $vcf --keep ${adme_samples} --recode --stdout | bgzip > ${data}
      tabix $data
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


