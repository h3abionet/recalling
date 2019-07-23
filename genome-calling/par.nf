#!/usr/bin/env nextflow

Channel.from( file(params.gvcf_file) )
        .set{ gvcf_file_ch }

ref_seq = Channel.fromPath(params.ref_seq).toList()


gvcf_file_ch = Channel.create()
gvcf_file_2ch = Channel.create()

Channel.fromFilePairs("${params.gvcf_file}")
       { file ->
         b = file.baseName
         m = b =~ /.*\.([0-9]+|X|Y|MT)\.g.*/
         return m[0][1].replaceAll("^0","")     // remove leading 0s
       }.separate ( gvcf_file_ch, gvcf_file_2ch ) { it -> [it,it] }



// NB: we remove leading 0s because sometimes chromsomes are listed 01 not 1 and when we do
// cross we do string not integer matching


// get list of chromosome numbers
chrom_nums = files(params.gvcf_file)
    .findAll {  it =~ /.*.vcf.gz$/ }
    .collect { 
         b = it.baseName
         m = b =~ /.*\.([0-9]+|X|Y|MT)\.g.*/
         return m[0][1].replaceAll("^0","")     // remove leading 0s
     }

println "We are processing the following chromosomes: $chrom_nums"


//  Work out which segments we need for each file
process getSegments {
  input:
    set val(chr), file(gvcf) from gvcf_file_ch
  output:
    stdout into seg_lines
  script:
     vcf = gvcf[0]
     """
     getsegs.py $vcf ${params.seg_size}
     """
}   


// Now we have the segments, match the chrom, vcf, tbi files to  the segments using cross
// -- the final map turns a list of this shape [[chrom, [vcf,tbi]], [chrom,padded_chrom region]] into [chrom, reg-start,reg-end, vcf, tbi]
gvcf_file_2ch.
    cross(seg_lines.splitText().map { it.trim ().split(",") }).map { a,b -> [a[0],b[1],b[2],a[1][0],a[1][1]]}.set { input_ch }


String getFormattedChrom(String chrom) {
       if (chrom == "MT")  // temporarily rename so we put chroms in right order 01 02 .. 22 X Y Z
           return "Z"
       else if (['X','Y'].contains(chrom))
            return chrom;
       else  
	 return chrom.padLeft(2,"0") // we want chromosome order to match lexicographic ordder
}

process run_genotype_gvcf_on_genome {
    tag { "${params.project_name}.${params.cohort_id}.${chrom}:${start}.rGGoG" }
    maxForks 70
    input:
       set val(chrom), val(start), val(end), file (vcf), file (tbi) from input_ch
       file (ref_seq)
    output:
       set val(chrom), file(outf), file("${outf}.tbi") into gg_vcf_set
       file(outf) into gg_vcf
    script:
       region = "${start}-${end}"
       ref_fa = "${ref_seq[0].simpleName}.fasta"
       out_chrom = getFormattedChrom (chrom)
       outf   = "${params.cohort_id}-chr${out_chrom}-${region}.vcf.gz"
       call_conf = 30 // set default
       if ( params.sample_coverage == "high" )  
         call_conf = 30
       else if ( params.sample_coverage == "low" )
         call_conf = 10
       """
       gatk --java-options "-Xmx${task.memory.toGiga()}g" \
        GenotypeGVCFs \
        -R $ref_fa \
        -V ${vcf} \
        -L ${chrom}:${region} \
        -stand-call-conf ${call_conf} \
        -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest \
        -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
        -O $outf
       """
}


gg_vcf.toList().set{ concat_ready  }



process run_concat_vcf {
     tag { "${params.project_name}.${params.cohort_id}.rCV" }
     label "bigmem"
     input:
     file(vcf) from concat_ready

     output:
       file("${params.cohort_id}.vcf.gz") into (combined_calls, combined_calls_1, index_in_ch)
     script:
       vcf_list="${params.cohort_id}.list"
       """
       hostname 
       ls *vcf.gz > $vcf_list
       /usr/bin/time -f "MaxRSS=%M time=%e" gatk --java-options "-Xmx${task.memory.toGiga()}g"  \
           GatherVcfs -I $vcf_list  -O ${params.cohort_id}.vcf.gz 
       """
}

process  indexVCF {
  time '36h'
  memory '8G'
  cpus 2
  tag "indexVCF"
  input:
    file(vcf) from index_in_ch
  output:
    file(tbi) into (index_out_ch, index_out_1_ch)
  script:
    tbi = "${vcf}.tbi"
    """
    bcftools index -t --threads 2 $vcf
    """
}





process run_vqsr_on_snps {
    tag { "${params.project_name}.${params.cohort_id}.rVoS" }
    label "vbigmem"
    input:
      file(vcf) from combined_calls
      file(tbi) from index_out_ch
      file ref_seq 
    output:
       set file("${params.cohort_id}.vcf.recal-SNP.recal"), file("${params.cohort_id}.vcf.recal-SNP.recal.idx"),\
           file("${params.cohort_id}.vcf.recal-SNP.tranches") into snps_vqsr_recal
    script:
    ref  = "${ref_seq[0].simpleName}.fasta"
    """
    hostname
    /usr/bin/time -f "MaxRSS=%M time=%e" gatk --java-options "-Xmx${task.memory.toGiga()}g" \
       VariantRecalibrator    -R $ref \
         -resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${params.hapmap} \
         -resource:omni,known=false,training=true,truth=true,prior=12.0 ${params.omni} \
         -resource:1000G,known=false,training=true,truth=false,prior=10.0 ${params.phase1_snps} \
         -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${params.dbsnp} \
         -an DP -an FS -an SOR -an MQ -an MQRankSum -an QD -an ReadPosRankSum \
         -mode SNP --max-gaussians "${params.max_gaussians_snps}" -V ${vcf} \
         -O "${params.cohort_id}.vcf.recal-SNP.recal" \
         --tranches-file "${params.cohort_id}.vcf.recal-SNP.tranches"
   """
}



process run_apply_vqsr_on_snps {
    tag { "${params.project_name}.${params.cohort_id}.rAVoS" }
    label "bigmem"
    publishDir "${params.out_dir}/${params.cohort_id}/genome-calling", mode: 'copy', overwrite: false
    input:
       file(vcf) from combined_calls_1
       file(tbi) from index_out_1_ch
       set  file(snp_recal), file(snp_recal_index), file(snp_tranches) from snps_vqsr_recal
       file ref_seq
       each chrom from chrom_nums
    output:
       set file(out), file("${out}.tbi") into snps_vqsr_vcf_chrom
    script:
      out_chrom = getFormattedChrom(chrom)
      out="${params.cohort_id}-vqsr-SNP-${out_chrom}.vcf.gz"
      ref = "${ref_seq[0].simpleName}.fasta"
      """
      hostname
     /usr/bin/time -f "MaxRSS=%M time=%e"  gatk --java-options "-Xmx${task.memory.toGiga()}g"	\
	 ApplyVQSR \
	 -R ${ref} \
         -L ${chrom} \
	 --recal-file ${snp_recal} \
	 --tranches-file ${snp_tranches} \
	 -mode SNP \
	 -ts-filter-level "${params.ts_filter_level_snps}" \
	 -V ${vcf} \
	 -O $out
      """
}




// Combine all the VQSR VCF files -- needed for indel file
process join_chrom_vcfs {
     tag { "${params.project_name}.${params.cohort_id}.join" }
     label "medmem"
     input:
     file(vcf) from snps_vqsr_vcf_chrom.flatten().toList()
     output:
        file("${params.cohort_id}.vcf.gz") into (snps_vqsr_vcf, vqsr_index_in_ch)
     script:
       vcf_list="${params.cohort_id}.list"
       """
       hostname 
       ls *vcf.gz > $vcf_list
       /usr/bin/time -f "MaxRSS=%M time=%e" gatk --java-options "-Xmx${task.memory.toGiga()}g"  \
           GatherVcfs -I $vcf_list  -O ${params.cohort_id}.vcf.gz 
       """
}


process  indexVQSRVCF {
  time '36h'
  memory '8G'
  cpus 2
  tag "indexVQSRVCF"
  input:
    file(vcf) from vqsr_index_in_ch
  output:
    file(tbi) into vqsr_index_out_ch
  script:
    tbi = "${vcf}.tbi"
    """
    bcftools index -t --threads 2 $vcf
    """
}




process run_vqsr_on_indels {
    tag { "${params.project_name}.${params.cohort_id}.rVoI" }
    label "vbigmem"
    input:
      file(vcf) from snps_vqsr_vcf
      file(vcf_index) from vqsr_index_out_ch
      file ref_seq
    output:
       set file(vcf), file(vcf_index), file("${params.cohort_id}.recal-SNP.vcf.recal-INDEL.recal"), \
           file("${params.cohort_id}.recal-SNP.vcf.recal-INDEL.recal.idx"), \
           file("${params.cohort_id}.recal-SNP.vcf.recal-INDEL.tranches") into indel_vqsr_recal
    script:
      ref = "${ref_seq[0].simpleName}.fasta"
      """
      hostname
      /usr/bin/time -f "MaxRSS=%M time=%e"     gatk --java-options "-Xmx${task.memory.toGiga()}g" \
        VariantRecalibrator \
	-R ${ref} \
	-resource:mills,known=false,training=true,truth=true,prior=12.0 ${params.golden_indels} \
	-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${params.dbsnp} \
	-an DP -an FS -an SOR -an MQ -an MQRankSum -an QD -an ReadPosRankSum \
	-mode INDEL \
	 --max-gaussians "${params.max_gaussians_indels}" \
	-V ${vcf} \
	-O "${params.cohort_id}.recal-SNP.vcf.recal-INDEL.recal" \
	--tranches-file "${params.cohort_id}.recal-SNP.vcf.recal-INDEL.tranches"
   """
}

process run_apply_vqsr_on_indels {
    tag { "${params.project_name}.${params.cohort_id}.rAVoI" }
    label "bigmem"
    publishDir "${params.out_dir}/${params.cohort_id}/genome-calling", mode: 'copy', overwrite: false
    input:
      set file(vcf), file(vcf_index), file(indel_recal), file(indel_recal_index), file(indel_tranches)\
            from indel_vqsr_recal
      file ref_seq
      each chrom from chrom_nums
    output:
      set file(out), \
          file("${out}.tbi") into indel_vqsr_vcf
    script:
      out_chrom = getFormattedChrom(chrom)
      out="${params.cohort_id}-vqsr-INDEL-${out_chrom}.vcf.gz"
      ref = "${ref_seq[0].simpleName}.fasta"
      """
      hostname
      /usr/bin/time -f "MaxRSS=%M time=%e"     gatk --java-options "-Xmx${task.memory.toGiga()}g" \
      ApplyVQSR \
	 -R $ref \
	 --recal-file ${indel_recal} \
	 -L $chrom \
	 --tranches-file ${indel_tranches} \
	 -mode INDEL \
	 -ts-filter-level "${params.ts_filter_level_indels}" \
	 -V ${vcf} \
	 -O  $out
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






