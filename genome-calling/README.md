# Intro

The Nextflow script calls on genome level (runs GenotypeGVCFs) and does calling of VSQR SNPs and INDELs on the whole genome (chr1 -> 22, X, Y and MT).

Please see `nextflow.conf` for GATK version and references databases used. Path to combined gVCF file is also specified in nextflow config.

## To run

For each dataset
1) Modify your `nextflow.config` to read the `gvcf_file` and specify the output directory e.g. `out_dir = "/global/scratch/gerrit/projects/adme/datasets/NA12878/nextflow-out"`
2) Run the workflow
```
nextflow -log nextflow.log run -w /global5/scratch/gerrit/projects/adme/datasets/NA12878-work -c /home/gerrit/code/recalling/genome-calling/nextflow.config.NA12878 /home/gerrit/code/recalling/genome-calling/main.nf -profile cbio -with-report NA12878.report.html -with-trace -with-timeline NA12878.timeline.html -resume
```


## Output

The output directory will contain
1. A per autosome and X, Y and MT VCF file
1. A recalibrated SNP genome VCF file
1. A recalibrated SNP and INDEL genome VCF file


## Genotyping in parallel

Large parts of the process can be parallelised, speeding up the process significantly. Use the `par.nf` script. Note, the following changes

* There are two new parameters that must be set: `seg_size` specified in millions of bases  (i.e., 20 means 20 million) and `max_forks` the number of genotyping steps that can be done in parallel (NB: although in theory these two should be directly correlated, because chromosome sizes are so unevenly distributed, there is a benefit in having relatively small segment size even if `max_forks` is low. 

* Only the final genotype files are published (the main script also puts intermediate files in the published directories)

