#!/usr/bin/env nextflow

include { VCF_FILTER } from '../modules/vcf_filter.nf'
include { VEP } from '../modules/vep.nf'
include { VEP_PROCESSING } from '../modules/vep_processing.nf'

workflow ANNOTATION {
    take:
    ch_mixed_vcf

    main:
    // Filter
    VCF_FILTER(ch_mixed_vcf)
        .set { ch_filtered_vcf }

    // Generate annotations via Ensembl VEP
    VEP(ch_filtered_vcf)
        .set { ch_anno_vcf }

    // Compress and build file indexes
    VEP_PROCESSING(ch_anno_vcf)
        .set { ch_processed_vcf }

    emit:
    vcf = ch_processed_vcf
}
