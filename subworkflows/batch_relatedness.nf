#!/usr/bin/env nextflow

include { BATCH_HOMOZYGOSITY } from '../modules/savvy.nf'
include { RELATEDNESS } from '../modules/vcftools_relatedness.nf'

workflow BATCH_RELATEDNESS {
    take:
    ch_mixed_vcf

    main:
    ch_mixed_vcf
        .collect(flat: false)
        .map { rows ->
            // Sort by ID
            rows.sort { it[0] }               
            tuple( rows.collect{ it[0] },    // List of IDs
                   rows.collect{ it[1] },    // List of Sexes
                   rows.collect{ it[2] },    // List of Families
                   rows.collect{ it[3] },    // List of Sample Counts
                   rows.collect{ it[4] },    // List of VCF Files
                   rows.collect{ it[5] } )   // List of CSI Index Files
        }
        .set { ch_vcf_collapsed } 

    // Conditional evaluation blocks
    if ( !file("${params.batchDir}/r04_metrics/relatedness2.csv").exists() || params.rerun_all ) {
        RELATEDNESS(ch_vcf_collapsed)
    }
    
    if ( !file("${params.batchDir}/r04_metrics/homozygosity.csv").exists() || params.rerun_all ) {
        BATCH_HOMOZYGOSITY(ch_vcf_collapsed)
    }
}