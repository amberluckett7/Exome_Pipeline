#!/usr/bin/env nextflow

include { EXTRACT_UNMAPPED } from '../modules/extract_unmapped.nf'
include { ALIGN_CYTOMEGALOVIRUS } from '../modules/align_cytomegalovirus.nf'
include { CYTOMEGALOVIRUS_COUNT } from '../modules/cytomegalovirus_count.nf'

workflow CYTOMEGALOVIRUS {
    take:
    ch_final_bam
    ch_check_fastq

    main:

    // Extract unmapped reads from BAM files
    EXTRACT_UNMAPPED(ch_final_bam)
        .set { ch_unmapped }

    // Align those reads, then match up via .join()
    ALIGN_CYTOMEGALOVIRUS(ch_unmapped)
        .join(ch_check_fastq)
        .map { row -> 
            def id         = row[0]
            def sex_sam    = row[1]
            def family_sam = row[2]
            def famCt_sam  = row[3]
            def sam        = row[4]
            
            // Mapped exactly to the array positions
            def sex_cfq    = row[5]
            def family_cfq = row[6]
            def famCt_cfq  = row[7]
            def cfq        = row[8]

            // Workflow checks
            assert sex_sam.toString().toLowerCase() == sex_cfq.toString().toLowerCase() : "Sex mismatch for sample ${id}: ${sex_sam} vs ${sex_cfq}"
            assert family_sam.toString().toLowerCase() == family_cfq.toString().toLowerCase() : "Family mismatch for sample ${id}: ${family_sam} vs ${family_cfq}"
            assert famCt_sam.toString() == famCt_cfq.toString() : "Sample count mismatch for sample ${id}: ${famCt_sam} vs ${famCt_cfq}"
      
            tuple(id, sam, cfq)
        }
        .collect(flat: false)
        .map { rows ->
            // Sort by ID for the downstream java tool
            rows.sort { it[0] }              
            tuple( rows.collect{ it[0] },   
                   rows.collect{ it[1] },   
                   rows.collect{ it[2] } )  
        }
        .set { ch_sam_cfq_collapsed } 

    // Conditionally run the multi-sample counting batch step
    if ( !file("${params.batchDir}/r04_cytomegalovirus/stats").exists() || params.rerun_all ) {
        CYTOMEGALOVIRUS_COUNT(ch_sam_cfq_collapsed)
    }
}
