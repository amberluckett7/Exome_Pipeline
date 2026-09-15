#!/usr/bin/env nextflow

include { DENOVOCNN } from '../modules/denovocnn.nf'
include { DENOVO_LARGE_INSERTS } from '../modules/denovo_large_inserts.nf'
include { buildStatusById; parseLineToMeta; parseLineToTuple; parseLineToTupleSpring } from '../lib/helpers.nf'

workflow TRIO_DE_NOVO {
    take:
    ch_final_bam
    ch_single_sample_vcf
    ch_ref_fasta
    ch_ref_fai
    ch_gnomad_common
    ch_gnomad_common_idx
    statusById

    main:

    // Structural grouping
    ch_final_bam
        .filter { row -> row[3] != 'NA' }
        .map { row ->
            def (id, sex, family, trio, famSampleCount, bam_file, bai_file) = row
            def key = [family:family, trio:trio, famSampleCount:famSampleCount]        
            tuple(key, [id:id, sex:sex, bam:bam_file, bai:bai_file])
        }
        .groupTuple(size: 3)
        .map { key, bam_bai -> 
            def meta = key
            def orderedIds = meta.trio.tokenize('-')
            def idToIndex = orderedIds.collectEntries { famId ->
                [(famId): orderedIds.indexOf(famId)]
            }
            def sortedGroup = bam_bai.sort { rec -> idToIndex[rec.id] }
            
            def sortedIds  = sortedGroup*.id
            def sortedSex  = sortedGroup*.sex         
            def sortedbams = sortedGroup*.bam
            def sortedbais = sortedGroup*.bai
            
            // Produces output required to process input validation blocks
            tuple(sortedIds, sortedSex, meta.trio.tokenize('-'), sortedbams, sortedbais)
        }  
        .set { ch_trios_bam }

    // Linking BAM metadata channels with structural variants (VCFs)
    ch_final_bam
        .filter { row -> row[3] != 'NA' }
        .join(ch_single_sample_vcf)
        .map { row -> 
            def (id, sex_bam, family_bam, trio, famCt_bam, bam, bai, sex_vcf, family_vcf, famCt_vcf, vcf, csi) = row

            assert sex_bam == sex_vcf
            assert family_bam == family_vcf
            assert famCt_bam == famCt_vcf

            def key = [trio:trio]        
            tuple(key, [id:id, bam:bam, bai:bai, vcf:vcf, csi:csi])
        }
        .groupTuple(size: 3)
        .map { key, samples -> 
            def meta = key
            def orderedIds = meta.trio.tokenize('-')
            def idToIndex = orderedIds.collectEntries { famId ->
                [(famId): orderedIds.indexOf(famId)]
            }
            def sortedGroup = samples.sort { rec -> idToIndex[rec.id] }
            
            def sortedIds  = sortedGroup*.id        
            def sortedbams = sortedGroup*.bam
            def sortedbais = sortedGroup*.bai
            def sortedvcfs = sortedGroup*.vcf
            def sortedcsis = sortedGroup*.csi
            
            tuple(sortedIds, meta.trio.tokenize('-'), sortedbams, sortedbais, sortedvcfs, sortedcsis)
        }  
        .set { ch_trios_bam_vcf_concat }
    
    // Filtering blocks 
    ch_trios_bam_vcf_concat
        .filter { ids, trio, bam, bai, vcf, csi ->
            def rep = (ids instanceof List ? ids[0] : ids)
            return statusById[rep]?.denovocnn_needed
        }
        .set { ch_trios_bam_vcf_concat_for_cnn }   

    DENOVOCNN(ch_trios_bam_vcf_concat_for_cnn, ch_ref_fasta, ch_ref_fai, ch_gnomad_common, ch_gnomad_common_idx)

    ch_trios_bam
        .filter { ids, sex, trio, bam, bai ->
            def rep = (ids instanceof List ? ids[0] : ids)
            return statusById[rep]?.denovoLI_needed
        }
        .set { ch_trios_bam_for_large_inserts }   

    DENOVO_LARGE_INSERTS(ch_trios_bam_for_large_inserts)
}