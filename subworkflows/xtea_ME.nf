#!/usr/bin/env nextflow

include { XTEA } from '../modules/xtea.nf'

workflow MOBILE_ELEMENTS {
    take:
    bam_files
    ch_ref_fasta
    ch_ref_fai
    ch_ref_gff

    main:

    bam_files
        .map { row ->
            def (id, sex, family, trio, famSampleCount, bam, bai) = row
            def key = [family:family, famSampleCount:famSampleCount]
            def gKey = groupKey(key, famSampleCount)          
            tuple(gKey, [id:id, sex:sex, bam:bam, bai:bai])
        }
        .groupTuple() 
        .map { key, bams -> 
            def meta = key.getGroupTarget()
            def orderedIds = meta.family.tokenize('-')
            def idToIndex = orderedIds.collectEntries { famId ->
                [(famId): orderedIds.indexOf(famId)]
            }
            def sortedGroup = bams.sort { rec ->
                idToIndex[rec.id]
            }
            def sortedIds   = sortedGroup*.id
            def sortedSex   = sortedGroup*.sex        
            def sortedBams  = sortedGroup*.bam
            def sortedBais  = sortedGroup*.bai
            
            tuple(sortedIds, sortedSex, meta.family, meta.famSampleCount, sortedBams, sortedBais)
        }
        .set { family_bams }

    XTEA (family_bams, ch_ref_fasta, ch_ref_fai, ch_ref_gff)
}