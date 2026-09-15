include { BWA_MEM } from '../modules/bwa_mem.nf'
include { KEVLAR_COUNT } from '../modules/kevlar_count.nf'
include { MERGE_BAMS } from '../modules/merge_sam.nf'
include { FIXMATE } from '../modules/fix_mate.nf'
include { MARKDUP } from '../modules/mark_duplicates.nf'
workflow FASTQ_TO_BAM {
    take:
    ch_parsed

    main:
    BWA_MEM (ch_parsed)
        .set { ch_bwa_out }

    ch_bwa_out
        .filter { row -> row[4] >= 2 }
        .map { row ->
            def (id, sex, family, trio, laneCount, famSampleCount, sam_file) = row
            def key = [id:id, sex:sex, family:family, trio:trio, laneCount:laneCount, famSampleCount:famSampleCount]
            def gKey = groupKey(key, laneCount)          
            tuple(gKey, sam_file)
        }
        .groupTuple() 
        .map { key, sam_file -> tuple(key.getGroupTarget(), sam_file) }
        .set  { ch_multilane_bams }

    MERGE_BAMS (ch_multilane_bams)
        .set  { ch_merged_bams }

    ch_bwa_out
        .filter { row -> row[4] < 2 }
        .mix(ch_merged_bams)
        .set { ch_mix_bams }

    FIXMATE (ch_mix_bams)
        .set  { ch_fixmate }
  
    MARKDUP (ch_fixmate,)

    MARKDUP.out.cram
        .set { ch_final_cram }

    emit:
    ch_final_cram
}