include { MANTA } from '../modules/manta.nf'
include { SURVINDEL2 } from '../modules/survindel2.nf'
include { CNVKIT } from '../modules/cnvkit.nf'
include { buildStatusById; parseLineToMeta; parseLineToTuple; parseLineToTupleSpring } from '../lib/helpers.nf'

workflow CNV_CALLING {
    take:
    ch_final_bam
    ch_ref_fasta
    ch_ref_fai
    statusById
    ch_target_bed
    ch_target_bed_tbi

    main:

    ch_final_bam
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
            def sortedGroup = bams.sort { rec -> idToIndex[rec.id] }
            def sortedIds   = sortedGroup*.id
            def sortedSex   = sortedGroup*.sex
            def sortedBams  = sortedGroup*.bam
            def sortedBais  = sortedGroup*.bai
            tuple(sortedIds, sortedSex, meta.family, meta.famSampleCount, sortedBams, sortedBais)
        }
        .set { family_bams }

    family_bams
        .filter { ids, sex, family, n, bam, bai ->
            def rep = (ids instanceof List ? ids[0] : ids)
            statusById[rep].cnv_needed
        }
        .set { ch_manta_bam }

    MANTA (ch_manta_bam, ch_ref_fasta, ch_ref_fai, ch_target_bed, ch_target_bed_tbi)

    ch_final_bam
        .filter { row -> row[0] in statusById.keySet() && statusById[row[0]].survindel_needed }
        .map { row ->
            def (id, sex, family, trio, famSampleCount, bam, bai) = row
            tuple(id, sex, trio, family, famSampleCount, bam, bai)
        }
        .set { ch_survindel_bam }

    SURVINDEL2 (ch_survindel_bam, ch_ref_fasta, ch_ref_fai, ch_target_bed)

    ch_final_bam
        .filter { row -> row[0] in statusById.keySet() && statusById[row[0]].cnv_needed }
        .set { ch_cnvkit_bam }

    CNVKIT (ch_cnvkit_bam, ch_ref_fasta, ch_target_bed)

    emit:
    manta = MANTA.out
}
