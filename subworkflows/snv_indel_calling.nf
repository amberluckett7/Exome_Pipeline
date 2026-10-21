include { DEEP_VARIANT } from '../modules/deep_variant.nf'
include { GLNEXUS } from '../modules/glnexus.nf'
include { ANCESTRY } from '../modules/ancestry.nf'

workflow SNV_INDEL_CALLING {
    take:
    ch_final_bam
    ch_ref_fasta
    ch_ref_fai
    statusById
    metaById
    ch_target_bed

    main:
    DEEP_VARIANT (ch_final_bam, ch_ref_fasta, ch_ref_fai, ch_target_bed)

    DEEP_VARIANT.out.vcf.set { ch_vcf }
    DEEP_VARIANT.out.gvcf.set { ch_gvcf }

    def ch_gvcf_existing = Channel
        .from(metaById.keySet().toList())
        .filter { id -> !statusById[id].snv_needed }
        .map { id ->
            def meta = metaById[id]
            tuple(meta.id, meta.sex, meta.family, meta.famSampleCount, statusById[id].gvcf, statusById[id].gvcf_csi)
        }

    def ch_final_gvcf = ch_gvcf_existing.mix(ch_gvcf)

    ch_final_gvcf
        .filter { row -> row[3] >= 2 }
        .map { row ->
                def (id, sex, family, famSampleCount, gvcf, gvcfcsi) = row
                def key = [family:family, famSampleCount:famSampleCount]
                def gKey = groupKey(key, famSampleCount)          
                tuple(gKey, [id:id, sex:sex, gvcf:gvcf, gvcfcsi:gvcfcsi])
            }
        .groupTuple() 
        .map { key, gvcfs -> 
                def meta = key.getGroupTarget()
                def orderedIds = meta.family.tokenize('-')
                def idToIndex = orderedIds.collectEntries { famId ->
                    [(famId): orderedIds.indexOf(famId)]
                }
                def sortedGroup = gvcfs.sort { rec ->
                    idToIndex[rec.id]
                }
                def sortedIds      = sortedGroup*.id
                def sortedSex      = sortedGroup*.sex        
                def sortedGvcfs    = sortedGroup*.gvcf
                def sortedGvcfCsis = sortedGroup*.gvcfcsi
                tuple(sortedIds, sortedSex, meta.family, meta.famSampleCount, sortedGvcfs, sortedGvcfCsis)
            }
        .set { ch_family_gvcf }

    ch_family_gvcf
    .filter { ids, sex, family, n, gvcf, gvcfcsi ->
        def rep = (ids instanceof List ? ids[0] : ids)
        statusById[rep].cnv_needed
    }
    .set { ch_family_gvcf_needed }  

    GLNEXUS (ch_family_gvcf_needed)
        .set { ch_family_vcf }

    // Maps 6 elements onto 7 by inserting "NA" if analysing a singleton and not a trio
    ch_final_gvcf
        .filter { row -> row[0] in statusById.keySet() && statusById[row[0]].ancestry_needed }
        .map { id, sex, family, count, gvcf, gvcf_csi ->
            tuple(id, sex, family, "NA", count, gvcf, gvcf_csi)
        }
        .set { ch_ancestry } 

    ANCESTRY(ch_ancestry)

    emit:
    single_sample = ch_vcf
    family = ch_family_vcf
}
