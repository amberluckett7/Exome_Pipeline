include { CHECK_FASTQ } from '../modules/check_fastq.nf'

workflow QC {
    take:
    ch_parsed

    main:
    ch_parsed
        .map { row ->
            def (id, platform, sex, family, trio, flowcell, laneCount, famSampleCount, fastq1, fastq2) = row
            def key = [id:id, sex:sex, family:family, trio:trio, famSampleCount:famSampleCount]
            def fastqFiles = [f1:fastq1,f2:fastq2]
            def gKey = groupKey(key, laneCount)          
            tuple(gKey, fastqFiles)
        }
        .groupTuple() 
        .map { key, fastqFiles -> tuple(key.getGroupTarget(), fastqFiles.f1,fastqFiles.f2) }
        .set  { ch_grouped_fq }

    CHECK_FASTQ (ch_grouped_fq)
        .set { ch_check_fastq }

    emit:
    ch_check_fastq
}