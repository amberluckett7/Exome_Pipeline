include { DEPTH_OF_COVERAGE } from '../modules/depth_of_coverage.nf'
include { INDEX_COVERAGE } from '../modules/index_coverage.nf'
include { COVERAGE_REPORT } from '../modules/coverage_report.nf'
include { COVERAGE_BINNER } from '../modules/coverage_binner.nf'
include { XY_COVERAGE } from '../modules/xy_coverage.nf'
include { INSERT_SIZES } from '../modules/insert_sizes.nf'
include { CLIP_RATE } from '../modules/clip_rate.nf'

workflow COVERAGE {
    take:
    ch_final_bam
    ch_control
    ch_target_bed

    main:
    COVERAGE_BINNER(ch_final_bam, ch_target_bed)
    INSERT_SIZES(ch_final_bam, ch_target_bed)
    CLIP_RATE(ch_final_bam)

    DEPTH_OF_COVERAGE(ch_final_bam, ch_target_bed)

    INDEX_COVERAGE(DEPTH_OF_COVERAGE.out.coverage, ch_control, ch_target_bed)

    DEPTH_OF_COVERAGE.out.coverage
        .map { id, sex, cov_file -> cov_file }
        .collect()
        .set { ch_collected_coverage_files }

    COVERAGE_REPORT(ch_collected_coverage_files, ch_target_bed)

    // Collapse id, sex, file into batch arrays for XY_COVERAGE
    COVERAGE_BINNER.out
        .collect(flat: false)
        .map { rows ->
            rows.sort { it[0] }
            tuple( rows.collect{ it[0] },    // ids
                   rows.collect{ it[1] },    // sex
                   rows.collect{ it[2] } )   // coverageBinner files
        }
        .set { ch_cb_collapsed }

    XY_COVERAGE(ch_cb_collapsed)
}
