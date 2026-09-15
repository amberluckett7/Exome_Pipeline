include { DETECT_MITO_CONTIG } from '../modules/detect_mito_contig.nf'
include { MITO_EXTRACT       } from '../modules/mito_extract.nf'
include { MITO_COVERAGE      } from '../modules/mito_coverage.nf'
include { MITO_VARIANTS      } from '../modules/mito_variants.nf'
include { MITO_SALT          } from '../modules/mito_salt.nf'

workflow MITO_ANALYSIS {
    take:
    ch_final_bam
    ch_ref_fasta
    ch_ref_fai
    ch_ref_dict

    main:
    // Detect the mito contig name (chrM / chrMT / MT / M) once per sample,
    // from the header of the final BAM
    DETECT_MITO_CONTIG(ch_final_bam)

    // Join detected contig name onto the final BAM channel using sample id
    ch_final_bam_with_contig = ch_final_bam.join(DETECT_MITO_CONTIG.out.contig)

    // Extract mitochondrial reads using the detected contig
    MITO_EXTRACT(ch_final_bam_with_contig)

    // Coverage and variant calling run on the extracted mito BAM, which
    // carries the contig name through in its output tuple
    MITO_COVERAGE(MITO_EXTRACT.out)
    MITO_VARIANTS(MITO_EXTRACT.out, ch_ref_fasta, ch_ref_fai, ch_ref_dict)

    // MitoSALT runs on the original BAM (it handles its own MT read
    // extraction) but still needs the correct contig name
    MITO_SALT(ch_final_bam_with_contig)
}
