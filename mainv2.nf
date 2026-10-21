#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.batchDir = params.batchDir ?: '.'
params.control  = params.control  ?: "${params.batchDir}/control"
params.target_bed = "" // EXOME: Parameter pointing to the BED file for exome regions

include { CONTROL_PARSER } from './modules/parse_control.nf'
include { CONTAM_SMALL } from './modules/clean_call_contamination.nf'
include { EXPANSION_HUNTER_DE_NOVO } from './modules/expansionHunterDeNovo.nf'
include { MOBILE_ELEMENTS } from './subworkflows/xtea_ME.nf'
include { TRIO_DE_NOVO } from './subworkflows/trio_de_novo.nf'
include { FASTQ_TO_BAM } from './subworkflows/fastqtobam.nf'
include { COVERAGE } from './subworkflows/coverage.nf'
include { MITO_ANALYSIS } from './subworkflows/mito_analysis.nf'
include { BCFTOOLS_FILTER_SV; VEP_SV } from './modules/vep_sv.nf'
include { SNV_INDEL_CALLING } from './subworkflows/snv_indel_calling.nf'
include { HOMOZYGOSITY_AND_HAPLOTYPES } from './subworkflows/homozygosity_and_haplotypes.nf'
// include { HLA_REGION } from './subworkflows/hla_region.nf'
include { QC } from './subworkflows/qc.nf'
include { CNV_CALLING } from './subworkflows/cnv_calling.nf'
include { ANNOTATION } from './subworkflows/annotation.nf'
include { CYTOMEGALOVIRUS } from './subworkflows/cytomegalovirus.nf'
include { BATCH_RELATEDNESS } from './subworkflows/batch_relatedness.nf'
include { CRAM2BAM } from './modules/cram_to_bam.nf'
include { SPRING2FQ } from './modules/spring_to_fastq.nf'
include { DUPMETRICS } from './modules/duplicates_metrics.nf'
include { GENE_PANEL_REPORT } from './modules/gene_panel_report.nf'
include { buildStatusById; parseLineToMeta; parseLineToTuple; parseLineToTupleSpring; makeParsedLines } from './lib/helpers.nf'

workflow {

      if (!params.target_bed) {
          error "CRITICAL ERROR: You must provide an exome target regions footprint file using --target_bed <file.bed>"
      }

      spring_dir = file(params.batchDir + "/spring")

      ch_ref_fasta = file(params.referenceFasta)
      ch_ref_fai = file(params.referenceFasta + ".fai")
      ch_ref_dict = file(params.referenceFasta.replaceAll('.fasta$', '.dict'))
      ch_ref_gff = file( params.referenceGFF )
      ch_gnomad_common = file ( params.gnomadCommon )
      ch_gnomad_common_idx = file ( params.gnomadCommon + ".csi" )
      ch_target_bed = Channel.value(file(params.target_bed))
      ch_target_bed_tbi = Channel.value(file("${params.target_bed}.tbi"))

      def ch_control = file(params.control)
      def controlFile = new File(ch_control.toString())
      def controlLines = makeParsedLines(controlFile)

      def metaById = [:]
      controlLines.each { line ->
        def m = parseLineToMeta(line)
        if( !metaById.containsKey(m.id as String) )
            metaById[m.id as String] = m
      }

      def rowsForStatus = controlLines.collect { line ->
        def m = parseLineToMeta(line)
        [ m.id, null, null, m.family, m.trio, m.bamFile ]
      }
      def statusById = buildStatusById(rowsForStatus)

      // ---- DEBUG: print per-sample flags ----
      println "\n=== STATUS FLAGS (per sample) ==="
      metaById.keySet().sort().each { id ->
          def s = statusById[id]
          println sprintf(
              "%-20s cram_needed=%-5s qc_needed=%-5s bam_needed=%-5s bam_direct=%-5s snv_needed=%-5s cov_needed=%-5s cnv_needed=%-5s sv_vep_needed=%-5s",
              id,
              s.cram_needed, s.qc_needed, s.bam_needed, s.bam_direct,
              s.snv_needed, s.cov_needed, s.cnv_needed, s.sv_vep_needed
          )
      }
      println "=== END STATUS FLAGS ===\n"

      // BAM: channel for samples supplied as pre-existing BAMs — inject directly into ch_final_bam
      def ch_bam_direct = Channel
        .from(metaById.keySet().toList())
        .filter { id -> statusById[id].bam_direct }
        .map { id ->
            def meta = metaById[id]
            tuple(meta.id, meta.sex, meta.family, meta.trio, meta.famSampleCount,
                  statusById[id].bam_direct_file, statusById[id].bai_direct_file)
        }

      def ch_cram_existing = Channel
        .from(metaById.keySet().toList())
        .filter { id -> !statusById[id].cram_needed && !statusById[id].bam_direct }
        .map { id ->
            def meta = metaById[id]
            tuple(meta.id,meta.sex,meta.family,meta.trio,meta.famSampleCount, statusById[id].cram, statusById[id].crai)
        }

      def needFQIds = metaById.keySet().findAll { id ->
          (statusById[id].cram_needed || statusById[id].qc_needed) && !statusById[id].bam_direct
      } as Set

      if (spring_dir.isDirectory() && !params.forcefq) {
        def springRows = controlLines
            .collect { parseLineToTupleSpring(it) }
            .findAll { row -> 
                def id = row[0] as String
                needFQIds.contains(id) && !statusById[id].bam_direct 
            }

            Channel.from(springRows).set { ch_spring_need }

            SPRING2FQ (ch_spring_need)
            
            SPRING2FQ.out
                .filter { row ->
                    def id = (row[0] instanceof Map) ? row[0].id : row[0] as String
                    !statusById[id]?.bam_direct
                }
                .set { ch_fq }

      } else {
            Channel.from(controlLines)
            .map { parseLineToTuple(it) }
            .filter { row -> 
                def id = row[0] as String
                needFQIds.contains(id) && !statusById[id].bam_direct 
            }
            .set { ch_fq }
      }

      def ch_fq_existing = Channel
        .from(metaById.keySet().toList())
        .filter { id -> !statusById[id].qc_needed && !statusById[id].bam_direct }
        .map { id ->
            def meta = metaById[id]
            tuple(meta.id,meta.sex,meta.family,meta.famSampleCount, statusById[id].fq_qc)
        }

      ch_fq
        .filter { row -> 
            def id = (row[0] instanceof Map) ? row[0].id : row[0] as String
            id in statusById.keySet() && statusById[id].qc_needed && !statusById[id].bam_direct
        }
        .unique { row -> (row[0] instanceof Map) ? row[0].id : row[0] }
        .set { ch_qc_fq }

      QC(ch_qc_fq)
      QC.out.set { ch_check_fastq }
      
      def ch_final_fq_qc = ch_fq_existing.mix(ch_check_fastq)

      ch_fq
        .filter { row -> 
            def id = (row[0] instanceof Map) ? row[0].id : row[0] as String
            id in statusById.keySet() && statusById[id].cram_needed && !statusById[id].bam_direct
        }
        .set { ch_cram_fq }

      FASTQ_TO_BAM (ch_cram_fq)
      FASTQ_TO_BAM.out.set { ch_cram_new }
      
      def ch_final_cram = ch_cram_existing.mix(ch_cram_new)

      ch_final_cram
        .filter { row -> row[0] in statusById.keySet() && statusById[row[0]].bam_needed }
        .set { ch_cram_to_bam }

      CRAM2BAM (ch_cram_to_bam)
      CRAM2BAM.out.set { ch_cram_converted_bam }

      // BAM: merge direct BAM input with CRAM-converted BAMs
      def ch_final_bam = ch_cram_converted_bam.mix(ch_bam_direct)

      ch_final_bam
        .filter { row -> row[0] in statusById.keySet() && statusById[row[0]].dup_metrics_needed }
        .set { ch_dup_metrics }

      DUPMETRICS (ch_dup_metrics)

      CNV_CALLING (ch_final_bam, ch_ref_fasta, ch_ref_fai, statusById, ch_target_bed, ch_target_bed_tbi)

      // VEP annotation of Manta SV calls
      def ch_manta_existing = Channel
        .from(metaById.keySet().toList())
        .filter { id -> !statusById[id].cnv_needed && statusById[id].sv_vep_needed }
        .map { id ->
            def meta = metaById[id]
            def bDir = file(params.batchDir).toAbsolutePath().normalize().toString()
            def mantaVcf = file("${bDir}/r04_manta/${meta.family}.SV.vcf.gz", checkIfExists: true)
            def mantaTbi = file("${bDir}/r04_manta/${meta.family}.SV.vcf.gz.tbi", checkIfExists: true)
            tuple(meta.id, meta.sex, meta.family, meta.famSampleCount, mantaVcf, mantaTbi)
        }

      // 2. Newly called Manta output from this execution
      CNV_CALLING.out.manta
        .filter { id, sex, family, n, vcf, tbi ->
            def rep = (id instanceof List ? id[0] : id)
            statusById[rep].sv_vep_needed
        }
        .set { ch_manta_new }

      // 3. Mix both existing and new channels into ch_sv_vep
      ch_manta_new
        .mix(ch_manta_existing)
        .set { ch_sv_vep }

      BCFTOOLS_FILTER_SV(ch_sv_vep)
      VEP_SV(BCFTOOLS_FILTER_SV.out)

      ch_final_bam
        .filter { row -> row[0] in statusById.keySet() && statusById[row[0]].cymegv_needed }
        .set { ch_cymegv_bam }

      ch_final_fq_qc
        .filter { row -> row[0] in statusById.keySet() && statusById[row[0]].cymegv_needed }
        .set { ch_cymegv_fq_qc }

      CYTOMEGALOVIRUS (ch_cymegv_bam, ch_cymegv_fq_qc)

      ch_final_bam
        .filter { row -> row[0] in statusById.keySet() && statusById[row[0]].te_needed }
        .set { ch_te_bam }

      MOBILE_ELEMENTS (ch_te_bam, ch_ref_fasta, ch_ref_fai, ch_ref_gff)

      ch_final_bam
        .filter { row -> row[0] in statusById.keySet() && statusById[row[0]].exphunter_needed }
        .set { ch_exphunter_bam }

      EXPANSION_HUNTER_DE_NOVO (ch_exphunter_bam)

      ch_final_bam
        .filter { row -> row[0] in statusById.keySet() && statusById[row[0]].contam_needed }
        .set { ch_contam_bam }

      CONTAM_SMALL (ch_contam_bam)

      ch_final_bam
        .filter { row -> row[0] in statusById.keySet() && statusById[row[0]].cov_needed }
        .set { ch_cov_bam }

      COVERAGE (ch_cov_bam, ch_control, ch_target_bed)

      ch_final_bam
        .filter { row -> row[0] in statusById.keySet() && statusById[row[0]].mito_needed }
        .set { ch_mito_bam }

      MITO_ANALYSIS(ch_mito_bam, ch_ref_fasta, ch_ref_fai, ch_ref_dict)

      ch_final_bam
        .filter { row -> row[0] in statusById.keySet() && statusById[row[0]].snv_needed }
        .set { ch_snv_bam }

      SNV_INDEL_CALLING(ch_snv_bam, ch_ref_fasta, ch_ref_fai, statusById, metaById, ch_target_bed)
      
      SNV_INDEL_CALLING.out.single_sample.filter { row -> row[3] < 2 }
        .map { row -> tuple( [row[0]], [row[1]], row[2], row[3], row[4], row[5] )}
        .set { ch_singleton_vcf }

      SNV_INDEL_CALLING.out.family.mix(ch_singleton_vcf)
        .set { ch_combined_new_vcf }

      def families = metaById.values()
        .collect { it.family as String }
        .unique()

      def ch_vcf_existing = Channel
        .from(families)
        .filter { fam ->
            def rep = fam.tokenize('-')[0]
            !statusById[rep].fam_vcf_needed
        }
        .map { fam ->
            def ids = fam.tokenize('-')
            def rep = ids[0]
            def sexList = ids.collect { sid -> metaById[sid]?.sex }
            tuple(
                ids,
                sexList,
                fam,
                ids.size(),
                statusById[rep].fam_vcf,
                statusById[rep].fam_vcf_csi
            )
        }

      def ch_combined_vcf = ch_vcf_existing.mix(ch_combined_new_vcf)

      SNV_INDEL_CALLING.out.single_sample
        .set { ch_single_sample_vcf }

      BATCH_RELATEDNESS (ch_combined_vcf)

      HOMOZYGOSITY_AND_HAPLOTYPES(ch_combined_vcf,statusById)

      ch_combined_vcf
      .filter { ids, sexList, fam, n, vcf, csi ->
        def rep = (ids instanceof List ? ids[0] : ids)
        statusById[rep].vep_needed
      }
      .set { ch_vep_vcf }

      ANNOTATION(ch_vep_vcf)

      // --- Gene panel report ---
      // SNV and SV VEP outputs may each independently be "fresh this run"
      // or "already existed on disk" — a family can need only one of the
      // two to rerun. Resolve each side separately (live channel output if
      // it ran this time, otherwise from disk if it already existed),
      // THEN join them — rather than requiring both sides to come from the
      // same source, which silently dropped families needing only one side
      // to rerun.
      ch_gene_panel_tsv    = Channel.value(file(params.genePanelTsv))
      ch_gene_panel_script = Channel.value(file(params.genePanelScript))

      def bDir = file(params.batchDir).toAbsolutePath().normalize().toString()

      // --- SNV side: live output from ANNOTATION this run ---
      ANNOTATION.out.vcf
          .map { ids, sexList, family, n, vcf, csi ->
              tuple(family, ids, sexList, n, vcf, csi)
          }
          .set { ch_snv_new }

      // --- SNV side: fallback to disk for families that didn't need to rerun ---
      def ch_snv_existing = Channel
          .from(families)
          .filter { fam ->
              def rep = fam.tokenize('-')[0]
              !statusById[rep].vep_needed
          }
          .map { fam ->
              def ids = fam.tokenize('-')
              def sexList = ids.collect { sid -> metaById[sid]?.sex }
              def snvVcf = file("${bDir}/r04_vep/${fam}_vep_annotated.vcf.gz")
              def snvIdx = file("${bDir}/r04_vep/${fam}_vep_annotated.vcf.gz.csi")
              tuple(fam, ids, sexList, ids.size(), snvVcf, snvIdx)
          }

      def ch_snv_keyed = ch_snv_new.mix(ch_snv_existing)

      // --- SV side: live output from VEP_SV this run ---
      VEP_SV.out.sv_vcf
          .map { id, sex, family, n, vcf, tbi ->
              tuple(family, vcf, tbi)
          }
          .set { ch_sv_new }

      // --- SV side: fallback to disk for families that didn't need to rerun ---
      def ch_sv_existing = Channel
          .from(families)
          .filter { fam ->
              def rep = fam.tokenize('-')[0]
              !statusById[rep].sv_vep_needed
          }
          .map { fam ->
              def svVcf = file("${bDir}/r04_sv_vep/${fam}_vep_annotated.SV.vcf.gz")
              def svTbi = file("${bDir}/r04_sv_vep/${fam}_vep_annotated.SV.vcf.gz.tbi")
              tuple(fam, svVcf, svTbi)
          }

      def ch_sv_keyed = ch_sv_new.mix(ch_sv_existing)

      // Join SNV and SV sides by family — each side independently resolved
      // above, so a family needing only one side to rerun is still covered.
      // Filtered to only families whose gene panel report is actually
      // missing/stale.
      ch_snv_keyed
          .join(ch_sv_keyed)
          .filter { family, ids, sexList, n, snvVcf, snvIdx, svVcf, svTbi ->
              def rep = ids instanceof List ? ids[0] : ids
              statusById[rep].gene_panel_needed
          }
          .map { family, ids, sexList, n, snvVcf, snvIdx, svVcf, svTbi ->
              tuple(
                  tuple(ids, sexList, family, n, snvVcf, snvIdx),
                  tuple(ids, sexList, family, n, svVcf,  svTbi)
              )
          }
          .set { ch_report_input }

      GENE_PANEL_REPORT(
          ch_report_input.map { snv, sv -> snv },
          ch_report_input.map { snv, sv -> sv },
          ch_gene_panel_tsv,
          ch_gene_panel_script
      )

      TRIO_DE_NOVO (ch_final_bam, ch_single_sample_vcf, ch_ref_fasta, ch_ref_fai,ch_gnomad_common,ch_gnomad_common_idx,statusById)
}

workflow.onComplete { wf ->
    if (params.send_mail) {
        try {
            def email = "${System.getenv('USER') ?: 'unknown'}@exeter.ac.uk"
            def msg = """\
Workflow COMPLETED
Project : ${workflow.projectDir}
Run name: ${workflow.runName}
Time    : ${new Date()}
""".stripIndent()

            Mailer.send(email, "Nextflow workflow COMPLETED: ${workflow.runName}", msg)
        } catch (Exception e) {
            log.warn "Mail send failed: ${e.message}"
        }
    }
}
