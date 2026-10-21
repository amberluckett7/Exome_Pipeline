process BCFTOOLS_FILTER_SV {
    tag "${family}"
    cpus 2
    container 'docker://quay.io/biocontainers/bcftools:1.17--haef29d1_0'

    input:
      tuple val(id), val(sex), val(family), val(famSampleCount), file(vcf), file(tbi)

    output:
      tuple val(id), val(sex), val(family), val(famSampleCount), file("${family}.filtered.SV.vcf.gz"), file("${family}.filtered.SV.vcf.gz.tbi")

    script:
      """
      bcftools view -i 'INFO/SVTYPE!="BND"' -Oz -o ${family}.filtered.SV.vcf.gz ${vcf}
      bcftools index -t ${family}.filtered.SV.vcf.gz
      """
}

process VEP_SV {
    tag "${family}"
    cpus 4
    container 'docker://ensemblorg/ensembl-vep:release_114.1'
    containerOptions "-B ${params.vepData}:/data -B ${params.batchDir} -B ${params.rsync}"

    input:
      tuple val(id), val(sex), val(family), val(famSampleCount), file(vcf), file(tbi)

    output:
      tuple val(id), val(sex), val(family), val(famSampleCount), file("${family}_vep_annotated.SV.vcf.gz"), file("${family}_vep_annotated.SV.vcf.gz.tbi"), emit: sv_vcf

    script:
      """
      vep --dir /data \\
        --fork ${task.cpus} \\
        --cache \\
        --offline \\
        --assembly GRCh38 \\
        --synonyms /data/homo_sapiens/114_GRCh38/chr_synonyms.txt \\
        --regulatory \\
        --hgvsg \\
        --hgvs \\
        --mane \\
        --fasta /data/homo_sapiens/114_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \\
        --format vcf \\
        --vcf \\
        --force_overwrite \\
        --input_file ${vcf} \\
        --output_file ${family}_vep_annotated.SV.vcf || true

      [[ -s "${family}_vep_annotated.SV.vcf" ]] || { echo "VEP produced no output" >&2; exit 1; }

      bgzip -f ${family}_vep_annotated.SV.vcf
      tabix -p vcf ${family}_vep_annotated.SV.vcf.gz

      DEST_DIR="${params.batchDir}/r04_sv_vep"
      mkdir -p "\$DEST_DIR"

      for f in ${family}_vep_annotated.SV.vcf.gz ${family}_vep_annotated.SV.vcf.gz.tbi; do
          TMP_FILE=".\${f}.partial.\$\$"
          attempts=5
          delay=10
          published=0

          for i in \$(seq 1 \$attempts); do
              if ${params.rsync} -a --no-g --checksum --delay-updates --partial \\
                  "\$f" "\$DEST_DIR/\$TMP_FILE"; then
                  mv -f "\$DEST_DIR/\$TMP_FILE" "\$DEST_DIR/\$f"
                  published=1
                  break
              fi
              echo "rsync failed (attempt \$i/\$attempts) — retrying in \${delay}s" >&2
              sleep "\$delay"
          done

          if [[ \$published -eq 0 ]]; then
              echo "Publish failed for \$f" >&2
              exit 1
          fi
      done

      exit 0
      """
}
