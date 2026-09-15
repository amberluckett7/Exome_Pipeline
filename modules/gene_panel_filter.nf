process GENE_PANEL_FILTER {
    tag { id instanceof List ? id.join("-") : id }
    cpus 2
    container 'docker://quay.io/biocontainers/bcftools:1.17--haef29d1_0'
    containerOptions "-B ${params.batchDir} -B ${params.resourcesDir} -B ${params.rsync}"

    input:
      tuple val(id), val(sex), val(family), val(famSampleCount), file(vcf), file(tbi)
      path gene_panel_bed

    output:
      tuple val(id), file("${prefix}.gene_panel.tsv.gz"), optional: true

    script:
    prefix = id instanceof List ? id.join("-") : id
    """
    bcftools view -R ${gene_panel_bed} -Ou ${vcf} | \\
    bcftools +split-vep -d -A tab \\
      -f "%FILTER\\t%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT][\\t%GQ][\\t%DP]\\t%CSQ\\n" \\
      -o - | \\
    awk -F"\\t" -v OFS="\\t" -v sid="${prefix}" -v mindp="${params.genePanelMinDP}" -v mingq="${params.genePanelMinGQ}" '
      {
        filt=\$1; gq=\$7; dp=\$8
        if (filt != "PASS" && filt != ".") next
        if (gq == "." || gq == "") next
        if (dp == "." || dp == "") next
        if (gq+0 < mingq) next
        if (dp+0 < mindp) next
        print sid, \$0
      }
    ' > ${prefix}.gene_panel.raw.tsv

    {
      echo -e "sample_id\\tFILTER\\tCHROM\\tPOS\\tREF\\tALT\\tGT\\tGQ\\tDP\\tCSQ"
      cat ${prefix}.gene_panel.raw.tsv
    } | bgzip -c > ${prefix}.gene_panel.tsv.gz

    DEST_DIR="${params.batchDir}/gene_panel"
    mkdir -p "\$DEST_DIR"

    if [[ -s "${prefix}.gene_panel.tsv.gz" ]]; then
        TMP_FILE=".${prefix}.gene_panel.tsv.gz.partial.\$\$"
        attempts=5
        delay=10
        published=0

        for i in \$(seq 1 \$attempts); do
            if ${params.rsync} -a --no-g --checksum --delay-updates --partial \\
                "${prefix}.gene_panel.tsv.gz" "\$DEST_DIR/\$TMP_FILE"; then
                mv -f "\$DEST_DIR/\$TMP_FILE" "\$DEST_DIR/${prefix}.gene_panel.tsv.gz"
                published=1
                break
            fi
            echo "rsync failed (attempt \$i/\$attempts) — retrying in \${delay}s" >&2
            sleep "\$delay"
        done

        if [[ \$published -eq 0 ]]; then
            echo "Publish failed for ${prefix}.gene_panel.tsv.gz" >&2
            exit 1
        fi
    fi

    exit 0
    """
}
