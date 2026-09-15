process DETECT_MITO_CONTIG {
    tag "${id}"
    cpus 1
    memory '1 GB'
    container 'docker://mgibio/samtools:v1.21-noble'

    input:
      tuple val(id), val(sex), val(family), val(trio), val(famSampleCount), file(bam), file(bai)

    output:
      tuple val(id), env(MITO_CONTIG), emit: contig

    script:
      """
      MITO_CONTIG=\$(samtools view -H ${bam} \\
          | awk -F'\\t' '/^@SQ/ {for(i=1;i<=NF;i++) if(\$i ~ /^SN:/) print substr(\$i,4)}' \\
          | grep -iE '^(chrM|chrMT|MT|M)\$' | head -n1)

      if [[ -z "\$MITO_CONTIG" ]]; then
          echo "ERROR: no mitochondrial contig found in header of ${bam}" >&2
          exit 1
      fi

      echo "Detected mito contig for ${id}: \$MITO_CONTIG" >&2
      """
}
