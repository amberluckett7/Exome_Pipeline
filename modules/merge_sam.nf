process MERGE_BAMS {
    cpus 16
    tag "${meta.id}"
    container 'mgibio/samtools:v1.21-noble'

    input:
      tuple val(meta), file(bams)

    output:
      tuple val(meta.id), val(meta.sex), val(meta.family), val(meta.trio), val(meta.laneCount), val(meta.famSampleCount), file("${meta.id}_merged_ns.bam")

    script:
      """
      echo "[${meta.id}] merging ${bams.size()} lane(s)"
      samtools merge -n -@${task.cpus} -l 1 \
        ${meta.id}_merged_ns.bam ${bams.join(' ')}
      """
}