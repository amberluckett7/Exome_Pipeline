#!/usr/bin/env nextflow
process FIXMATE {
    cpus 16
    container 'mgibio/samtools:v1.21-noble'
    tag "${id}"

    input:
      tuple val(id), val(sex), val(family), val(trio), val(laneCount), val(famSampleCount), file(bam)

    output:
      tuple val(id), val(sex), val(family), val(trio),
            val(famSampleCount), file("${id}_fx.bam")

    script:
      """
      echo "[${id}] samtools fixmate"
      samtools fixmate -m -@ ${task.cpus} ${bam} - |
      samtools view -@ ${task.cpus} -b -e 'exists([ms])' - |
      samtools sort -@ ${task.cpus} -l 1 -O bam -o ${id}_fx.bam
      """
}