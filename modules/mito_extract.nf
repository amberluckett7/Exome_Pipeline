process MITO_EXTRACT {
    tag "${id}"
    cpus 4
    container 'docker://mgibio/samtools:v1.21-noble'
    containerOptions "-B ${params.resourcesDir}"

    input:
      tuple val(id), val(sex), val(family), val(trio), val(famSampleCount), file(bam), file(bai), val(mito_contig)

    output:
      tuple val(id), val(sex), val(family), val(trio), val(famSampleCount), file("${id}.mito.bam"), file("${id}.mito.bam.bai"), val(mito_contig)

    script:
      """
      samtools view -@ ${task.cpus} -b ${bam} ${mito_contig} -o ${id}.mito.bam
      samtools index -@ ${task.cpus} ${id}.mito.bam
      """
}
