#!/usr/bin/env nextflow
process KEVLAR_COUNT {
    tag "${id}"
    container 'jamesrusssilsby/exekevlar:latest'
    containerOptions "-B ${params.base}/resources/kevlar"
    input:
      tuple val(id), val(platform), val(sex), val(family), val(trio), val(flowcell), val(laneCount), val(famSampleCount), val(fastqFiles)
    output:
      tuple val(id), val(sex), val(family), val(trio), val(laneCount), val(famSampleCount),file("${id}_${flowcell}.ct")
    script:
      """
      echo "Running kevlar count for sample ${id} on flowcell ${flowcell}"
      kevlar count --memory 70000M --mask ${params.kevlarMask} ${fastqFiles.R1} ${fastqFiles.R2} > ${id}_${flowcell}.ct
      """
}