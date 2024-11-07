#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MAG_POPGEN } from './workflows/mag_popgen'

workflow {
    MAG_POPGEN()
}
