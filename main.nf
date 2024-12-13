#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.runFile = file(params.runFile) // Input RunFile
params.outputDir = file(params.outputDir) // Output directory

workflow {
    main:
        runScript(params.runFile, params.outputDir)
}

process runScript {
    input:
        path runFile
        path outputDir

    output:
        path "${outputDir}/output_file.html"

    script:
        def runDir = "${task.workDir}/run_${task.hash}"
        """
        mkdir -p ${runDir}
        cd ${runDir}
        cp ${runFile} ./RunFile.R
        Rscript RunFile.R
        mv Reports/*.html ${outputDir}/output_file.html
        """
}

