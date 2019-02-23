#!/usr/bin/env nextflow

Channel
    .fromPath(params.input)
    .splitCsv(header:true, sep:'\t')
    .map{ row-> tuple(row.sampleId, file(row.read1), file(row.read2)) }
    .into { samples_fastqc_ch; samples_shovill_ch; samples_kraken2_ch }


def summary = [:]
summary['Pipeline Name']  = 'cpo'
summary['Input']        = params.input
summary['Output dir']   = params.outdir
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile


summary.each{ k, v -> println "${k}: ${v}" }

//process printInput {
//    tag $sampleId
//    publishDir "$params.outdir/$sampleId", mode: 'copy'
//    input:
//	set sampleId, file(read1), file(read2) from print_input_ch
//    output:
//	file 'output.txt'
//    script:
//    """
//    echo $sampleId $read1 $read2 > output.txt
//    """
//}

/*
 * FastQC
 */
process fastqc {
    tag "$sampleId"
    cpus 4
    conda 'fastqc=0.11.8'
    container 'quay.io/biocontainers/fastqc:0.11.8--1'
    publishDir "${params.outdir}/${sampleId}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
	set sampleId, file(read1), file(read2) from samples_fastqc_ch

    output:
	file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    fastqc -q $read1 $read2
    """
}


/*
 * shovill
 */
process shovill {
    cpus 8
    memory '12 GB'
    tag "$sampleId"
    conda 'shovill=1.0.1'
    container 'quay.io/biocontainers/shovill:1.0.1--0'
    publishDir "${params.outdir}/${sampleId}/shovill", mode: 'copy',
	saveAs: {filename -> 'contigs.fa'}

    input:
	set val(sampleId), file(read1), file(read2) from samples_shovill_ch

    output:
	set val("$sampleId"), file("output/contigs.fa") into assembly_mlst_ch

    script:
    """
    shovill --cpus 8 --outdir output --R1 $read1 --R2 $read2 
    """
}

process kraken2 {
    cpus 8
    memory '64 GB'
    tag "${sampleId}"
    conda 'kraken2=2.0.7_beta'
    container 'quay.io/biocontainers/kraken2:2.0.7_beta--pl526h2d50403_0'
    publishDir "${params.outdir}/${sampleId}/kraken2", mode: 'copy',
	saveAs: {filename -> 'kraken2_report.txt'}
    input:
	set sampleId, file(read1), file(read2) from samples_kraken2_ch
    output:
	file "report.txt" into kraken2_report_ch
    script:
    """
    kraken2 --db /data/ref_databases/kraken2/minikraken2_v2_8GB --threads 8 --output "-" --report report.txt --paired ${read1} ${read2}
    """
}

process mlst {
    cpus 1
    memory '2 GB'
    tag "$sampleId"
    conda 'mlst=2.15.1'
    container 'quay.io/biocontainers/mlst:2.15.1--0'
    publishDir "${params.outdir}/${sampleId}/mlst", mode: 'copy',
	saveAs: {filename -> 'mlst_report.txt'}
    input:
	set sampleId, file(contig) from assembly_mlst_ch
    output:
	file "report.txt" into mlst_report_ch
    script:
    """
    mlst $contig > report.txt
    """
}
