nextflow.enable.types = true

include { BwaAligner ; bwaOptions ; BWA_INDEX ; BWA_MEM } from './modules/bwa.nf'
include { COMPARE_SAMS } from './modules/compare_sams.nf'

// Demo: the same two processes run with every supported aligner.
workflow {
    def fasta = file("${projectDir}/data/ref.fa")
    def read1 = file("${projectDir}/data/r1.fq")
    def read2 = file("${projectDir}/data/r2.fq")

    // args go straight to the chosen aligner, so they are aligner specific
    runs = channel.of(
        record(id: 'bwa', options: bwaOptions(BwaAligner.BWA), args: '-M'),
        record(id: 'bwa-mem2', options: bwaOptions(BwaAligner.BWA_MEM2), args: '-M'),
        record(id: 'bwa-mem3', options: bwaOptions(BwaAligner.BWA_MEM3), args: '-M'),
        record(id: 'bwa-mem4', options: bwaOptions(BwaAligner.BWA_MEM4), args: '-M'),
        record(id: 'bwa-meme', options: bwaOptions(BwaAligner.BWA_MEME), args: '-M'),
        record(id: 'bwameth', options: bwaOptions(BwaAligner.BWAMETH), args: '--do-not-penalize-chimeras'),
    )

    indexes = BWA_INDEX(runs.map { r -> record(id: r.id, fasta: fasta, options: r.options) })

    // join the index back to its run by id
    sams = BWA_MEM(
        runs.join(indexes, by: 'id').map { r ->
            record(
                id: r.id,
                read1: read1,
                read2: read2,
                index: r.index,
                fasta_name: fasta.name,
                options: r.options,
                args: r.args,
            )
        }
    )

    // compare every aligner against bwa
    COMPARE_SAMS(sams.map { r -> r.sam }.collect(), 'bwa').view()
}
