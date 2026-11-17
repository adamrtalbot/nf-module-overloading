nextflow.enable.types = true

include { BwaAligner ; BWA_INDEX ; BWA_MEM } from './modules/bwa.nf'

// Demo: the same two processes run with three different aligners.
workflow {
    def fasta = file("${projectDir}/data/ref.fa")
    def read1 = file("${projectDir}/data/r1.fq")
    def read2 = file("${projectDir}/data/r2.fq")

    // args are passed straight to the chosen aligner, so they are aligner specific
    runs = channel.of(
        record(id: 'bwa', aligner: BwaAligner.BWA, args: '-M'),
        record(id: 'bwa-mem2', aligner: BwaAligner.BWA_MEM2, args: '-M'),
        record(id: 'bwameth', aligner: BwaAligner.BWAMETH, args: '--do-not-penalize-chimeras'),
    )

    indexes = BWA_INDEX(runs.map { r -> record(id: r.id, aligner: r.aligner, fasta: fasta) })

    // join the index back to its run by id
    sams = BWA_MEM(
        runs.join(indexes, by: 'id').map { r ->
            record(
                id: r.id,
                aligner: r.aligner,
                read1: read1,
                read2: read2,
                index: r.index,
                fasta_name: fasta.name,
                args: r.args,
            )
        }
    )

    sams.view { r -> "${r.aligner}: ${r.sam.name} ${r.sam.countLines()} lines" }
}
