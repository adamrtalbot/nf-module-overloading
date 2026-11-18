nextflow.enable.types = true

include { BwaAligner ; bwaOptions ; BWA_INDEX ; BWA_MEM } from './modules/bwa.nf'

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

process COMPARE_SAMS {
    container 'community.wave.seqera.io/library/samtools:1.22.1--1b65784f3468426d'

    input:
    sams: Bag<Path>
    baseline: String

    output:
    stdout()

    script:
    // each SAM is named <aligner>.sam. Primary records become "read/mate chr pos mapq cigar"
    // rows, which are joined to the baseline's rows on read and mate.
    """
    export LC_ALL=C
    for sam in *.sam; do
        id=\$(basename \$sam .sam)
        samtools view -F 0x900 \$sam \
            | awk '{print \$1 (int(\$2/128)%2 ? "/2" : "/1"), \$3, \$4, \$5, \$6}' \
            | sort -k1,1 > \$id.tsv
    done

    echo "Alignments compared to ${baseline} (primary records matched by read and mate)"
    printf '%-10s %9s %9s %9s %9s %9s %9s\\n' aligner records mapped chr= pos= cigar= mapq=
    for sam in *.sam; do
        id=\$(basename \$sam .sam)
        records=\$(samtools view -c -F 0x900 \$sam)
        mapped=\$(samtools view -c -F 0x904 \$sam)
        join ${baseline}.tsv \$id.tsv \
            | awk -v id=\$id -v r=\$records -v m=\$mapped \
                '{chr+=(\$2==\$6); pos+=(\$3==\$7); mapq+=(\$4==\$8); cigar+=(\$5==\$9)}
                 END {printf "%-10s %9d %9d %9d %9d %9d %9d\\n", id, r, m, chr, pos, cigar, mapq}'
    done
    """
}
