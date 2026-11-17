nextflow.enable.types = true

// Everything that differs between aligners is data in this record. The
// processes read these fields and contain no per-aligner branching.
record BwaOptions {
    container: String
    index_cmd: String
    mem_cmd: String
    rg_flag: String
    ref_flag: String
    index_marker: String
}

enum BwaAligner {
    BWA,
    BWA_MEM2,
    BWA_MEM3,
    BWA_MEM4,
    BWA_MEME,
    BWAMETH,
}

// Pick the options for a known aligner. Callers can also build their own
// BwaOptions record to use a tool that is not listed here.
def bwaOptions(aligner: BwaAligner) -> BwaOptions {
    def all = [(BwaAligner.BWA): record(
        container: 'community.wave.seqera.io/library/bwa:0.7.19--abc3b25d25d1ccb1',
        index_cmd: 'bwa index',
        mem_cmd: 'bwa mem',
        rg_flag: '-R',
        ref_flag: '',
        index_marker: '.bwt',
    ), (BwaAligner.BWA_MEM2): record(
        container: 'community.wave.seqera.io/library/bwa-mem2:2.3--2ed086f836954b3a',
        index_cmd: 'bwa-mem2 index',
        mem_cmd: 'bwa-mem2 mem',
        rg_flag: '-R',
        ref_flag: '',
        index_marker: '.bwt.2bit.64',
    ), (BwaAligner.BWA_MEM3): record(
        container: 'community.wave.seqera.io/library/bwa-mem3:0.17.0--3450d6c4d620ffe3',
        index_cmd: 'bwa-mem3 index',
        mem_cmd: 'bwa-mem3 mem',
        rg_flag: '-R',
        ref_flag: '',
        index_marker: '.bwt.2bit.64',
    ), (BwaAligner.BWA_MEM4): record(
        container: 'community.wave.seqera.io/library/bwa-mem4:4.4.1--945b5c33469f0b7d',
        index_cmd: 'bwa-mem4 index',
        mem_cmd: 'bwa-mem4 mem',
        rg_flag: '-R',
        ref_flag: '',
        index_marker: '.bwt.2bit.64',
    ), (BwaAligner.BWA_MEME): record(
        container: 'community.wave.seqera.io/library/bwa-meme:1.0.6--30966d045013b854',
        index_cmd: 'bwa-meme index',
        mem_cmd: 'bwa-meme mem',
        rg_flag: '-R',
        ref_flag: '',
        index_marker: '.bwt.2bit.64',
    ), (BwaAligner.BWAMETH): record(
        container: 'community.wave.seqera.io/library/bwameth:0.2.10--0e9f5f44844aec75',
        index_cmd: 'bwameth.py index',
        mem_cmd: 'bwameth.py',
        rg_flag: '--read-group',
        ref_flag: '--reference',
        index_marker: '.bwameth.c2t.bwt',
    )]
    def options = all[aligner]
    if (!options) {
        throw new IllegalArgumentException("Unsupported aligner: ${aligner}")
    }
    return options
}

process BWA_INDEX {
    container options.container

    input:
    record(
        id: String,
        fasta: Path,
        options: BwaOptions
    )

    output:
    record(id: id, options: options, index: file('index'))

    script:
    // every tool writes the index next to the fasta, so build inside the output dir
    """
    mkdir index
    cp ${fasta} index/${fasta.name}
    cd index
    ${options.index_cmd} ${fasta.name}
    """
}

process BWA_MEM {
    container options.container

    input:
    record(
        id: String,
        read1: Path,
        read2: Path,
        index: Path,
        fasta_name: String,
        options: BwaOptions,
        args: String?
    )

    output:
    record(id: id, sam: file("${id}.sam"))

    script:
    def rg = "@RG\\tID:${id}\\tSM:${id}"
    def ref = "${index}/${fasta_name}"
    def ref_arg = options.ref_flag ? "${options.ref_flag} ${ref}" : ref
    """
    test -e ${ref}${options.index_marker} || { echo "${index} is not an index for ${options.mem_cmd}" >&2; exit 1; }
    ${options.mem_cmd} -t ${task.cpus} ${options.rg_flag} '${rg}' ${args ?: ''} ${ref_arg} ${read1} ${read2} > ${id}.sam
    """
}
