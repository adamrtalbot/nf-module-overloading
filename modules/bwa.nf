nextflow.enable.types = true

process BWA_INDEX {
    container bwaContainer(aligner)

    input:
    record(
        id: String,
        aligner: BwaAligner,
        fasta: Path
    )

    output:
    record(id: id, aligner: aligner, index: file('index'))

    script:
    def is_bwa = aligner == BwaAligner.BWA
    def is_bwa_mem2 = aligner == BwaAligner.BWA_MEM2
    def cmd = is_bwa ? 'bwa index' : is_bwa_mem2 ? 'bwa-mem2 index' : 'bwameth.py index'
    // all three write the index next to the fasta, so build inside the output dir
    """
    mkdir index
    cp ${fasta} index/${fasta.name}
    cd index
    ${cmd} ${fasta.name}
    """
}

process BWA_MEM {
    container bwaContainer(aligner)

    input:
    record(
        id: String,
        aligner: BwaAligner,
        read1: Path,
        read2: Path,
        index: Path,
        fasta_name: String,
        args: String?
    )

    output:
    record(id: id, aligner: aligner, sam: file("${id}.sam"))

    script:
    def is_bwa = aligner == BwaAligner.BWA
    def is_bwa_mem2 = aligner == BwaAligner.BWA_MEM2
    def extra = args ?: ''
    def rg = "@RG\\tID:${id}\\tSM:${id}"
    def ref = "${index}/${fasta_name}"
    // one index file that only this aligner produces
    def marker = is_bwa ? '.bwt' : is_bwa_mem2 ? '.bwt.2bit.64' : '.bwameth.c2t.bwt'
    def align = is_bwa
        ? "bwa mem -t ${task.cpus} -R '${rg}' ${extra} ${ref} ${read1} ${read2}"
        : is_bwa_mem2
            ? "bwa-mem2 mem -t ${task.cpus} -R '${rg}' ${extra} ${ref} ${read1} ${read2}"
            : "bwameth.py -t ${task.cpus} --read-group '${rg}' ${extra} --reference ${ref} ${read1} ${read2}"
    """
    test -e ${ref}${marker} || { echo "${index} is not a ${aligner} index" >&2; exit 1; }
    ${align} > ${id}.sam
    """
}

def bwaContainer(aligner: BwaAligner) -> String {
    if (aligner == BwaAligner.BWA) {
        return 'community.wave.seqera.io/library/bwa:0.7.19--abc3b25d25d1ccb1'
    }
    if (aligner == BwaAligner.BWA_MEM2) {
        return 'community.wave.seqera.io/library/bwa-mem2:2.3--2ed086f836954b3a'
    }
    if (aligner == BwaAligner.BWAMETH) {
        return 'community.wave.seqera.io/library/bwameth:0.2.10--0e9f5f44844aec75'
    }
    throw new IllegalArgumentException("Unsupported aligner: ${aligner}")
}

// One module covering several BWA flavours. The `aligner` field of the input
// record selects the container, CLI and expected index layout.
enum BwaAligner {
    BWA,
    BWA_MEM2,
    BWAMETH,
}
