# nf-module-overloading

One Nextflow module that wraps several BWA flavours. An input record picks which one runs.

Supported: `bwa`, `bwa-mem2`, `bwa-mem3`, `bwa-mem4`, `bwa-meme`, `bwameth`.

Requires Nextflow 26.10 (record types, `nextflow.enable.types = true`) and Docker.

```bash
nextflow run adamrtalbot/nf-module-overloading
```

## The idea

These tools share almost the same CLI: index a fasta, then run `mem` on two read files. The differences are small and are all data: the container, the command names, the read group flag, whether the reference needs a flag, and which file shows the index is the right kind.

So the differences live in a record, and the process reads its fields:

```nextflow
record BwaOptions {
    container: String
    index_cmd: String
    mem_cmd: String
    rg_flag: String
    ref_flag: String
    index_marker: String
}

process BWA_MEM {
    container options.container

    input:
    record(id: String, read1: Path, read2: Path, index: Path,
           fasta_name: String, options: BwaOptions, args: String?)

    script:
    """
    ${options.mem_cmd} -t ${task.cpus} ${options.rg_flag} '${rg}' ${args ?: ''} ${ref_arg} ${read1} ${read2}
    """
}
```

The script has no per-aligner branching. An enum and one lookup function build the record for known tools:

```nextflow
runs = channel.of(
    record(id: 'bwa',     options: bwaOptions(BwaAligner.BWA),     args: '-M'),
    record(id: 'bwameth', options: bwaOptions(BwaAligner.BWAMETH), args: '--do-not-penalize-chimeras'),
)
```

Because the process only sees a `BwaOptions` record, you can add a new BWA-like tool by passing your own record, with no change to the module.

- **Container.** `container options.container` reads a process input.
- **Index check.** `BWA_MEM` checks that the index contains `options.index_marker` and fails if not, so a `bwa` index passed to `bwa-mem2` is an error. `bwa-mem2`, `bwa-mem3`, `bwa-mem4` and `bwa-meme` all write `.bwt.2bit.64`, so the check cannot tell those four apart.
- **Args.** `args` is optional and passed verbatim to the chosen aligner, so it is tool specific.

`bwa-meme` runs with its default index, not the learned index (`-a meme` and `-7`). `bwa-fastalign` is not included: it crashed building an index under amd64 emulation when I tried it.

## Layout

- `modules/bwa.nf`: `BwaOptions`, the `BwaAligner` enum, `bwaOptions()`, `BWA_INDEX` and `BWA_MEM`
- `main.nf`: demo workflow running every aligner on `data/`
- `data/`: a 20 kb synthetic reference and 200 read pairs
