# nf-module-overloading

One Nextflow module that wraps several BWA flavours (`bwa`, `bwa-mem2`, `bwameth`). An input picks which one runs.

Requires Nextflow 26.10 (record types, `nextflow.enable.types = true`) and Docker.

```bash
nextflow run adamrtalbot/nf-module-overloading
```

## The idea

The three tools have nearly the same CLI: index a fasta, then `mem`-style align two read files. Instead of three modules, the process takes an enum in its input record:

```nextflow
enum BwaAligner { BWA, BWA_MEM2, BWAMETH }

process BWA_MEM {
    container bwaContainer(aligner)

    input:
    record(id: String, aligner: BwaAligner, read1: Path, read2: Path,
           index: Path, fasta_name: String, args: String?)
    ...
}
```

The `aligner` field switches three things:

1. **Container.** `bwaContainer(aligner)` maps the enum to a pinned Seqera Container. The `container` directive can read process inputs.
2. **CLI.** The script block sets a flag per tool and builds the command from it:

   ```nextflow
   def is_bwa      = aligner == BwaAligner.BWA
   def is_bwa_mem2 = aligner == BwaAligner.BWA_MEM2
   def align = is_bwa      ? "bwa mem ..."
             : is_bwa_mem2 ? "bwa-mem2 mem ..."
             :               "bwameth.py --reference ..."
   ```

3. **Index layout.** Each tool writes different index files. `BWA_MEM` checks for one file only the matching tool produces and fails if it is missing, so passing a `bwa` index to `bwa-mem2` is an error.

`args` is optional and passed verbatim to the chosen aligner, so it is tool specific (`-M` for the BWA family, `--do-not-penalize-chimeras` for bwameth).

## Layout

- `modules/bwa.nf`: the enum, `BWA_INDEX` and `BWA_MEM`
- `main.nf`: demo workflow running all three aligners on `data/`
- `data/`: a 20 kb synthetic reference and 200 read pairs
