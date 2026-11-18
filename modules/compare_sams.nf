nextflow.enable.types = true

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
