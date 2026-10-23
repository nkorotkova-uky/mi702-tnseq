"""Tn-seq read-processing pipeline for the MI 702 homework.

Pure Python, no external aligner, so it runs anywhere Colab runs.

Read architecture of ENA run ERR2683524 (recovered empirically; the archive
metadata does not document it):

    [0:6]      6 bp inline Illumina index (sample barcode)
    [6:j-2]    14-15 bp genomic DNA
    [j-2:j]    TA  <- the Himar1/Krmit insertion site
    [j:]       ACAGGTTGGATGATAAGTCCCCGG... (transposon repeat + adapter)

The insertion site is the TA at the 3' END of the genomic tag, not the read
start.
"""
import gzip
import collections

ITR = "ACAGGTTGGATGATAAGTCCCCGG"
BARCODE_LENGTH = 6

#: Barcode -> sample. Recovered by correlating per-barcode gene profiles
#: against the published per-replicate RPKM columns of S2 Table; all eight
#: best matches formed a bijection in TruSeq index order.
BARCODE_TO_SAMPLE = {
    "ATCACG": "control_rep1",   # TruSeq index 1
    "CGATGT": "control_rep2",   # TruSeq index 2
    "TTAGGC": "control_rep3",   # TruSeq index 3
    "GCCAAT": "control_rep4",   # TruSeq index 6
    "CAGATC": "hGIIA_rep1",     # TruSeq index 7
    "ACTTGA": "hGIIA_rep2",     # TruSeq index 8
    "GATCAG": "hGIIA_rep3",     # TruSeq index 9
    "GGCTAC": "hGIIA_rep4",     # TruSeq index 11
}

SAMPLES = ["control_rep1", "control_rep2", "control_rep3", "control_rep4",
           "hGIIA_rep1", "hGIIA_rep2", "hGIIA_rep3", "hGIIA_rep4"]
CONTROLS = SAMPLES[:4]
TREATED = SAMPLES[4:]


def revcomp(seq):
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def read_fastq(path, limit=None):
    """Yield (header, sequence, quality) from a plain or gzipped FASTQ."""
    opener = gzip.open if str(path).endswith(".gz") else open
    n = 0
    with opener(path, "rt") as fh:
        while True:
            head = fh.readline()
            if not head:
                return
            seq = fh.readline().rstrip("\n")
            fh.readline()
            qual = fh.readline().rstrip("\n")
            if not qual:
                return
            yield head.rstrip("\n"), seq, qual
            n += 1
            if limit and n >= limit:
                return


def load_genome(path):
    """Return the genome as one uppercase string."""
    parts = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.startswith(">"):
                parts.append(line.strip())
    return "".join(parts).upper()


def find_ta_sites(genome):
    """0-based positions of the T of every TA dinucleotide."""
    out = []
    i = genome.find("TA")
    while i != -1:
        out.append(i)
        i = genome.find("TA", i + 1)
    return out


def describe_reads(path, n=50000):
    """Summarize read structure: lengths, barcode shares, junction dinucleotides.

    Use this to SEE the architecture rather than taking it on trust.
    """
    reads = [s for _, s, _ in read_fastq(path, limit=n)]
    if not reads:
        raise ValueError("no reads read from %s" % path)
    bc = collections.Counter(s[:BARCODE_LENGTH] for s in reads)
    junction = collections.Counter()
    offsets = collections.Counter()
    for s in reads:
        j = s.find(ITR)
        if j >= 2:
            junction[s[j - 2:j]] += 1
            offsets[j] += 1
    known = {b: c for b, c in bc.items() if b in BARCODE_TO_SAMPLE}
    return {
        "reads_examined": len(reads),
        "read_lengths": collections.Counter(len(s) for s in reads).most_common(3),
        "known_barcode_fraction": round(sum(known.values()) / len(reads), 4),
        "reads_per_barcode": {BARCODE_TO_SAMPLE[b]: c for b, c in sorted(known.items())},
        "transposon_found_fraction": round(sum(offsets.values()) / len(reads), 4),
        "transposon_start_positions": offsets.most_common(4),
        "junction_dinucleotides": junction.most_common(4),
    }


def dissect_read(seq):
    """Split one read into its four parts. Returns None if it has no transposon."""
    j = seq.find(ITR)
    if j < BARCODE_LENGTH + 2:
        return None
    return {
        "barcode": seq[:BARCODE_LENGTH],
        "sample": BARCODE_TO_SAMPLE.get(seq[:BARCODE_LENGTH], "unknown barcode"),
        "genomic_dna": seq[BARCODE_LENGTH:j - 2],
        "insertion_dinucleotide": seq[j - 2:j],
        "transposon_and_adapter": seq[j:],
        "tag_for_mapping": seq[BARCODE_LENGTH:j],
    }


def demultiplex(src, outdir, per_sample_cap=None, progress_every=500000):
    """Split a multiplexed FASTQ by inline barcode into <outdir>/<sample>.fastq.gz.

    Measured at about 80 s for the 2.4 M-read course file; allow up to roughly
    twice that on a slower machine. Progress is printed so you can see it
    working. Returns (tallies, reads per sample).
    """
    import os
    import sys
    os.makedirs(outdir, exist_ok=True)
    handles = {b: gzip.open(os.path.join(outdir, "%s.fastq.gz" % s), "wt")
               for b, s in BARCODE_TO_SAMPLE.items()}
    tally = collections.Counter()
    written = collections.Counter()
    try:
        for head, seq, qual in read_fastq(src):
            tally["reads_read"] += 1
            if progress_every and tally["reads_read"] % progress_every == 0:
                print("  %.1f M reads processed" % (tally["reads_read"] / 1e6))
                sys.stdout.flush()
            b = seq[:BARCODE_LENGTH]
            if b not in BARCODE_TO_SAMPLE:
                tally["unrecognized_barcode"] += 1
                continue
            if per_sample_cap and written[b] >= per_sample_cap:
                tally["over_cap"] += 1
                continue
            handles[b].write("%s\n%s\n+\n%s\n" % (head, seq, qual))
            written[b] += 1
            tally["assigned"] += 1
    finally:
        for h in handles.values():
            h.close()
    return dict(tally), {BARCODE_TO_SAMPLE[b]: n for b, n in sorted(written.items())}


def trim(src, dst, min_tag=15, max_tag=20):
    """Remove the barcode and everything from the transposon onward, keeping TA.

    Returns (reads in, reads kept). A retention near 100% means the transposon
    was NOT removed; around 73-76% is expected here.
    """
    import os
    parent = os.path.dirname(dst)
    if parent:
        os.makedirs(parent, exist_ok=True)
    total = kept = 0
    with gzip.open(dst, "wt") as out:
        for head, seq, qual in read_fastq(src):
            total += 1
            j = seq.find(ITR)
            if j < BARCODE_LENGTH + 2:
                continue
            if seq[j - 2:j] != "TA":
                continue
            tag = seq[BARCODE_LENGTH:j]
            if not (min_tag <= len(tag) <= max_tag) or "N" in tag:
                continue
            out.write("%s\n%s\n+\n%s\n" % (head, tag, qual[BARCODE_LENGTH:j]))
            kept += 1
    return total, kept


def build_ta_index(genome, ta_sites, k=12):
    """Index the last k bases of every tag that could end at a TA site.

    Insertions can only occur at TA sites, so indexing the ~133k TA contexts
    (both strands) instead of all ~1.8M genomic k-mers is about nine times
    smaller and four times faster to build, and it handles tags shorter than
    the k-mer length.
    """
    idx = {}
    n = len(genome)
    for p in ta_sites:
        if p + 2 - k >= 0:
            idx.setdefault(genome[p + 2 - k:p + 2], []).append((p, "+"))
        if p + k <= n:
            idx.setdefault(revcomp(genome[p:p + k]), []).append((p, "-"))
    return idx


def map_insertions(trimmed_fastq, genome, index, k=12):
    """Match trimmed tags to TA sites on both strands.

    Returns (Counter of 0-based TA positions, tallies). A unique rate far below
    ~90% means the barcode boundary is wrong, not that the reference is wrong.
    """
    sites = collections.Counter()
    tally = collections.Counter()
    n = len(genome)
    for _, tag, _ in read_fastq(trimmed_fastq):
        L = len(tag)
        if L < k:
            tally["too_short"] += 1
            continue
        hits = []
        for p, strand in index.get(tag[-k:], ()):
            if strand == "+":
                s = p + 2 - L
                if s >= 0 and genome[s:p + 2] == tag:
                    hits.append(p)
            else:
                if p + L <= n and revcomp(genome[p:p + L]) == tag:
                    hits.append(p)
        if len(hits) == 1:
            sites[hits[0]] += 1
            tally["unique"] += 1
        elif hits:
            tally["ambiguous"] += 1
        else:
            tally["unmapped"] += 1
    return sites, dict(tally)


def read_wig(path):
    """Load a TRANSIT-format .wig into {0-based TA position: count}."""
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith(("#", "variableStep", "track")):
                continue
            parts = line.split()
            if len(parts) >= 2:
                out[int(parts[0]) - 1] = int(parts[1])
    return collections.Counter(out)


def write_wig(site_counts, ta_sites, path, track_name, chrom="CP008776.1"):
    """Write one row per genomic TA site, 1-based position and count."""
    with open(path, "w") as fh:
        fh.write("# %s\n" % track_name)
        fh.write("variableStep chrom=%s\n" % chrom)
        for p in ta_sites:
            fh.write("%d %d\n" % (p + 1, site_counts.get(p, 0)))
    return path


def count_per_gene(site_counts, genes):
    """Sum insertions per gene.

    `genes` is a DataFrame with start_90pct / stop_90pct columns (1-based,
    inclusive) -- the published analysis discarded the final 10% of each gene
    because insertions there often leave a functional protein.
    """
    import numpy as np
    pos = np.fromiter(sorted(site_counts), dtype=np.int64)
    val = np.fromiter((site_counts[p] for p in sorted(site_counts)), dtype=np.int64)
    pos1 = pos + 1
    cum = np.concatenate([[0], np.cumsum(val)])
    lo = np.searchsorted(pos1, genes["start_90pct"].to_numpy(float), side="left")
    hi = np.searchsorted(pos1, genes["stop_90pct"].to_numpy(float), side="right")
    return (cum[hi] - cum[lo]).astype(int)
