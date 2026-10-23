"""Dataset download and loading for the MI 702 Tn-seq homework.

The data lives as release assets on the course GitHub repository, so a student
session is one download and no authentication.
"""
import os
import zipfile
import urllib.request

#: Set by the course instructor. Release assets are fetched from
#: https://github.com/<REPO>/releases/download/<RELEASE_TAG>/<asset>
REPO = "nkorotkova-uky/mi702-tnseq"
RELEASE_TAG = "v1.0"

ASSETS = {
    "core": "mi702-tnseq-core.zip",     # reference, annotation, tables, wig files
    "reads": "mi702-tnseq-reads.zip",   # multiplexed, demultiplexed and trimmed FASTQ
}

DATA_DIR = "mi702_data"


def asset_url(name):
    return "https://github.com/%s/releases/download/%s/%s" % (
        REPO, RELEASE_TAG, ASSETS[name])


def fetch(name="core", data_dir=DATA_DIR, force=False):
    """Download and unpack one release asset. Returns the data directory.

    'core' (~13 MB) is needed by every part of the homework.
    'reads' (~157 MB) is only needed for Part B.
    """
    os.makedirs(data_dir, exist_ok=True)
    marker = os.path.join(data_dir, ".%s_ready" % name)
    if os.path.exists(marker) and not force:
        print("%s data already present in %s/" % (name, data_dir))
        return data_dir
    url = asset_url(name)
    local = os.path.join(data_dir, ASSETS[name])
    print("downloading %s ..." % ASSETS[name])
    urllib.request.urlretrieve(url, local)
    size_mb = os.path.getsize(local) / 1048576
    with zipfile.ZipFile(local) as z:
        z.extractall(data_dir)
    os.remove(local)
    open(marker, "w").close()
    print("unpacked %.1f MB into %s/" % (size_mb, data_dir))
    return data_dir


def path(*parts, data_dir=DATA_DIR):
    """Build a path inside the data directory, with a clear error if absent."""
    p = os.path.join(data_dir, *parts)
    if not os.path.exists(p):
        raise FileNotFoundError(
            "%s not found. Run mi702tnseq.data.fetch('core') and, for read "
            "files, fetch('reads')." % p)
    return p


def load_published_table(data_dir=DATA_DIR):
    """Published full-depth per-gene RPKM, 1,080 genes x 8 replicates.

    This is the table deposited with the paper. It has NO fold changes or
    p-values -- computing those is Part C.
    """
    import pandas as pd
    return pd.read_csv(path("GAS_hGIIA_gene_RPKM.csv", data_dir=data_dir))


def load_subset_counts(data_dir=DATA_DIR):
    """Per-gene insertion counts computed from the 300k-read-per-sample subset."""
    import pandas as pd
    return pd.read_csv(path("SUBSET_gene_counts.csv", data_dir=data_dir))


def load_ta_sites(data_dir=DATA_DIR):
    """0-based positions of every TA dinucleotide in the genome."""
    out = []
    with open(path("TA_sites_CP008776.txt", data_dir=data_dir)) as fh:
        for line in fh:
            if not line.startswith("#"):
                out.append(int(line) - 1)
    return out


def load_annotation(data_dir=DATA_DIR):
    """All 1,760 annotated CDS from the RefSeq record, as a DataFrame.

    Deliberately larger than the results table -- the difference is Part D.
    """
    import pandas as pd
    rows = []
    with open(path("GAS5448.prot_table", data_dir=data_dir)) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 9:
                rows.append({"product": f[0], "start": int(f[1]), "end": int(f[2]),
                             "strand": f[3], "aa_len": int(f[4]),
                             "gene": f[7], "locus_tag": f[8]})
    return pd.DataFrame(rows)
