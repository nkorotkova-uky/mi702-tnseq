"""Figures for the MI 702 Tn-seq homework."""
import numpy as np
import matplotlib.pyplot as plt

GRAY = "#6b7280"
DOWN = "#1b6ca8"
UP = "#c1440e"


def volcano(log2_fold_change, qvalues, labels=None, alpha=0.05, title=None,
            label_genes=None, ax=None):
    """Volcano plot: effect size against significance.

    Returns (fig, ax). Genes passing `alpha` are colored by direction;
    `label_genes` is a list of gene names to annotate.
    """
    lfc = np.asarray(log2_fold_change, dtype=float)
    q = np.asarray(qvalues, dtype=float)
    y = -np.log10(np.clip(q, 1e-12, None))
    sig = q < alpha
    up, dn = sig & (lfc > 0), sig & (lfc < 0)

    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 5))
    else:
        fig = ax.figure
    ax.scatter(lfc[~sig], y[~sig], s=10, c=GRAY, alpha=0.4, lw=0,
               label="not significant")
    ax.scatter(lfc[dn], y[dn], s=40, c=DOWN, lw=0.6, edgecolor="white", zorder=3,
               label="depleted, n=%d" % int(dn.sum()))
    ax.scatter(lfc[up], y[up], s=40, c=UP, lw=0.6, edgecolor="white", zorder=3,
               label="enriched, n=%d" % int(up.sum()))
    ax.axhline(-np.log10(alpha), color=GRAY, lw=0.8, ls="--", zorder=1)

    if label_genes and labels is not None:
        labels = np.asarray(labels, dtype=object)
        for name in label_genes:
            hit = np.where(labels == name)[0]
            if hit.size:
                i = int(hit[0])
                ax.annotate(name, (lfc[i], y[i]), textcoords="offset points",
                            xytext=(8, 6), fontsize=8, style="italic",
                            arrowprops=dict(arrowstyle="-", lw=0.5, color=GRAY,
                                            shrinkA=0, shrinkB=2))
    ax.set_xlabel("log$_2$ fold change (hGIIA vs control)")
    ax.set_ylabel("$-$log$_{10}$ corrected p")
    if title:
        ax.set_title(title, loc="left")
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    ax.margins(0.1)
    fig.tight_layout()
    return fig, ax


def insertion_track(site_counts_by_sample, start, end, genes=None, ax=None,
                    gene_label_col="gene_name"):
    """Insertion counts across a genomic window, one lane per sample.

    This is the notebook equivalent of looking at the data in a genome browser.
    `site_counts_by_sample` maps sample name -> {0-based TA position: count}.
    """
    names = list(site_counts_by_sample)
    n = len(names)
    if ax is None:
        fig, axes = plt.subplots(n + (1 if genes is not None else 0), 1,
                                 figsize=(10, 1.0 * n + 1.4), sharex=True,
                                 gridspec_kw={"hspace": 0.15})
    else:
        fig, axes = ax.figure, np.atleast_1d(ax)
    axes = np.atleast_1d(axes)

    top = 0
    for name in names:
        c = site_counts_by_sample[name]
        vals = [v for p, v in c.items() if start <= p + 1 <= end]
        if vals:
            top = max(top, max(vals))
    top = max(top, 1)

    for i, name in enumerate(names):
        a = axes[i]
        c = site_counts_by_sample[name]
        pts = [(p + 1, v) for p, v in c.items() if start <= p + 1 <= end and v > 0]
        if pts:
            xs, ys = zip(*pts)
            a.vlines(xs, 0, ys, lw=1.0,
                     color=DOWN if name.startswith("control") else UP)
        a.set_ylim(0, top * 1.1)
        a.set_ylabel(name.replace("_rep", " "), rotation=0, ha="right",
                     va="center", fontsize=7)
        a.spines[["top", "right"]].set_visible(False)
        a.tick_params(labelsize=7)

    if genes is not None:
        a = axes[-1]
        sub = genes[(genes["stop_90pct"] >= start) & (genes["start_90pct"] <= end)]
        for r in sub.itertuples(index=False):
            a.add_patch(plt.Rectangle((r.start_90pct, 0.25),
                                      r.stop_90pct - r.start_90pct, 0.5,
                                      facecolor="#94a3b8", edgecolor="white"))
            nm = getattr(r, gene_label_col, None)
            if isinstance(nm, str) and nm.strip():
                a.text((r.start_90pct + r.stop_90pct) / 2, 0.5, nm, ha="center",
                       va="center", fontsize=6.5, style="italic")
        a.set_ylim(0, 1)
        a.set_yticks([])
        a.set_ylabel("genes", rotation=0, ha="right", va="center", fontsize=7)
        a.spines[["top", "right", "left"]].set_visible(False)

    axes[-1].set_xlim(start, end)
    axes[-1].set_xlabel("position on CP008776.1 (bp)")
    # Not tight_layout(): the gene panel uses raw Rectangle patches, which it
    # cannot measure, and it warns instead of laying out correctly.
    fig.subplots_adjust(left=0.13, right=0.98, top=0.97, bottom=0.10, hspace=0.15)
    return fig, axes
