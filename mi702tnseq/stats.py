"""Statistics for the MI 702 Tn-seq homework."""
import math
import itertools

import numpy as np


def rpkm(counts, gene_length_bp, library_size=None):
    """Reads per kilobase per million. counts may be 1-D or (n_samples, n_genes)."""
    counts = np.asarray(counts, dtype=float)
    length_kb = np.asarray(gene_length_bp, dtype=float) / 1000.0
    if counts.ndim == 1:
        lib = float(counts.sum()) if library_size is None else float(library_size)
        return counts / (lib / 1e6) / length_kb
    lib = counts.sum(axis=1, keepdims=True) if library_size is None \
        else np.asarray(library_size, dtype=float).reshape(-1, 1)
    return counts / (lib / 1e6) / length_kb[None, :]


def signed_fold_change(mean_treated, mean_control):
    """The published convention: a ratio, but values below 1 are written as the
    negative reciprocal, so -4 means fourfold depleted."""
    mt = np.asarray(mean_treated, dtype=float)
    mc = np.asarray(mean_control, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(mc > 0, mt / mc, np.nan)
        return np.where(ratio >= 1, ratio, -1.0 / np.where(ratio > 0, ratio, np.nan))


def benjamini_hochberg(pvalues):
    """BH-corrected q-values, same order as the input."""
    p = np.asarray(pvalues, dtype=float)
    p = np.where(np.isfinite(p), p, 1.0)
    n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(n)
    q[order] = np.clip(ranked, 0, 1)
    return q


def permutation_floor(n_a, n_b, n_tests=1, n_at_floor=1, alpha=0.05):
    """Resolution limits of a permutation test with n_a vs n_b replicates.

    Parameters
    ----------
    n_a, n_b : int
        Replicates in each group.
    n_tests : int
        Number of genes tested, used for the Benjamini-Hochberg correction.
    n_at_floor : int
        How many genes share the smallest p-value. Count this from your own
        results; it changes the achievable q.
    alpha : float
        Significance threshold.

    Returns
    -------
    dict with `distinct_splits`, `min_possible_p`, `genes_tied_at_floor`,
    `q_for_those_genes`, `genes_needed_at_floor` and
    `reachable_with_this_many_ties`.
    """
    total = math.comb(n_a + n_b, n_a)
    distinct = total // 2 if n_a == n_b else total
    min_p = 1.0 / distinct
    q_at_floor = min(1.0, min_p * n_tests / max(1, n_at_floor))
    needed = min_p * n_tests / alpha
    return {
        "distinct_splits": distinct,
        "min_possible_p": min_p,
        "genes_tied_at_floor": n_at_floor,
        "q_for_those_genes": q_at_floor,
        "genes_needed_at_floor": math.ceil(needed),
        "reachable_with_this_many_ties": bool(q_at_floor < alpha),
    }


def resampling_test(control, treated, exact=True, n_perm=10000, seed=0):
    """Permutation test on the difference of means, per gene.

    `control` and `treated` are (n_replicates, n_genes) normalized counts.
    With small equal group sizes `exact=True` enumerates every distinct split,
    which is both faster and exact; Monte-Carlo sampling cannot go below the
    same floor and only adds noise.
    """
    c = np.asarray(control, dtype=float)
    t = np.asarray(treated, dtype=float)
    data = np.vstack([c, t])
    n_c, n_t = c.shape[0], t.shape[0]
    observed = np.abs(t.mean(0) - c.mean(0))

    if exact:
        combos = [s for s in itertools.combinations(range(n_c + n_t), n_c) if 0 in s]
        diffs = np.empty((len(combos), data.shape[1]))
        for i, s in enumerate(combos):
            other = [j for j in range(n_c + n_t) if j not in s]
            diffs[i] = data[other].mean(0) - data[list(s)].mean(0)
        p = (np.abs(diffs) >= observed[None, :]).sum(0) / len(combos)
        n_used = len(combos)
    else:
        rng = np.random.default_rng(seed)
        idx = np.arange(n_c + n_t)
        ge = np.zeros(data.shape[1])
        for _ in range(n_perm):
            pr = rng.permutation(idx)
            d = data[pr[n_c:]].mean(0) - data[pr[:n_c]].mean(0)
            ge += np.abs(d) >= observed
        p = (ge + 1) / (n_perm + 1)
        n_used = n_perm
    return {"p": p, "q": benjamini_hochberg(p), "permutations_used": n_used}


def welch_or_student_test(control, treated, log_transform=True, equal_var=True):
    """Per-gene t-test with BH correction.

    Published Tn-seq screens usually use a regularized test (Cyber-T) that
    borrows variance across genes. A plain t-test on the same table recovers a
    noticeably different gene list at four replicates -- that difference is a
    result, not a mistake.
    """
    import warnings
    from scipy import stats
    c = np.asarray(control, dtype=float)
    t = np.asarray(treated, dtype=float)
    if log_transform:
        c, t = np.log2(c + 1), np.log2(t + 1)
    # Genes with no insertions in either arm have zero variance, so the t-test
    # is undefined and SciPy warns about precision loss. That is expected here,
    # not a problem: those genes get p = 1 just below.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        stat, p = stats.ttest_ind(t, c, axis=0, equal_var=equal_var)
    p = np.where(np.isfinite(p), p, 1.0)
    return {"t": stat, "p": p, "q": benjamini_hochberg(p)}
