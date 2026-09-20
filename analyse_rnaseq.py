# Reusable differential-expression scaffold with an independent verification pass.
import argparse
import os
import numpy as np
import pandas as pd


def log_cpm(counts):
    lib = counts.sum(axis=0)
    cpm = counts.divide(lib, axis=1) * 1e6
    return np.log2(cpm + 1.0)


def differential_expression(counts, groups):
    # Agent-drafted path: simple two-group comparison on log-CPM.
    lcpm = log_cpm(counts)
    g = groups["group"]
    levels = sorted(g.unique())
    a_cols = g[g == levels[0]].index
    b_cols = g[g == levels[1]].index
    a = lcpm[a_cols].values
    b = lcpm[b_cols].values
    mean_a = a.mean(axis=1)
    mean_b = b.mean(axis=1)
    lfc = mean_b - mean_a
    va = a.var(axis=1, ddof=1)
    vb = b.var(axis=1, ddof=1)
    na, nb = a.shape[1], b.shape[1]
    se = np.sqrt(va / na + vb / nb) + 1e-8
    t = lfc / se
    from math import erf
    # two-sided p using a normal approximation to keep dependencies light
    p = 2.0 * (1.0 - np.array([0.5 * (1.0 + erf(abs(x) / np.sqrt(2))) for x in t]))
    res = pd.DataFrame({
        "gene": counts.index,
        "logFC": lfc,
        "stat": t,
        "pvalue": p,
    }).set_index("gene")
    return res, levels, a_cols, b_cols


def verify(counts, groups, res, levels, a_cols, b_cols):
    # Independent verification: recompute direction from raw mean counts.
    raw_a = counts[a_cols].mean(axis=1) + 1.0
    raw_b = counts[b_cols].mean(axis=1) + 1.0
    raw_dir = np.sign(np.log2(raw_b / raw_a))
    model_dir = np.sign(res["logFC"].values)
    disagree = (raw_dir.values != model_dir) & (res["pvalue"].values < 0.05)
    res["verified_direction"] = ~disagree
    n_flag = int(disagree.sum())
    print("Verification: %d significant genes with mismatched direction" % n_flag)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True, help="CSV, genes x samples")
    ap.add_argument("--groups", required=True, help="CSV with columns sample,group")
    ap.add_argument("--out", default="results/de_results.csv")
    args = ap.parse_args()

    counts = pd.read_csv(args.counts, index_col=0)
    groups = pd.read_csv(args.groups, index_col=0)
    res, levels, a_cols, b_cols = differential_expression(counts, groups)
    res = verify(counts, groups, res, levels, a_cols, b_cols)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    res.sort_values("pvalue").to_csv(args.out)
    print("Wrote %s" % args.out)


if __name__ == "__main__":
    main()
