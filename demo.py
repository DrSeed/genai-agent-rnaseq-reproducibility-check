# Self-contained demo: simulate counts, run DE, verify, and plot a volcano.
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from analyse_rnaseq import differential_expression, verify


def simulate(n_genes=800, n_per_group=4, n_true=40, seed=0):
    rng = np.random.default_rng(seed)
    base = rng.uniform(20, 400, size=n_genes)
    samples = ["ctrl_%d" % i for i in range(n_per_group)]
    samples += ["treat_%d" % i for i in range(n_per_group)]
    mat = np.zeros((n_genes, len(samples)))
    effect = np.zeros(n_genes)
    true_idx = rng.choice(n_genes, size=n_true, replace=False)
    effect[true_idx] = rng.choice([-1, 1], size=n_true) * rng.uniform(1.5, 3.0, size=n_true)
    for j, s in enumerate(samples):
        is_treat = s.startswith("treat")
        mu = base * (2.0 ** (effect if is_treat else 0.0))
        mat[:, j] = rng.poisson(np.clip(mu, 1, None))
    genes = ["gene_%04d" % i for i in range(n_genes)]
    counts = pd.DataFrame(mat, index=genes, columns=samples)
    groups = pd.DataFrame({
        "sample": samples,
        "group": ["ctrl"] * n_per_group + ["treat"] * n_per_group,
    }).set_index("sample")
    return counts, groups, set(genes[i] for i in true_idx)


def main():
    os.makedirs("figures", exist_ok=True)
    os.makedirs("results", exist_ok=True)

    counts, groups, truth = simulate()
    counts.to_csv("results/demo_counts.csv")
    groups.to_csv("results/demo_groups.csv")

    res, levels, a_cols, b_cols = differential_expression(counts, groups)
    res = verify(counts, groups, res, levels, a_cols, b_cols)

    res["is_true"] = [g in truth for g in res.index]
    res["neglog10p"] = -np.log10(res["pvalue"].clip(lower=1e-300))

    sig = res[res["pvalue"] < 0.05]
    n_true_found = int(sig["is_true"].sum())
    summary = pd.DataFrame({
        "metric": ["n_genes", "n_planted_true", "n_significant", "n_true_recovered", "n_direction_flags"],
        "value": [len(res), len(truth), len(sig), n_true_found, int((~res["verified_direction"]).sum())],
    })
    summary.to_csv("results/summary.csv", index=False)

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(res["logFC"], res["neglog10p"], s=8, c="#bbbbbb", label="other genes")
    tp = res[res["is_true"]]
    ax.scatter(tp["logFC"], tp["neglog10p"], s=18, c="#cc3333", label="planted true DE")
    ax.axhline(-np.log10(0.05), color="#333333", ls="--", lw=1)
    ax.set_xlabel("log2 fold change (treat vs ctrl)")
    ax.set_ylabel("-log10 p-value")
    ax.set_title("Volcano from simulated RNA-seq (verified DE)")
    ax.legend()
    fig.tight_layout()
    fig.savefig("figures/demo.png", dpi=120)
    print("Saved figures/demo.png and results/summary.csv")


if __name__ == "__main__":
    main()
