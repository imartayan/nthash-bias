import argparse
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

parser = argparse.ArgumentParser()
parser.add_argument(
    "-f", "--format", nargs="*", default=["png"], help="output formats (default: png)"
)
parser.add_argument(
    "-l",
    "--label",
    type=str,
    default="plot-bias",
    help="output file prefix; saves label-abs.fmt and label-rel.fmt (default: show)",
)
args = parser.parse_args()

# Each element of `data` is a dict with keys r, k, seed, hist.
# hist[i][j] = # hash pairs where the first had i leading zeros and the second had j.
data = json.loads(input())


def make_fig(entry, relative):
    fig, ax = plt.subplots(1, 1, layout="constrained", figsize=(5.5, 4.5))

    mat = np.array(entry["hist"], dtype=float)
    mat[mat <= 5] = 0

    # Trim all-zero rows and columns to focus on the populated region.
    rows = np.where(mat.any(axis=1))[0]
    cols = np.where(mat.any(axis=0))[0]
    if rows.size == 0:
        plt.close(fig)
        return None
    r0, r1 = int(rows[0]), int(rows[-1])
    c0, c1 = int(cols[0]), int(cols[-1])
    mat = mat[r0 : r1 + 1, c0 : c1 + 1]

    # Row-normalize: P(next #LZ | current #LZ).
    row_sums = mat.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    cond = mat / row_sums

    if relative:
        # Expected marginal under independence: geometric P(lz=j) = 2^-(j+1),
        # with P(lz=32) = 2^-32 for the zero hash.  Renormalize over the
        # displayed column range so both cond and expected sum to 1.
        j_vals = np.arange(c0, c1 + 1)
        expected = np.where(j_vals < 32, 0.5 ** (j_vals + 1), 0.5**32)
        expected /= expected.sum()
        plot_data = cond / expected
        vmax = np.abs(plot_data).max()
        norm = mcolors.TwoSlopeNorm(vmin=0, vcenter=1, vmax=vmax)
        cmap = "bwr"
        cbar_label = "observed probability / expected probability"
    else:
        plot_data = cond
        pos = cond[cond > 0]
        norm = mcolors.LogNorm(vmin=pos.min(), vmax=pos.max()) if pos.size else None
        cmap = None
        cbar_label = "P(next #LZ | current #LZ)"

    im = ax.imshow(
        plot_data,
        origin="upper",
        aspect="auto",
        interpolation="nearest",
        norm=norm,
        cmap=cmap,
    )
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label=cbar_label)

    H, W = cond.shape
    ax.set_yticks(range(H))
    ax.set_yticklabels(range(r0, r1 + 1), fontsize=9)
    ax.set_xticks(range(W))
    ax.set_xticklabels(range(c0, c1 + 1), rotation=90, ha="center", fontsize=9)
    ax.set_xlabel("# leading zeros of next hash")
    ax.set_ylabel("# leading zeros of current hash")

    title = (
        "Bias of leading-zero transition (vs geometric dist.)"
        if relative
        else "Leading-zero transition P(next #LZ | current #LZ)"
    )
    fig.suptitle(f"{title}\n$R={entry["r"]}$, $k={entry["k"]}$")
    return fig


for entry in data:
    for relative, suffix in [(False, "abs"), (True, "rel")]:
        fig = make_fig(entry, relative)
        if fig is None:
            continue
        for fmt in args.format:
            seed = entry["seed"]
            seed_label = f"-s{seed}" if seed else ""
            out = f"{args.label}-R{entry["r"]}-k{entry["k"]}{seed_label}-{suffix}.{fmt}"
            print(f"Saving {out}")
            fig.savefig(out, bbox_inches="tight", dpi=300)
        plt.close(fig)
