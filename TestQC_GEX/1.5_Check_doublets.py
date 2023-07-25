# python3 1.5_Check_doublets.py --path_to_counts_matrix /home/sevastopol/data/gserranos/H1_kaust/Data/WT_GEX_experiment/outs/filtered_feature_bc_matrix --plots True --sample_name WT &
# python3 1.5_Check_doublets.py --path_to_counts_matrix /home/sevastopol/data/gserranos/H1_kaust/Data/H1X_GEX_experiment/outs/filtered_feature_bc_matrix --plots True --sample_name H1X &

import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse


parser = argparse.ArgumentParser(description="Detects the doublets in the SC data")
parser.add_argument("--path_to_counts_matrix", type=str, help="path to counts matrix")
parser.add_argument("--plots", type=bool, default=False, help="Generate the plots")
parser.add_argument(
    "--sample_name", type=str, default="Sample", help="Name of the sample"
)
args = parser.parse_args()

# Set the graphical parameters
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = "Arial"
plt.rc("font", size=14)
plt.rcParams["pdf.fonttype"] = 42


input_dir = args.path_to_counts_matrix
print(f"Seeking doublets int: {input_dir}")

counts_matrix = scipy.io.mmread(os.path.join(input_dir, "matrix.mtx.gz")).T.tocsc()
genes = np.array(
    scr.load_genes(os.path.join(input_dir, "features.tsv"), delimiter="\t", column=1)
)

print(
    "Counts matrix shape: {} rows, {} columns".format(
        counts_matrix.shape[0], counts_matrix.shape[1]
    )
)

print("Number of genes in gene list: {}".format(len(genes)))
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)

doublet_scores, predicted_doublets = scrub.scrub_doublets(
    min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30
)


if args.plots:
    PLOT_PATH = os.path.join(os.getcwd(), "Plots", "Doublets", args.sample_name)
    if not os.path.exists(PLOT_PATH):
        os.makedirs(PLOT_PATH)
    scrub.plot_histogram()
    plt.title("Histogram scores")
    plt.savefig(os.path.join(PLOT_PATH, "Hist.pdf"), bbox_inches="tight")
    plt.close()

    print("Running UMAP...")
    scrub.set_embedding("UMAP", scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
    scrub.plot_embedding("UMAP", order_points=True)
    plt.title("UMAP")
    plt.savefig(os.path.join(PLOT_PATH, "UMAP.png"), bbox_inches="tight")
    plt.close()

    print("Running tSNE...")
    scrub.set_embedding("tSNE", scr.get_tsne(scrub.manifold_obs_, angle=0.9))
    scrub.plot_embedding("tSNE", order_points=True)
    plt.title("tSNE")
    plt.savefig(os.path.join(PLOT_PATH, "tSNE.png"), bbox_inches="tight")
    plt.close()

    # print("Running ForceAtlas2...")
    # scrub.set_embedding(
    #     "FA", scr.get_force_layout(scrub.manifold_obs_, n_neighbors=5, n_iter=1000)
    # )
    # scrub.plot_embedding("FA", order_points=True)
    # plt.title("ForceAtlas2")
    # plt.savefig(os.path.join(PLOT_PATH, "FA.pdf"), bbox_inches="tight")
    # plt.close()
    print("Done.")
