#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script: plot_tree_identity_heatmap.py

Description:
    Generate an integrated visualization combining:
        1. phylogenetic tree
        2. pairwise sequence identity heatmap
        3. cluster annotation bar

Input:
    - Newick tree file
    - Pairwise identity matrix in Excel format
    - Directory containing cluster_*.txt files

Output:
    - tree_heatmap_cluster.png
    - tree_heatmap_cluster.pdf
    - tree_only.png / .pdf
    - heatmap_only.png / .pdf
    - cluster_only.png / .pdf
    - gene_id_order.txt
    - cluster_assignment.csv
    - reordered_identity_matrix.csv

Usage:
    python plot_tree_identity_heatmap.py \
        --tree tree.nwk \
        --identity identity.xlsx \
        --cluster_dir clusters_90 \
        --output_dir output_tree_heatmap
"""

import argparse
import os
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Rectangle
from Bio import Phylo


# ---------------------------------------------------------------------
# Cluster colors
# Cluster1-Cluster10 use the same colors as the final figure.
# Additional colors are retained as fallback colors if more clusters
# are present in future analyses.
# ---------------------------------------------------------------------

CLUSTER_COLORS = [
    "#EFD3AC",  # Cluster1
    "#C9B4C7",  # Cluster2
    "#EED1CC",  # Cluster3
    "#D4E3DD",  # Cluster4
    "#D6D9B9",  # Cluster5
    "#A1B0AD",  # Cluster6
    "#F4EEAC",  # Cluster7
    "#C9DCC4",  # Cluster8
    "#E69191",  # Cluster9
    "#A09952",  # Cluster10
    "#A8D5BA",
    "#6EC4B6",
    "#B5D3E7",
    "#89BFD7",
    "#6E9AB6",
    "#97A7D5",
    "#CAB7D5",
    "#E7C1B5",
    "#D5B8A8",
    "#B69A85",
]


# ---------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Plot phylogenetic tree, pairwise sequence identity heatmap, "
            "and cluster annotation."
        )
    )

    parser.add_argument(
        "--tree",
        required=True,
        help="Input phylogenetic tree file in Newick format."
    )

    parser.add_argument(
        "--identity",
        required=True,
        help="Pairwise sequence identity matrix in Excel format."
    )

    parser.add_argument(
        "--cluster_dir",
        required=True,
        help="Directory containing cluster files."
    )

    parser.add_argument(
        "--output_dir",
        required=True,
        help="Output directory."
    )

    parser.add_argument(
        "--tree_scale",
        type=float,
        default=2.0,
        help="Scale factor for tree branch lengths in the plot. Default: 2.0"
    )

    parser.add_argument(
        "--vmin",
        type=float,
        default=40,
        help="Minimum value for identity heatmap color scale. Default: 40"
    )

    parser.add_argument(
        "--vmax",
        type=float,
        default=100,
        help="Maximum value for identity heatmap color scale. Default: 100"
    )

    return parser.parse_args()


# ---------------------------------------------------------------------
# Cluster file handling
# ---------------------------------------------------------------------

def read_cluster_files(cluster_dir):
    """
    Read cluster files and return a gene-to-cluster dictionary.

    Supported file names include:
        cluster1.txt
        cluster_1.txt
        Cluster1.txt
        Cluster_1.txt

    Cluster names are normalized as:
        Cluster1
        Cluster2
        ...
        Cluster10
    """

    cluster_dir = Path(cluster_dir)

    if not cluster_dir.exists():
        raise FileNotFoundError(
            f"Cluster directory not found: {cluster_dir}"
        )

    clusters = {}
    cluster_files = []

    # Search explicitly for Cluster1-Cluster10 so that numerical order
    # is preserved and Cluster10 is not placed between Cluster1/Cluster2.
    for cluster_num in range(1, 11):

        possible_names = [
            f"cluster{cluster_num}.txt",
            f"cluster_{cluster_num}.txt",
            f"Cluster{cluster_num}.txt",
            f"Cluster_{cluster_num}.txt",
        ]

        for filename in possible_names:

            cluster_file = cluster_dir / filename

            if cluster_file.exists():
                cluster_files.append(
                    (cluster_num, cluster_file)
                )
                break

    if not cluster_files:
        raise FileNotFoundError(
            f"No cluster files were found in: {cluster_dir}"
        )

    for cluster_num, cluster_file in cluster_files:

        with open(
            cluster_file,
            "r",
            encoding="utf-8"
        ) as f:

            genes = [
                line.strip()
                for line in f
                if line.strip()
            ]

        cluster_name = f"Cluster{cluster_num}"

        for gene in genes:
            clusters[gene] = cluster_name

        print(
            f"Loaded {cluster_name}: "
            f"{len(genes)} genes"
        )

    return clusters


# ---------------------------------------------------------------------
# Tree handling
# ---------------------------------------------------------------------

def get_tree_leaf_order(tree_file):
    """
    Read a Newick tree and return:
        1. leaf order
        2. Bio.Phylo tree object
    """

    with open(
        tree_file,
        "r",
        encoding="utf-8"
    ) as f:

        tree_content = f.read()

    tree = Phylo.read(
        StringIO(tree_content),
        "newick"
    )

    leaves = [
        leaf.name
        for leaf in tree.get_terminals()
    ]

    return leaves, tree


def scale_tree_branches(tree, scale_factor):
    """
    Scale tree branch lengths for visualization only.
    """

    for clade in tree.find_clades():

        if clade.branch_length:
            clade.branch_length *= scale_factor

    return tree


# ---------------------------------------------------------------------
# Matrix reordering
# ---------------------------------------------------------------------

def reorder_matrix_by_tree(
    identity_df,
    tree_leaves,
    clusters
):
    """
    Reorder identity matrix according to phylogenetic tree leaf order.
    """

    available_leaves = [
        leaf
        for leaf in tree_leaves
        if (
            leaf in identity_df.index
            and leaf in identity_df.columns
        )
    ]

    if not available_leaves:
        raise ValueError(
            "No tree leaf IDs were found in the identity matrix. "
            "Please check whether sequence IDs match."
        )

    reordered_df = identity_df.loc[
        available_leaves,
        available_leaves
    ]

    cluster_info = [
        clusters.get(
            leaf,
            "No_Cluster"
        )
        for leaf in available_leaves
    ]

    print(
        f"Tree leaves: {len(tree_leaves)}"
    )

    print(
        "Matched leaves in identity matrix: "
        f"{len(available_leaves)}"
    )

    print(
        "Genes with cluster assignment: "
        f"{sum(c != 'No_Cluster' for c in cluster_info)}"
    )

    return (
        reordered_df,
        cluster_info,
        available_leaves
    )


# ---------------------------------------------------------------------
# Tree plotting
# ---------------------------------------------------------------------

def draw_tree_on_axis(
    tree,
    ax,
    leaves_order
):
    """
    Draw phylogenetic tree on a matplotlib axis.
    """

    Phylo.draw(
        tree,
        axes=ax,
        do_show=False,
        label_func=lambda x: "",
        show_confidence=False,
        branch_labels=lambda x: "",
    )

    ax.set_xlim(
        0,
        tree.total_branch_length() * 1.1
    )

    ax.set_ylim(
        -0.5,
        len(leaves_order) - 0.5
    )

    ax.invert_yaxis()

    ax.set_xticks([])
    ax.set_yticks([])

    for spine in ax.spines.values():
        spine.set_visible(False)

    for line in ax.lines:
        line.set_linewidth(0.5)


# ---------------------------------------------------------------------
# Cluster color handling
# ---------------------------------------------------------------------

def cluster_sort_key(cluster_name):
    """
    Sort cluster names numerically.

    Example:
        Cluster1
        Cluster2
        ...
        Cluster10

    rather than:
        Cluster1
        Cluster10
        Cluster2
    """

    digits = "".join(
        filter(
            str.isdigit,
            cluster_name
        )
    )

    if digits:
        return int(digits)

    return float("inf")


def build_cluster_color_map(cluster_info):
    """
    Assign fixed colors to clusters in numerical cluster order.
    """

    unique_clusters = sorted(
        (
            cluster
            for cluster in set(cluster_info)
            if cluster != "No_Cluster"
        ),
        key=cluster_sort_key,
    )

    cluster_color_map = {
        cluster: CLUSTER_COLORS[
            i % len(CLUSTER_COLORS)
        ]
        for i, cluster
        in enumerate(unique_clusters)
    }

    return cluster_color_map


# ---------------------------------------------------------------------
# Save gene order
# ---------------------------------------------------------------------

def save_gene_order(
    leaves_order,
    output_dir
):
    """
    Save the gene order used in the final figure.
    """

    output_file = (
        Path(output_dir)
        / "gene_id_order.txt"
    )

    with open(
        output_file,
        "w",
        encoding="utf-8"
    ) as f:

        for gene_id in leaves_order:
            f.write(
                f"{gene_id}\n"
            )

    print(
        "Gene order saved to: "
        f"{output_file}"
    )


# ---------------------------------------------------------------------
# Heatmap grid helper
# ---------------------------------------------------------------------

def add_heatmap_grid(
    ax,
    n_genes
):
    """
    Add thin white grid lines between heatmap cells.
    """

    ax.set_xticks(
        np.arange(
            -0.5,
            n_genes,
            1
        ),
        minor=True
    )

    ax.set_yticks(
        np.arange(
            -0.5,
            n_genes,
            1
        ),
        minor=True
    )

    ax.grid(
        which="minor",
        color="white",
        linestyle="-",
        linewidth=0.1
    )

    ax.tick_params(
        which="minor",
        size=0
    )


# ---------------------------------------------------------------------
# Save separate panels
# ---------------------------------------------------------------------

def save_separate_plots(
    tree,
    identity_matrix,
    cluster_info,
    leaves_order,
    cluster_color_map,
    output_dir,
    vmin,
    vmax,
):
    """
    Save tree, heatmap, and cluster annotation separately.
    """

    output_dir = Path(output_dir)

    n_genes = len(
        leaves_order
    )

    # -------------------------------------------------------------
    # Heatmap color scale
    # -------------------------------------------------------------

    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        "identity_cmap",
        [
            "#F2F7FB",
            "#9CB6DD",
            "#4A5989",
        ],
        N=256,
    )

    # -------------------------------------------------------------
    # Tree only
    # -------------------------------------------------------------

    fig_tree = plt.figure(
        figsize=(
            8,
            max(
                4,
                n_genes * 0.05
            )
        )
    )

    ax_tree = fig_tree.add_subplot(
        111
    )

    draw_tree_on_axis(
        tree,
        ax_tree,
        leaves_order
    )

    ax_tree.set_title(
        "Phylogenetic Tree",
        fontsize=14,
        pad=20
    )

    fig_tree.savefig(
        output_dir / "tree_only.png",
        dpi=300,
        bbox_inches="tight",
        facecolor="white"
    )

    fig_tree.savefig(
        output_dir / "tree_only.pdf",
        bbox_inches="tight",
        facecolor="white"
    )

    plt.close(
        fig_tree
    )

    # -------------------------------------------------------------
    # Heatmap only
    # -------------------------------------------------------------

    heatmap_size = min(
        30,
        max(
            6,
            n_genes * 0.03
        )
    )

    fig_heatmap = plt.figure(
        figsize=(
            heatmap_size,
            heatmap_size
        )
    )

    ax_heatmap = fig_heatmap.add_subplot(
        111
    )

    heatmap_data = (
        identity_matrix
        .values
        .astype(float)
    )

    im = ax_heatmap.imshow(
        heatmap_data,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        aspect="auto",
        interpolation="nearest",
    )

    ax_heatmap.set_xticks([])
    ax_heatmap.set_yticks([])

    add_heatmap_grid(
        ax_heatmap,
        n_genes
    )

    ax_heatmap.set_title(
        "Sequence Identity Matrix",
        fontsize=14,
        pad=20
    )

    cbar = plt.colorbar(
        im,
        ax=ax_heatmap,
        fraction=0.046,
        pad=0.04
    )

    cbar.set_label(
        "Sequence Identity (%)",
        rotation=270,
        labelpad=20
    )

    fig_heatmap.savefig(
        output_dir / "heatmap_only.png",
        dpi=300,
        bbox_inches="tight",
        facecolor="white"
    )

    fig_heatmap.savefig(
        output_dir / "heatmap_only.pdf",
        bbox_inches="tight",
        facecolor="white"
    )

    plt.close(
        fig_heatmap
    )

    # -------------------------------------------------------------
    # Cluster annotation only
    # -------------------------------------------------------------

    fig_cluster = plt.figure(
        figsize=(
            2,
            max(
                4,
                n_genes * 0.05
            )
        )
    )

    ax_cluster = fig_cluster.add_subplot(
        111
    )

    for i, cluster in enumerate(
        cluster_info
    ):

        if cluster == "No_Cluster":
            facecolor = "none"
        else:
            facecolor = cluster_color_map[
                cluster
            ]

        rect = Rectangle(
            (0, i),
            1,
            1,
            facecolor=facecolor,
            edgecolor="white",
            linewidth=0.3,
        )

        ax_cluster.add_patch(
            rect
        )

    ax_cluster.set_xlim(
        0,
        1
    )

    ax_cluster.set_ylim(
        0,
        len(cluster_info)
    )

    ax_cluster.set_xticks([])
    ax_cluster.set_yticks([])

    ax_cluster.invert_yaxis()

    ax_cluster.set_title(
        "Cluster Assignment",
        fontsize=14,
        pad=20
    )

    fig_cluster.savefig(
        output_dir / "cluster_only.png",
        dpi=300,
        bbox_inches="tight",
        facecolor="white"
    )

    fig_cluster.savefig(
        output_dir / "cluster_only.pdf",
        bbox_inches="tight",
        facecolor="white"
    )

    plt.close(
        fig_cluster
    )


# ---------------------------------------------------------------------
# Combined plot
# ---------------------------------------------------------------------

def create_combined_plot(
    tree,
    identity_matrix,
    cluster_info,
    leaves_order,
    output_dir,
    vmin,
    vmax,
):
    """
    Create integrated:

        phylogenetic tree
        +
        sequence identity heatmap
        +
        cluster annotation

    figure.
    """

    output_dir = Path(
        output_dir
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    n_genes = len(
        leaves_order
    )

    # -------------------------------------------------------------
    # Figure dimensions
    # -------------------------------------------------------------

    gene_size = 0.1

    heatmap_size = (
        n_genes
        * gene_size
    )

    fig_width = (
        heatmap_size
        * 1.5
    )

    fig_height = (
        heatmap_size
    )

    max_size = 50

    if (
        fig_width > max_size
        or fig_height > max_size
    ):

        scale_factor = (
            max_size
            / max(
                fig_width,
                fig_height
            )
        )

        fig_width *= (
            scale_factor
        )

        fig_height *= (
            scale_factor
        )

    fig = plt.figure(
        figsize=(
            fig_width,
            fig_height
        )
    )

    # -------------------------------------------------------------
    # Layout
    # -------------------------------------------------------------

    gs = fig.add_gridspec(
        1,
        3,
        width_ratios=[
            0.8,
            2,
            0.1
        ],
        wspace=0.02
    )

    # -------------------------------------------------------------
    # Tree
    # -------------------------------------------------------------

    ax_tree = fig.add_subplot(
        gs[0]
    )

    ax_tree.set_title(
        "Phylogenetic Tree",
        fontsize=12,
        pad=10
    )

    draw_tree_on_axis(
        tree,
        ax_tree,
        leaves_order
    )

    # -------------------------------------------------------------
    # Heatmap
    # -------------------------------------------------------------

    ax_heatmap = fig.add_subplot(
        gs[1]
    )

    ax_heatmap.set_title(
        "Sequence Identity Matrix",
        fontsize=12,
        pad=10
    )

    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        "identity_cmap",
        [
            "#F2F7FB",
            "#9CB6DD",
            "#4A5989",
        ],
        N=256,
    )

    heatmap_data = (
        identity_matrix
        .values
        .astype(float)
    )

    im = ax_heatmap.imshow(
        heatmap_data,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        aspect="auto",
        interpolation="nearest",
    )

    ax_heatmap.set_xticks([])
    ax_heatmap.set_yticks([])

    add_heatmap_grid(
        ax_heatmap,
        n_genes
    )

    # -------------------------------------------------------------
    # Cluster annotation
    # -------------------------------------------------------------

    ax_cluster = fig.add_subplot(
        gs[2]
    )

    ax_cluster.set_title(
        "Cluster",
        fontsize=10,
        pad=10
    )

    cluster_color_map = (
        build_cluster_color_map(
            cluster_info
        )
    )

    for i, cluster in enumerate(
        cluster_info
    ):

        if cluster == "No_Cluster":
            facecolor = "none"
        else:
            facecolor = cluster_color_map[
                cluster
            ]

        rect = Rectangle(
            (0, i),
            1,
            1,
            facecolor=facecolor,
            edgecolor="white",
            linewidth=0.1,
        )

        ax_cluster.add_patch(
            rect
        )

    ax_cluster.set_xlim(
        0,
        1
    )

    ax_cluster.set_ylim(
        0,
        len(cluster_info)
    )

    ax_cluster.set_xticks([])
    ax_cluster.set_yticks([])

    ax_cluster.invert_yaxis()

    # -------------------------------------------------------------
    # Heatmap colorbar
    # -------------------------------------------------------------

    cbar_ax = fig.add_axes(
        [
            0.94,
            0.15,
            0.02,
            0.7
        ]
    )

    cbar = plt.colorbar(
        im,
        cax=cbar_ax
    )

    cbar.set_label(
        "Sequence Identity (%)",
        rotation=270,
        labelpad=15,
        fontsize=10
    )

    # -------------------------------------------------------------
    # Cluster legend
    # -------------------------------------------------------------

    if cluster_color_map:

        sorted_clusters = sorted(
            cluster_color_map.keys(),
            key=cluster_sort_key
        )

        legend_elements = [
            Rectangle(
                (0, 0),
                1,
                1,
                facecolor=cluster_color_map[
                    cluster
                ],
                edgecolor="black",
                label=cluster,
            )
            for cluster
            in sorted_clusters
        ]

        n_clusters = len(
            legend_elements
        )

        n_cols = min(
            5,
            n_clusters
        )

        n_rows = (
            n_clusters
            + n_cols
            - 1
        ) // n_cols

        legend_height = (
            0.05
            * n_rows
        )

        bottom_margin = (
            0.05
            + legend_height
        )

        legend_ax = fig.add_axes(
            [
                0.25,
                0.01,
                0.5,
                legend_height
            ]
        )

        legend_ax.axis(
            "off"
        )

        legend = legend_ax.legend(
            handles=legend_elements,
            loc="center",
            ncol=n_cols,
            fontsize=8,
            frameon=True,
            handlelength=1.5,
            handleheight=1.5,
        )

        legend.set_title(
            "Cluster Legend",
            prop={
                "size": 10
            }
        )

    else:

        bottom_margin = 0.05

    # -------------------------------------------------------------
    # Final layout
    # -------------------------------------------------------------

    plt.subplots_adjust(
        left=0.05,
        right=0.93,
        top=0.95,
        bottom=bottom_margin,
        wspace=0.02,
    )

    # -------------------------------------------------------------
    # Save combined figure
    # -------------------------------------------------------------

    fig.savefig(
        output_dir
        / "tree_heatmap_cluster.png",
        dpi=300,
        bbox_inches="tight",
        facecolor="white"
    )

    fig.savefig(
        output_dir
        / "tree_heatmap_cluster.pdf",
        bbox_inches="tight",
        facecolor="white"
    )

    plt.close(
        fig
    )

    # -------------------------------------------------------------
    # Save individual panels
    # -------------------------------------------------------------

    save_separate_plots(
        tree=tree,
        identity_matrix=identity_matrix,
        cluster_info=cluster_info,
        leaves_order=leaves_order,
        cluster_color_map=cluster_color_map,
        output_dir=output_dir,
        vmin=vmin,
        vmax=vmax,
    )

    return cluster_color_map


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():

    args = parse_args()

    tree_file = Path(
        args.tree
    )

    identity_file = Path(
        args.identity
    )

    cluster_dir = Path(
        args.cluster_dir
    )

    output_dir = Path(
        args.output_dir
    )

    # -------------------------------------------------------------
    # Check files
    # -------------------------------------------------------------

    if not tree_file.exists():

        raise FileNotFoundError(
            f"Tree file not found: "
            f"{tree_file}"
        )

    if not identity_file.exists():

        raise FileNotFoundError(
            f"Identity matrix not found: "
            f"{identity_file}"
        )

    if not cluster_dir.exists():

        raise FileNotFoundError(
            f"Cluster directory not found: "
            f"{cluster_dir}"
        )

    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    # -------------------------------------------------------------
    # Load cluster information
    # -------------------------------------------------------------

    print(
        "Loading cluster information..."
    )

    clusters = read_cluster_files(
        cluster_dir
    )

    # -------------------------------------------------------------
    # Load tree
    # -------------------------------------------------------------

    print(
        "Loading phylogenetic tree..."
    )

    tree_leaves, tree = (
        get_tree_leaf_order(
            tree_file
        )
    )

    tree = scale_tree_branches(
        tree,
        scale_factor=args.tree_scale
    )

    # -------------------------------------------------------------
    # Load identity matrix
    # -------------------------------------------------------------

    print(
        "Loading identity matrix..."
    )

    identity_df = pd.read_excel(
        identity_file,
        index_col=0
    )

    print(
        "Identity matrix shape: "
        f"{identity_df.shape}"
    )

    # -------------------------------------------------------------
    # Reorder matrix
    # -------------------------------------------------------------

    print(
        "Reordering identity matrix "
        "according to tree order..."
    )

    (
        reordered_matrix,
        cluster_info,
        available_leaves
    ) = reorder_matrix_by_tree(
        identity_df=identity_df,
        tree_leaves=tree_leaves,
        clusters=clusters,
    )

    # -------------------------------------------------------------
    # Save metadata
    # -------------------------------------------------------------

    print(
        "Saving gene order and metadata..."
    )

    save_gene_order(
        available_leaves,
        output_dir
    )

    cluster_info_df = pd.DataFrame(
        {
            "Gene": available_leaves,
            "Cluster": cluster_info,
        }
    )

    cluster_info_df.to_csv(
        output_dir
        / "cluster_assignment.csv",
        index=False
    )

    reordered_matrix.to_csv(
        output_dir
        / "reordered_identity_matrix.csv"
    )

    # -------------------------------------------------------------
    # Report cluster composition
    # -------------------------------------------------------------

    cluster_counts = (
        pd.Series(
            cluster_info
        )
        .value_counts()
    )

    print(
        "Cluster distribution:"
    )

    for cluster in sorted(
        cluster_counts.index,
        key=lambda x: (
            cluster_sort_key(x)
            if x != "No_Cluster"
            else float("inf")
        )
    ):

        print(
            f"  {cluster}: "
            f"{cluster_counts[cluster]}"
        )

    # -------------------------------------------------------------
    # Create visualization
    # -------------------------------------------------------------

    print(
        "Creating visualization..."
    )

    cluster_color_map = (
        create_combined_plot(
            tree=tree,
            identity_matrix=reordered_matrix,
            cluster_info=cluster_info,
            leaves_order=available_leaves,
            output_dir=output_dir,
            vmin=args.vmin,
            vmax=args.vmax,
        )
    )

    # -------------------------------------------------------------
    # Report cluster colors
    # -------------------------------------------------------------

    if cluster_color_map:

        print(
            "Cluster color mapping:"
        )

        for cluster in sorted(
            cluster_color_map.keys(),
            key=cluster_sort_key
        ):

            print(
                f"  {cluster}: "
                f"{cluster_color_map[cluster]}"
            )

    print(
        "Analysis completed."
    )

    print(
        "Results saved to: "
        f"{output_dir}"
    )


if __name__ == "__main__":
    main()
