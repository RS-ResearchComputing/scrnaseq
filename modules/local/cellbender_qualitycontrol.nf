process CELLBENDER_QC {

    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::scanpy==1.10.2 conda-forge::matplotlib-base conda-forge::numpy"
    container "community.wave.seqera.io/library/scanpy:1.10.2--e83da2205b92a538"

    input:
    tuple val(meta), path(raw_matrix), path(cellbender_matrix)

    output:
    tuple val(meta), path("raw_vs_cellbender_cpg.pdf"), emit: plot
    path  "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3

    import os
    os.environ["NUMBA_CACHE_DIR"] = "."

    import platform
    import numpy as np
    import anndata
    import scanpy as sc
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt


    def format_yaml_like(data: dict, indent: int = 0) -> str:
        yaml_str = ""
        for key, value in data.items():
            spaces = "  " * indent
            if isinstance(value, dict):
                yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
            else:
                yaml_str += f"{spaces}{key}: {value}\\n"
        return yaml_str


    # Read matrices
    adata_raw = sc.read_h5ad("${raw_matrix}")
    adata_cb  = sc.read_h5ad("${cellbender_matrix}")

    # Align to common genes
    common_genes = adata_raw.var_names.intersection(adata_cb.var_names)
    adata_raw = adata_raw[:, common_genes]
    adata_cb  = adata_cb[:, common_genes]

    # Prefer gene_symbols column for labels, fall back to var_names (gene IDs)
    if "gene_symbols" in adata_raw.var.columns:
        gene_labels = adata_raw.var["gene_symbols"].values
    else:
        gene_labels = adata_raw.var_names.values

    # Sum counts per gene across all cells, then log1p-transform to reduce outlier influence
    raw_counts = np.asarray(adata_raw.X.sum(axis=0)).flatten()
    cb_counts  = np.asarray(adata_cb.X.sum(axis=0)).flatten()
    log_raw    = np.log1p(raw_counts)
    log_cb     = np.log1p(cb_counts)

    # Top 20 genes by raw count for labelling (shared across both plots)
    top_idx = np.argsort(raw_counts)[-20:]

    def add_labels(ax, x_vals, y_vals):
        for i in top_idx:
            ax.annotate(
                gene_labels[i],
                xy=(x_vals[i], y_vals[i]),
                xytext=(4, 2),
                textcoords='offset points',
                fontsize=5,
                color='dimgray',
                arrowprops=dict(arrowstyle='-', color='lightgray', lw=0.5),
            )

    from matplotlib.backends.backend_pdf import PdfPages

    with PdfPages("raw_vs_cellbender_cpg.pdf") as pdf:

        # Page 1: log1p-transformed counts
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.scatter(log_raw, log_cb, s=2, alpha=0.3, rasterized=True, color='steelblue')
        max_val = max(log_raw.max(), log_cb.max())
        ax.plot([0, max_val], [0, max_val], '--', color='gray', linewidth=1, label='y = x')
        add_labels(ax, log_raw, log_cb)
        ax.set_xlabel("log1p(Raw counts per gene)")
        ax.set_ylabel("log1p(CellBender counts per gene)")
        ax.set_title("${meta.id}: Raw vs CellBender counts per gene (log1p)")
        ax.legend(loc='upper left')
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # Page 2: raw counts
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.scatter(raw_counts, cb_counts, s=2, alpha=0.3, rasterized=True, color='steelblue')
        max_val = max(raw_counts.max(), cb_counts.max())
        ax.plot([0, max_val], [0, max_val], '--', color='gray', linewidth=1, label='y = x')
        add_labels(ax, raw_counts, cb_counts)
        ax.set_xlabel("Raw counts per gene")
        ax.set_ylabel("CellBender counts per gene")
        ax.set_title("${meta.id}: Raw vs CellBender counts per gene")
        ax.legend(loc='upper left')
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

    # Versions
    versions = {
        "${task.process}": {
            "python":   platform.python_version(),
            "scanpy":   sc.__version__,
            "anndata":  anndata.__version__,
            "numpy":    np.__version__,
            "matplotlib": matplotlib.__version__,
        }
    }
    with open("versions.yml", "w") as f:
        f.write(format_yaml_like(versions))
    """

    stub:
    """
    touch raw_vs_cellbender_cpg.pdf
    touch versions.yml
    """
}
