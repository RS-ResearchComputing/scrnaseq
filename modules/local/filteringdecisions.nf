process FILTERING_DECISIONS {

    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::scanpy==1.10.2 conda-forge::seaborn conda-forge::matplotlib-base conda-forge::numpy"
    container "community.wave.seqera.io/library/scanpy:1.10.2--e83da2205b92a538"

    input:
    tuple val(meta), path(raw_h5), path(filtered_h5)

    output:
    tuple val(meta), path("violin_plots_cell_gene.png"), emit: plot
    path  "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3

    import os
    os.environ["NUMBA_CACHE_DIR"] = "."

    import platform
    import warnings
    import numpy as np
    import anndata
    import scanpy as sc
    import seaborn as sns
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


    def load_and_compute_metrics(path):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            adata = sc.read_10x_h5(path)
        duplicates = adata.var_names.duplicated(keep='first')
        adata = adata[:, ~duplicates].copy()
        mito_genes = [g for g in adata.var_names if g.startswith('MT-') or g.startswith('mt-')]
        adata.obs['pct_counts_mt'] = (
            np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
        )
        adata.obs['n_counts']   = adata.X.sum(axis=1).A1
        adata.obs['log_counts'] = np.log(adata.obs['n_counts'])
        adata.obs['n_genes']    = (adata.X > 0).sum(1)
        adata.var['n_cells']      = np.array((adata.X > 0).sum(axis=0)).ravel()
        adata.var['total_counts'] = np.array(adata.X.sum(axis=0)).ravel()
        return adata


    matrices = [
        ("${raw_h5}",      "Raw"),
        ("${filtered_h5}", "Filtered"),
    ]

    col_titles = [
        'log_counts_per_cell',
        'n_genes_per_cell',
        'pct_counts_mt',
        'log1p_cells_per_gene',
        'log1p_counts_per_gene',
    ]

    fig, axes = plt.subplots(2, 5, figsize=(15, 8))
    fig.suptitle("${meta.id}", fontsize=14, y=1.01)

    for col, title in enumerate(col_titles):
        axes[0, col].set_title(title, fontsize=9)

    import matplotlib.collections as mc

    for row, (path, label) in enumerate(matrices):
        adata = load_and_compute_metrics(path)

        for col_idx, key in enumerate(['log_counts', 'n_genes', 'pct_counts_mt']):
            ax = axes[row, col_idx]
            sc.pl.violin(adata, [key], ax=ax, show=False, ylabel='', jitter=0.15, size=0.8)
            for coll in ax.collections:
                if isinstance(coll, mc.PathCollection):
                    coll.set_alpha(0.2)
                    coll.set_linewidths(0)

        sns.violinplot(y=np.log1p(adata.var['n_cells']),      inner='box', ax=axes[row, 3])
        sns.violinplot(y=np.log1p(adata.var['total_counts']), inner='box', ax=axes[row, 4])

        for col in range(5):
            ax = axes[row, col]
            ax.set_xticks([])
            ax.set_xlabel('')
            ax.set_facecolor('#e8f4fa' if col < 3 else '#fff4e6')
            ax.yaxis.grid(True, color='lightgray', linewidth=0.5, zorder=0)
            ax.set_axisbelow(True)

        # Set after violin calls so sc.pl.violin ylabel='' does not overwrite it
        axes[row, 0].set_ylabel(label, fontsize=10, labelpad=8)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig("violin_plots_cell_gene.png", dpi=150, bbox_inches='tight')
    plt.close()

    # Versions
    versions = {
        "${task.process}": {
            "python":     platform.python_version(),
            "scanpy":     sc.__version__,
            "anndata":    anndata.__version__,
            "numpy":      np.__version__,
            "seaborn":    sns.__version__,
            "matplotlib": matplotlib.__version__,
        }
    }
    with open("versions.yml", "w") as f:
        f.write(format_yaml_like(versions))
    """

    stub:
    """
    touch violin_plots_cell_gene.png
    touch versions.yml
    """
}
