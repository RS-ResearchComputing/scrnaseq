#!/usr/bin/env python

# Set numba chache dir to current working directory (which is a writable mount also in containers)
import os

os.environ["NUMBA_CACHE_DIR"] = "."

import scanpy as sc
import pandas as pd
from anndata import AnnData, concat as concat_ad
from scipy.sparse import csr_matrix
import platform
import glob
import numpy as np


def _mtx_to_adata(
    matrix: str,
    barcodes: str,
    features: str,
):
    """Load kallisto-formatted mtx files into AnnData"""
    adata = sc.read_mtx(matrix)
    adata.obs_names = pd.read_csv(barcodes, header=None, sep="\\t")[0].values
    adata.var_names = pd.read_csv(features, header=None, sep="\\t")[0].values
    return adata


def _add_metadata(adata: AnnData, t2g: str, sample: str):
    """Add var and obs metadata"""
    adata.obs["sample"] = sample

    txp2gene = pd.read_table(
        t2g, header=None, names=["gene_id", "gene_symbol"], usecols=[1, 2]
    )
    txp2gene = txp2gene.drop_duplicates(subset="gene_id").set_index("gene_id")
    adata.var = adata.var.join(txp2gene, how="left")

    # sanitize gene IDs into standard format
    # index are gene IDs and symbols are a column
    adata.var["gene_versions"] = adata.var.index
    adata.var.index = adata.var["gene_versions"].str.split(".").str[0].values
    adata.var_names_make_unique()  # in case user does not use ensembl references, names might not be unique


def format_yaml_like(data: dict, indent: int = 0) -> str:
    """Formats a dictionary to a YAML-like string.
    Args:
        data (dict): The dictionary to format.
        indent (int): The current indentation level.
    Returns:
        str: A string formatted as YAML.
    """
    yaml_str = ""
    for key, value in data.items():
        spaces = "  " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str


def dump_versions():
    versions = {
        "${task.process}": {
            "python": platform.python_version(),
            "scanpy": sc.__version__,
            "pandas": pd.__version__,
        }
    }

    with open("versions.yml", "w") as f:
        f.write(format_yaml_like(versions))


if __name__ == "__main__":
    # create the directory with the sample name
    os.makedirs("${meta.id}", exist_ok=True)

    # input_type comes from NF module
    if "${params.kb_workflow}" == "standard":
        adata = _mtx_to_adata(
            matrix=glob.glob("${inputs}/*.mtx")[0],
            barcodes=glob.glob("${inputs}/*.barcodes.txt")[0],
            features=glob.glob("${inputs}/*.genes.txt")[0],
        )

    else:
        spliced = _mtx_to_adata(
            matrix=glob.glob("${inputs}/spliced*.mtx")[0],
            barcodes=glob.glob("${inputs}/spliced*.barcodes.txt")[0],
            features=glob.glob("${inputs}/spliced*.genes.txt")[0],
        )
        unspliced = _mtx_to_adata(
            matrix=glob.glob("${inputs}/unspliced*.mtx")[0],
            barcodes=glob.glob("${inputs}/unspliced*.barcodes.txt")[0],
            features=glob.glob("${inputs}/unspliced*.genes.txt")[0],
        )

        # The barcodes of spliced / non-spliced are not necessarily the same.
        # We fill the missing barcodes with zeros
        all_barcodes = list(set(unspliced.obs_names) | set(spliced.obs_names))
        missing_spliced = list(set(unspliced.obs_names) - set(spliced.obs_names))
        missing_unspliced = list(set(spliced.obs_names) - set(unspliced.obs_names))
        ad_missing_spliced = AnnData(
            X=csr_matrix((len(missing_spliced), spliced.shape[1])),
            obs=pd.DataFrame(index=missing_spliced),
            var=spliced.var,
        )
        ad_missing_unspliced = AnnData(
            X=csr_matrix((len(missing_unspliced), spliced.shape[1])),
            obs=pd.DataFrame(index=missing_unspliced),
            var=unspliced.var,
        )

        spliced = concat_ad([spliced, ad_missing_spliced], join="outer")[
            all_barcodes, :
        ]
        unspliced = concat_ad([unspliced, ad_missing_unspliced], join="outer")[
            all_barcodes, :
        ]

        assert np.all(spliced.var_names == unspliced.var_names)

        adata = AnnData(
            X=spliced.X + unspliced.X,
            layers={"unspliced": unspliced.X, "spliced": spliced.X},
            obs=pd.DataFrame(index=all_barcodes),
            var=pd.DataFrame(index=spliced.var_names),
        )

    # out of the conditional: snippet for both standard and non-standard workflows

    # finalize generated adata object
    _add_metadata(adata, t2g="${txp2gene}", sample="${meta.id}")
    adata.write_h5ad("${meta.id}_${meta.input_type}_matrix.h5ad")

    # dump versions
    dump_versions()
