"""Query preprocessing and Symphony mapping for the Python (symphonypy) references of
Yasumizu et al. autoimmune single-cell atlases (CD4+ T cells and B cells).

Why this module exists
----------------------
Both reference atlases were built with

    normalize_total(target_sum=1e4) -> log1p
      -> regress_out([total_counts, pct_counts_mt, S_score, G2M_score]) -> scale -> PCA -> Harmony

so the gene statistics stored in ``reference.var['mean']`` and ``reference.var['std']`` are the
statistics of the *regression residuals*, not of log1p expression. ``symphonypy`` scales the query
with those statistics, so handing it plain log1p counts shifts the projection by a constant and
degrades label transfer (self-mapping of reference cells: 58% vs an 87% ceiling for B cells,
49% vs 82% for CD4+ T cells).

The reference objects distributed for the Python route therefore carry

* ``varm['regress_beta']`` - per-gene regression coefficients ``(n_genes, 5)``
  in the column order ``[intercept, total_counts, pct_counts_mt, S_score, G2M_score]``, and
* ``uns['harmony']`` - the Symphony compression object required to reproduce the per-batch
  Harmony correction (without it ``symphonypy`` silently falls back to a soft k-means cache),
* ``uns['symphony_query_preprocessing']`` - a machine-readable copy of the recipe below.

``prepare_query`` applies exactly the reference recipe to a query and returns the residual matrix
that ``symphonypy.tl.map_embedding`` expects. Reference PCA loadings, Harmony embedding, UMAP and
labels are never modified.

Recipe
------
1. raw counts, human gene symbols (duplicate symbols summed), a batch/library column
2. ``calculate_qc_metrics(qc_vars=['mt'])`` on **all** genes -> ``total_counts``, ``pct_counts_mt``
3. ``normalize_total(target_sum=1e4)`` -> ``log1p``   (1e4, not 1e5)
4. drop genes starting with ``IGKV`` / ``IGLV`` / ``IGHV`` / ``IGLC``
5. ``score_genes_cell_cycle`` with the Regev-lab list (first 43 entries = S, remainder = G2M)
6. for the reference genes: ``resid = log1p - [1, total_counts, pct_counts_mt, S_score, G2M_score] @ regress_beta.T``
   (genes absent from the query are left at 0, i.e. the reference mean)
7. ``map_embedding`` -> ``transfer_labels_kNN`` -> ``per_cell_confidence`` -> ``ingest``

Notes
-----
* Steps 2-5 must run in this order: ``total_counts`` and ``pct_counts_mt`` are computed on the full
  gene set *before* the immunoglobulin V/C genes are dropped, which is how the reference was built.
* ``key`` in ``map_query`` should name a query batch column (``key='sample'``) when the query has
  several libraries. ``key=None`` skips the Harmony correction of the query and costs a few points
  of label concordance (B cells: 85% -> 81%; CD4+ T cells: 77% -> 73%).
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sps

COVARIATES = ["total_counts", "pct_counts_mt", "S_score", "G2M_score"]
IG_PREFIXES = ("IGKV", "IGLV", "IGHV", "IGLC")
N_S_PHASE_GENES = 43


def load_cell_cycle_genes(path: str | Path) -> tuple[list[str], list[str]]:
    """Read the Regev-lab cell-cycle gene list (first 43 = S phase, remainder = G2/M)."""
    genes = [line.strip() for line in Path(path).read_text().splitlines() if line.strip()]
    return genes[:N_S_PHASE_GENES], genes[N_S_PHASE_GENES:]


def collapse_to_gene_symbols(adata: ad.AnnData, symbol_column: str | None = None) -> ad.AnnData:
    """Sum counts of duplicate gene symbols so that ``var_names`` are unique symbols."""
    if symbol_column is not None and symbol_column in adata.var.columns:
        symbols = adata.var[symbol_column].astype(str)
        symbols = symbols.where(
            symbols.notna() & symbols.ne("") & symbols.ne("nan"),
            pd.Series(adata.var_names, index=adata.var.index),
        )
    else:
        symbols = pd.Series(adata.var_names.astype(str), index=adata.var.index)

    codes, unique = pd.factorize(symbols, sort=False)
    if len(unique) == adata.n_vars:
        out = adata.copy()
        out.var_names = pd.Index(symbols.to_numpy(), name=adata.var.index.name)
        return out

    collapse = sps.csr_matrix(
        (np.ones(adata.n_vars), (np.arange(adata.n_vars), codes)),
        shape=(adata.n_vars, len(unique)),
    )
    counts = sps.csr_matrix(adata.X) @ collapse
    return ad.AnnData(
        counts.astype(np.float32),
        obs=adata.obs.copy(),
        var=pd.DataFrame(index=pd.Index(np.asarray(unique, dtype=str), name="gene_symbol")),
    )


def prepare_query(
    query: ad.AnnData,
    reference: ad.AnnData,
    cell_cycle_genes: str | Path,
    batch_key: str | None = "sample",
    symbol_column: str | None = None,
    verbose: bool = True,
) -> ad.AnnData:
    """Preprocess a raw-count query exactly like the reference atlas.

    Parameters
    ----------
    query
        Raw counts in ``X`` (integers), human gene symbols in ``var_names``, one cell per row,
        already restricted to the lineage of the reference (B cells or CD4+ T cells).
    reference
        Reference AnnData carrying ``varm['regress_beta']`` (and ``uns['harmony']`` for mapping).
    cell_cycle_genes
        Path to ``regev_lab_cell_cycle_genes.txt``.
    batch_key
        Column in ``query.obs`` identifying library/batch; copied to ``obs['sample']``.

    Returns
    -------
    AnnData with the residualised matrix on the reference gene set, ready for
    ``symphonypy.tl.map_embedding``.
    """
    import scanpy as sc

    if "regress_beta" not in reference.varm:
        raise KeyError(
            "reference lacks varm['regress_beta']; this query recipe only applies to the "
            "regress_out-based Python references (see module docstring)"
        )

    q = collapse_to_gene_symbols(query, symbol_column=symbol_column)
    q.var_names_make_unique()

    # Covariates are computed on the full gene set, before dropping immunoglobulin V/C genes.
    q.var["mt"] = q.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(q, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    sc.pp.normalize_total(q, target_sum=1e4)
    sc.pp.log1p(q)

    keep = ~np.logical_or.reduce([q.var_names.str.startswith(p) for p in IG_PREFIXES])
    q = q[:, keep].copy()

    s_genes, g2m_genes = load_cell_cycle_genes(cell_cycle_genes)
    sc.tl.score_genes_cell_cycle(
        q,
        s_genes=[g for g in s_genes if g in q.var_names],
        g2m_genes=[g for g in g2m_genes if g in q.var_names],
    )

    ref_genes = pd.Index(reference.var_names)
    shared = ref_genes.intersection(q.var_names)
    missing = ref_genes.difference(q.var_names)
    if verbose:
        print(
            f"reference genes: {len(ref_genes)} | found in query: {len(shared)} | "
            f"missing (set to the reference mean): {len(missing)}"
        )
        if len(missing):
            print("  first missing:", list(missing[:10]))

    x = np.zeros((q.n_obs, len(ref_genes)), dtype=np.float64)
    positions = ref_genes.get_indexer(shared)
    block = q[:, shared].X
    x[:, positions] = block.toarray() if sps.issparse(block) else np.asarray(block)

    design = np.column_stack(
        [np.ones(q.n_obs), q.obs[COVARIATES].to_numpy(dtype=np.float64)]
    )
    beta = np.asarray(reference.varm["regress_beta"], dtype=np.float64)  # (n_genes, 5)
    resid = x - design @ beta.T
    if len(missing):
        resid[:, ref_genes.get_indexer(missing)] = 0.0

    out = ad.AnnData(
        resid.astype(np.float32),
        obs=q.obs.copy(),
        var=pd.DataFrame(index=ref_genes.copy()),
    )
    out.uns["log1p"] = {"base": None}
    out.uns["screfmapping_query_preprocessing"] = {
        "version": "260918",
        "target_sum": 1e4,
        "covariates": COVARIATES,
        "n_reference_genes": int(len(ref_genes)),
        "n_genes_found": int(len(shared)),
        "n_genes_missing": int(len(missing)),
    }
    if batch_key is not None:
        if batch_key not in out.obs.columns:
            raise KeyError(f"batch_key {batch_key!r} not found in query.obs")
        out.obs["sample"] = out.obs[batch_key].astype(str)
    return out


def map_query(
    query_resid: ad.AnnData,
    reference: ad.AnnData,
    key: str | None = "sample",
    labels: tuple[str, ...] = ("cluster_L1", "cluster_L2"),
    n_neighbors: int = 5,
    ingest: bool = True,
) -> ad.AnnData:
    """Map a residualised query onto the reference and transfer labels.

    Adds ``obsm['X_pca_reference']``, ``obsm['X_pca_harmony']``, the transferred label columns,
    ``obs['symphony_dist']`` (per-cell confidence; higher = less confident) and, when
    ``ingest=True``, ``obsm['X_umap']`` in the published reference UMAP space.
    """
    import symphonypy as sp

    if "harmony" not in reference.uns:
        raise KeyError(
            "reference lacks uns['harmony']; symphonypy would silently fall back to a soft "
            "k-means cache and lose ~30 points of label concordance"
        )
    present = [lab for lab in labels if lab in reference.obs.columns]
    if not present:
        raise KeyError(f"none of {labels} found in reference.obs")

    sp.tl.map_embedding(query_resid, reference, key=key)
    sp.tl.transfer_labels_kNN(
        query_resid,
        reference,
        ref_labels=present,
        query_labels=present,
        n_neighbors=n_neighbors,
        weights="distance",
    )
    sp.tl.per_cell_confidence(query_resid, reference, obs="symphony_dist")
    if ingest:
        sp.tl.ingest(query_resid, reference)
    return query_resid
