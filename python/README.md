# Python route: mapping onto the atlases with symphonypy

The R pipeline in the repository root maps a query with Seurat + SCTransform + Symphony. This
directory is the Python equivalent, built on [`symphonypy`](https://github.com/potulabe/symphonypy),
for the AnnData/scanpy references of the same atlases.

| | R route (repository root) | Python route (this directory) |
|---|---|---|
| lineage extraction | Azimuth (`extract_cells_seuratobj`) | Azimuth or CellTypist, upstream and up to you |
| normalisation | SCTransform | CP10K + `log1p` + covariate regression |
| reference object | `ref_Reference_Mapping_20220525.RData` (CD4+ T) | `*_symphony_reference_*.h5ad` |
| UMAP | separate `.uwot` cache | stored inside the `h5ad` |
| output | `${prefix}_Reference_Mapping.csv` | labels, `symphony_dist`, UMAP in `query.obs` / `query.obsm` |

The two routes are independent implementations; use one or the other, not a mix.

## The one thing you must not skip

Both atlases were embedded as

```
normalize_total(target_sum=1e4) -> log1p
  -> regress_out([total_counts, pct_counts_mt, S_score, G2M_score]) -> scale -> PCA -> Harmony
```

so `reference.var['mean']` and `reference.var['std']` are statistics of the **regression
residuals**, not of log1p expression. `symphonypy` scales the query with those statistics, so a
query that is only `log1p`-normalised is projected with a constant offset. The query must be
residualised with the reference's own coefficients, stored in `reference.varm['regress_beta']`.
The reference must also carry `reference.uns['harmony']`; without it `symphonypy` falls back to a
soft k-means cache and loses roughly 30 points of label concordance.

Self-mapping of reference cells (kNN-5 `cluster_L2` agreement, 6,000 held-out reference cells
re-mapped as a query, own cell excluded from the neighbourhood):

| query preprocessing | reference | B cells | CD4+ T cells |
|---|---|---|---|
| stored embedding (ceiling) | - | 87.1% | 81.6% |
| `log1p` only | no `uns['harmony']` | 57.7% | 49.1% |
| `log1p` only, CP100K | no `uns['harmony']` | 71.0% | - |
| residuals | no `uns['harmony']` | 55.9% | - |
| **residuals** | **`uns['harmony']`, `key='sample'`** | **84.8%** | **77.5%** |
| residuals | `uns['harmony']`, `key=None` | 81.4% | 72.6% |

This is a property of how the references were normalised, not a change to them: PCA loadings,
Harmony embedding, UMAP and labels are the published ones.

## Files

- `symphonypy_reference_mapping.ipynb` - end-to-end notebook: raw counts to labels, per-cell
  confidence and UMAP, with the sanity checks worth looking at before trusting the labels.
- `screfmapping_symphonypy.py` - `prepare_query()` (the recipe above) and `map_query()`
  (`map_embedding` -> `transfer_labels_kNN` -> `per_cell_confidence` -> `ingest`).
- `data/regev_lab_cell_cycle_genes.txt` - cell-cycle gene list for `S_score` / `G2M_score`
  (first 43 entries are S phase, the remainder G2/M).

## Usage

```python
import scanpy as sc
import screfmapping_symphonypy as srm

reference = sc.read_h5ad("autoimmune_bcell_symphony_reference_v2.h5ad")
query_raw = sc.read_h5ad("query_raw_counts.h5ad")   # raw counts, gene symbols, one lineage only

query = srm.prepare_query(query_raw, reference,
                          cell_cycle_genes="data/regev_lab_cell_cycle_genes.txt",
                          batch_key="sample")
query = srm.map_query(query, reference, key="sample")

query.obs["cluster_L2"]      # transferred label
query.obs["symphony_dist"]   # per-cell confidence, higher = less confident
query.obsm["X_umap"]         # published reference UMAP coordinates
```

## Query input contract

Raw counts in `X`, human gene symbols in `var_names` (duplicate symbols are summed), one row per
cell, already restricted to the lineage of the reference. Extract B cells or CD4+ T cells upstream
with Azimuth or CellTypist; the mapping step does not filter lineages for you. A library/batch
column in `obs` is recommended and passed as `key=`; a single-library query can use `key=None`.

The regression coefficients are only transferable if `total_counts`, `pct_counts_mt`, `S_score`
and `G2M_score` are on the same scale as in the reference. The reference summary is stored in
`reference.uns['symphony_query_preprocessing']['reference_covariate_summary']`, and the notebook
prints both side by side. Very shallow or very deep libraries will be over- or under-corrected;
check `symphony_dist` before trusting the labels.

## References

- CD4+ T cell atlas: Yasumizu et al., *Cell Genomics* (2024),
  https://doi.org/10.1016/j.xgen.2023.100473. The R reference object is on Figshare,
  https://doi.org/10.6084/m9.figshare.25052648. A Symphony-ready AnnData version of this atlas is
  prepared but not yet distributed.
- B cell atlas: manuscript in preparation. The mapping-only AnnData reference is deposited as
  https://doi.org/10.6084/m9.figshare.33292215 and becomes available with that publication.

## Requirements

```
scanpy >= 1.10
symphonypy >= 0.2.2
anndata, numpy, pandas, scipy, scikit-learn, matplotlib
```

Validated with scanpy 1.11.4, symphonypy 0.2.2, anndata 0.12.2, numpy 2.3.5.
