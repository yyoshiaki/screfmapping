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

## Preprocess the query the way the reference was built

Symphony projects a query by scaling it with the reference's per-gene statistics and applying the
reference's PCA loadings, so the query has to go through the same transformations as the
reference. These atlases were built as

```
normalize_total(target_sum=1e4) -> log1p
  -> regress_out([total_counts, pct_counts_mt, S_score, G2M_score]) -> scale -> PCA -> Harmony
```

meaning `reference.var['mean']` and `['std']` describe the data **after** that covariate
regression. A query that stops at `log1p` sits on a different scale, and the projection carries
the difference through to the transferred labels.

The references therefore carry the coefficients they were regressed with
(`varm['regress_beta']`), the Symphony compression object (`uns['harmony']`) and the recipe
itself (`uns['symphony_query_preprocessing']`). `prepare_query()` replays all of it on your data.
The notebook runs the whole flow on a public dataset and ends by showing where the same cells
land when this step is skipped.

None of this changes the references: PCA loadings, Harmony embedding, UMAP and labels are the
published ones.

## Files

- `symphonypy_reference_mapping.ipynb` - end-to-end notebook on 10x PBMC 3k: query preprocessing,
  covariate check, mapping, labels on the atlas UMAP, marker genes, and the effect of skipping
  the covariate regression.
- `screfmapping_symphonypy.py` - `prepare_query()` (the recipe above) and `map_query()`
  (`map_embedding` -> `transfer_labels_kNN` -> `per_cell_confidence` -> `ingest`).
- `data/regev_lab_cell_cycle_genes.txt` - cell-cycle gene list for `S_score` / `G2M_score` (first
  43 entries are S phase, the remainder G2/M). From Tirosh et al., *Science* (2016), in the form
  distributed with the scanpy cell-cycle tutorial; vendored here because that upstream copy is no
  longer online.

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
column in `obs` is recommended and passed as `key=`, which lets Symphony correct the query batches
against the reference; a single-library query uses `key=None`.

The regression coefficients transfer cleanly when `total_counts`, `pct_counts_mt`, `S_score` and
`G2M_score` sit in a range similar to the reference. The reference summary is stored in
`reference.uns['symphony_query_preprocessing']['reference_covariate_summary']` and the notebook
prints both side by side. Reference genes missing from the query are held at the reference mean,
so a query with fewer genes still maps; `prepare_query()` reports how many were found.

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

Validated with scanpy 1.11.4, symphonypy 0.2.2, anndata 0.12.2.
