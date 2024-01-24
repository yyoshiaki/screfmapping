# screfmapping


## Overview

"screfmapping" is a pipeline that facilitates the extraction of CD4+ T cells from the single-cell RNA-seq raw data of peripheral blood mononuclear cells (PBMCs) and maps them to our annotated clusterL1, L2 atlas. The Azimuth pipeline is employed to extract CD4+ T cells, and Symphony mapping, which includes batch effect adjustment by Harmony, is used for mapping to our atlas.

We've included an example analysis in `example.R`.  
Below, we provide an overview of the function.


## Function

```
# 1st step, CD4T extraction
extract_cells_seuratobj(query = q,                   # query_SeuratObject
                        reference = reference,       # Azimuth_reference
                        prefix = prefix)             # output_file_path

# 2nd step, label transfer
reference_mapping_seuratobj(ref = ref,               # our_annotated_clusterL1,L2_data
                            query_obj = query_obj,   # extracted_CD4T_SeuratObject_with or without_metadata
                            refix = prefix)          # output_file_path
```


## Our annotated reference data
Please download the `ref_Reference_Mapping_20220525.RData` file from [here](https://doi.org/10.6084/m9.figshare.25052648) and place it in the `/screfmapping/data/` directory.　  
Users will need this file for `ref` in the `reference_mapping_seuratobj` function.


## Rscript

```
docker run --rm -it -v ${PWD}:/home/rstudio/autoimmune_10x  yyasumizu/screfmapping:0.0.1 Rscript example.R
```


  
## Output
### extract_cells_seuratobj
- ${prefix}_CD4T_MetaData.rds
- ${prefix}_CD4T_AssayData.rds : Input_data_for_ReferenceMapping

### reference_mapping_seuratobj
- ${prefix}_predict_clusterL1L2_Reference_Mapping.pdf
- ${prefix}_Reference_Mapping.csv : Symphony result


## Tips
### The number of neighbors (k) to use when finding anchors

Our "screfmapping" is expected to be used for PBMC datasets. However, some people may want to use it for CD4+ T cell-enriched datasets. In such cases, we have noticed that a proportion of CD4+ T cells tend to be misannotated as non-CD4+ T cells (approximately 4%). Empirically, we found that terminally differentiated effector memory T cells (Temra) tended to be annotated as CD8+ T cells or NK cells because those transcriptomes were similar.  
Users may be able to deal with this issue by optimizing the `k.anchor` values of the `FindTransferAnchors` in the `extract_cells_seuratobj` function in `ref_mapping_seuratobj.R`. The lower `k.anchor` values (for example, `k.anchor = 3`, compared to the default `k.anchor = 5`) worked well for CD4+ T cell-enriched datasets.

```
# lines 40-52 (in `ref_mapping_seuratobj.R`) should be replaced as below.

anchors <- FindTransferAnchors(reference = reference$map,
                               query = query,
                               k.anchor = 3, # change here
                               k.filter = NA,
                               reference.neighbors = "refdr.annoy.neighbors",
                               reference.assay = "refAssay",
                               query.assay = "refAssay",
                               reference.reduction = "refDR",
                               normalization.method = "SCT",
                               features = intersect(rownames(x = reference$map),
                                                    VariableFeatures(object = query)),
                               dims = 1:50,
                               n.trees = 20,
                               mapping.score.k = 100)
```

## Citation

Yasumizu, Y., Takeuchi, D., Morimoto, R., Takeshima, Y., Okuno, T., Kinoshita, M., Morita, T., Kato, Y., Wang, M., Motooka, D., et al. (2024). Single-cell transcriptome landscape of circulating CD4+ T cell populations in autoimmune diseases. Cell Genomics.  
[https://doi.org/10.1016/j.xgen.2023.100473](https://doi.org/10.1016/j.xgen.2023.100473)

