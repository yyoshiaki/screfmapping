# screfmapping


## Overview

screfmapping is a pipeline that facilitates the extraction of CD4+ T cells from the single-cell RNA-seq raw data of peripheral blood mononuclear cells (PBMCs) and maps them to our annotated clusterL1, L2 atlas. The Azimuth pipeline is employed to extract CD4+ T cells, and Symphony mapping, which includes batch effect adjustment by Harmony, is used for mapping to our atlas.

We've included an example analysis in example.R.  
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


## Adjustment
### The number of neighbors (k) to use when finding anchors

<p>
Our "screfmapping" is expected to be used for a PBMC dataset. However, some people may want to use it for a CD4+ T cells enriched dataset. In such cases, we have noticed that a proportion of CD4+ T cells tend to be misannotated as non-CD4+ T cells (approximately 4%). To address this issue, we optimized the `k.anchor` values. In conclusion, the `FindTransferAnchors` in the `extract_cells` function should be conducted with lower `k.anchor` values (for example, `k.anchor = 3`, compared to the default `k.anchor = 5`).<br>
As just a quick note, please modify `ref_mapping_seuratobj.R` if you want to analyse a CD4+ T cells enriched dataset. `k.anchor`'s option will be incorporated in the future revision.
</p>

<pre><code>
# line 40-52 (ref_mapping_seuratobj.R) should be replaced below.
anchors <- FindTransferAnchors(reference = reference$map,
                               query = query,
                               k.anchor = 3,
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
</code></pre>

## Citation

Yasumizu, Yoshiaki, Daiki Takeuchi, Reo Morimoto, Yusuke Takeshima, Tatsusada Okuno, Makoto Kinoshita, Takayoshi Morita, et al. 2023. “Single-Cell Transcriptome Landscape of Circulating CD4+ T Cell Populations in Human Autoimmune Diseases.” bioRxiv. [https://doi.org/10.1101/2023.05.09.540089](https://doi.org/10.1101/2023.05.09.540089).
