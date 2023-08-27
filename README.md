# screfmapping

## Overview

screfmapping is a cell-labeling tool that facilitates the extraction of CD4+ T cells from the single-cell RNA-seq raw data of peripheral blood mononuclear cells (PBMCs) and maps them to our annotated clusterL1, L2 atlas (Link###). The Azimuth pipeline is employed for the extraction of CD4+ T cells, and Symphony mapping, which includes batch effect adjustment by Harmony, is used for mapping to our atlas.

We've included an example analysis in example.R for your reference.
Below, we will provide an overview of the function.

## Function

```
extract_cells_seuratobj(query = q,                   :query_SeuratObject
                        reference = reference,       :Azimuth_reference
                        prefix = prefix)             :output_file_path


reference_mapping_seuratobj(ref = ref,               :our_annotated_clusterL1,L2_data
                            query_obj = query_obj,p  :extracted_CD4T_SeuratObject_with or without_metadata
                            refix = prefix)          :output_file_path
```


## Rscript

```
docker run --rm -it -v ${PWD}:/home/rstudio/autoimmune_10x  yyasumizu/screfmapping:0.0.1 Rscript example.R
```

## Output
### extract_cells_seuratobj()
- ${prefix}_CD4T_MetaData.rds
- ${prefix}_CD4T_AssayData.rds : Input_data_for_ReferenceMapping

### reference_mapping_seuratobj()
- ${prefix}_predict_clusterL1L2_Reference_Mapping.pdf
- ${prefix}_Reference_Mapping.csv : Symphony result
