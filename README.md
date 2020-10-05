# LiverAnnotation

## Reference Data
`Halpern_Human_Ortho_sig_profiles.rds`
Source: SupplementaryTable 3 https://www.nature.com/articles/nature21065
Filters: 5% FDR, has human ortholog with HDNC gene symbol.
Gene ids: Human gene symbols - duplicates from many-one relationships reduced to the gene with highest spatial patterning.

`Halpern_Mouse_sig_profiles.rds`
Source: Halpern at al. (2017) SupplementaryTable 3 https://www.nature.com/articles/nature21065
Filters: 5% FDR
Gene ids: Mouse MGI symbols

`Zheng_immune_profiles.rds`
Source: Zheng et al. (2017) https://support.10xgenomics.com/single-cell-gene-expression/datasets
Filters: 500 genes/cell, 750 UMI/cell, sorted immune cells only.
Processing: Normalized to 15,000 umi/cell, log2 transformed, pseudocount of 1.
Specificity calculation: 
0 if max mean expression < 0.5
(max(mean expression) - min(mean expression))/max(0.01, min(mean expression))

`Map2_scmap_minimal_reference.rds`
Source: MacParland et al. (2018) https://www.nature.com/articles/s41467-018-06318-7
Filters & Preprocessing: Same as original publication.

`LiverMap1_Markers.txt`
Source: MacParland et al. (2018) https://www.nature.com/articles/s41467-018-06318-7
Definition: 
Cut-point chosen to maximized the difference between "on" and "off" cell-types
This difference must be at least 0.3 difference in log mean expression and 0.1 difference in detection rate.
Markers are classified by specific or general cell-type according to which clusters were "on" and "off" to maximize the difference.
Markers unique to a cell-type had to be so both when considering mean expression and detection rate.





## Package Prequisites
SingleCellExperiment, scmap, SCINA, Seurat, proxy, gplots, ggplot2, RColourBrewer

## Scripts
`Liver_Colour_Scheme.R` : Matches standard cell types to a set of colours for consistency across plots.
`Autoannotation.R` : Wrappers for scmap & SCINA & a custom (hypergeometric test) autoannotation tool.
`Manual_annotation.R` : Script for making plots for manual annotation.
`Generic_functions.R` : Simple functions for calculating mean expression in across clusters.


