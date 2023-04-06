# Macro_osteo_comm

## Data Source
- 2022 Nature Communications [Msx1+ stem cells recruited by bioactive tissue engineering graft for bone regeneration](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9445030/)
## Pipeline
1. cellrange
2. Seurat
3. cellchat
## Original function
- cellPathwayHM:To address the limitations of CellChat, this function displays all signals sent or received by a specific cell population and presents them as a heatmap.

## Reproducibility
- 1_preprocess.R: Preprocess  Seurat object with Seurat and Harmony
- 2_cellcomm.R: Infer cell-cell communication with cellchat
