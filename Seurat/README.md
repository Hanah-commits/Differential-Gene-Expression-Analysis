
# Seurat
To demonstrate the DGE analysis of [sc-RNASeq data](https://www.10xgenomics.com/resources/datasets/20-k-mixture-of-nsclc-dt-cs-from-7-donors-3-v-3-1-3-1-standard-6-1-0) of non-small cell lung cancer (NSCLC) dissociated tumor cells obtained from 7 donors.

## Input

Input was obtained as an [HDF5 matrix](https://cf.10xgenomics.com/samples/cell-exp/6.1.2/20k_NSCLC_DTC_3p_nextgem_Multiplex/20k_NSCLC_DTC_3p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5) available on the 10x genomics platform.

## Installations

```
install("Seurat")
install("hdf5r")
```

## Visualization:

### Data Quality Before and After Filtering

![](https://i.imgur.com/hlaNGpt.png)

Positive linear relationship between read depth and gene count.

![](https://i.imgur.com/KczFa9F.png)

Most mitochondrial reads have low sequencing depth.

### Feature Selection

![](https://i.imgur.com/eFo5NQ8.png)

Top 20 highly varible genes are labelled

### Ranking Principle Components using Elbow Plot
![](https://i.imgur.com/Vy3V1gG.png)

Top 20 PCs were chosen.

### Cell Clusters Detected 
![](https://i.imgur.com/FvXHdGZ.png)

