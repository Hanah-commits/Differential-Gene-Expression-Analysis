# Limma

To demonstrate the DGE analysis of two samples with multiple replicates. 

## Input matrix


The count matrix used was obtained using featureCounts. Each replicate is labelled with its sample type and the ENCODE identifier.


```
             GENE ectodermalcell_ENCFF148VCG.bam ectodermalcell_ENCFF154POQ.bam ectodermalcell_ENCFF370OWM.bam
ENSG00000227232.5                             30                             39                             71
ENSG00000233750.3                              2                              5                              3
ENSG00000268903.1                             19                             31                             37
ENSG00000269981.1                             15                             27                             27
ENSG00000241860.7                             13                             25                             21
ENSG00000279457.4                             75                            111                            261
```

Raw counts of genes in 3 replicates of the ectodermalcell sample

## Normalization

The reads of the differentially expressed gene take up most of the library. This results in RNA compositional bias. edgeR and limma use *effective library size* to deal with compositional bias and difference in libraray sizes between samples. It is obtained as a product of library size and the normalization factor. The raw counts of the genes are normalized by dividing the count for each sample by its effective library size.

```r
calcNormFactors(object, method=c("TMM","RLE","upperquartile","none"))
```
The normalization factor can be any one of the following:
* TMM / trimmed mean of M-values (default)
* RLE / relative log expression
* upperquartile
* none: normalization factors are set to 1

## Mean-Variance Trend


![](https://i.imgur.com/V8KiF60.png)

Variance decreases as the expression level increases. This can plot can be used to determine whether further processing of data is required before DGE analysis.



