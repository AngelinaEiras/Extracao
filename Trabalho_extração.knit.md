
<!-- rnb-text-begin -->

---
title: "Grupo 6 - Trabalho prático de Extração de Conhecimento de Bases de Dados Biológicas"
author: Angelina Eiras - PG42861, Beatriz Silva - PG, Mariana Gonçalves - PG, Quitéria Pinheiro - PG
date: 30/05/22
output: html_notebook
# output: html_document
---


<!-- rnb-text-end -->



<!-- rnb-text-begin -->


# Explicação dos dados, a sua origem e relevância

...


## Obtenção dos dados


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuR0RTNTgyNiA8LSBnZXRHRU8oJ0dEUzU4MjYnLCBkZXN0ZGlyPVwiLlwiKVxuYGBgIn0= -->

```r
GDS5826 <- getGEO('GDS5826', destdir=".")
```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiVXNpbmcgbG9jYWxseSBjYWNoZWQgdmVyc2lvbiBvZiBHRFM1ODI2IGZvdW5kIGhlcmU6XG4uL0dEUzU4MjYuc29mdC5neiBcbiJ9 -->

```
Using locally cached version of GDS5826 found here:
./GDS5826.soft.gz 
```



<!-- rnb-output-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


## Obtenção de packages do BioConductor


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjpbImlmICghcmVxdWlyZU5hbWVzcGFjZShcIkJpb2NNYW5hZ2VyXCIsIHF1aWV0bHkgPSBUUlVFKSApIiwiICBpbnN0YWxsLnBhY2thZ2VzKFwiQmlvY01hbmFnZXJcIikiLCIgIEJpb2NNYW5hZ2VyOjppbnN0YWxsKCkiLCJCaW9jTWFuYWdlcjo6aW5zdGFsbChcIkJpb2Jhc2VcIiwgZm9yY2UgPSBUUlVFKSIsIkJpb2NNYW5hZ2VyOjppbnN0YWxsKFwiZ2VuZWZpbHRlclwiLCBmb3JjZSA9IFRSVUUpIiwiQmlvY01hbmFnZXI6Omluc3RhbGwoXCJvcmcuSHMuZWcuZGJcIiwgZm9yY2UgPSBUUlVFKSIsIkJpb2NNYW5hZ2VyOjppbnN0YWxsKFwiREVTZXEyXCIsIGZvcmNlID0gVFJVRSkiLCJCaW9jTWFuYWdlcjo6aW5zdGFsbChcImZnc2VhXCIsIGZvcmNlID0gVFJVRSkiLCJCaW9jTWFuYWdlcjo6aW5zdGFsbChcImxpbW1hXCIsIGZvcmNlID0gVFJVRSkiLCJCaW9jTWFuYWdlcjo6aW5zdGFsbChcIkFMTFwiLCBmb3JjZSA9IFRSVUUpIl19 -->

```r
if (!requireNamespace("BiocManager", quietly = TRUE) )
  install.packages("BiocManager")
  BiocManager::install()
BiocManager::install("Biobase", force = TRUE)
BiocManager::install("genefilter", force = TRUE)
BiocManager::install("org.Hs.eg.db", force = TRUE)
BiocManager::install("DESeq2", force = TRUE)
BiocManager::install("fgsea", force = TRUE)
BiocManager::install("limma", force = TRUE)
BiocManager::install("ALL", force = TRUE)
```



<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


## Carregamento de packages


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxubGlicmFyeShjYXIpXG5saWJyYXJ5KHJlYWRyKVxubGlicmFyeShzdW1tYXJ5dG9vbHMpXG5saWJyYXJ5KGxpbW1hKVxubGlicmFyeShCaW9iYXNlKVxubGlicmFyeShnZW5lZmlsdGVyKVxubGlicmFyeShkcGx5cilcbmxpYnJhcnkob3JnLkhzLmVnLmRiKVxubGlicmFyeShERVNlcTIpXG5saWJyYXJ5KHBoZWF0bWFwKVxubGlicmFyeShmZ3NlYSlcbmxpYnJhcnkoZ2dwbG90MilcbmxpYnJhcnkoQUxMKVxubGlicmFyeShSdHNuZSlcbmxpYnJhcnkoY2xhc3MpXG5saWJyYXJ5KGUxMDcxKVxubGlicmFyeShwYXJ0eSlcbmxpYnJhcnkobm5ldClcbmxpYnJhcnkocmFuZG9tRm9yZXN0KVxubGlicmFyeShkYXRhLnRhYmxlKVxubGlicmFyeShjYXJldClcbmBgYCJ9 -->

```r
library(car)
library(readr)
library(summarytools)
library(limma)
library(Biobase)
library(genefilter)
library(dplyr)
library(org.Hs.eg.db)
library(DESeq2)
library(pheatmap)
library(fgsea)
library(ggplot2)
library(ALL)
library(Rtsne)
library(class)
library(e1071)
library(party)
library(nnet)
library(randomForest)
library(data.table)
library(caret)
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



# Pré-processamento dos dados

...


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuc3VtKGlzLm5hKEdEUzU4MjYpKSAjIHNvbWEgbyBudW1lcm8gZGUgTkFzXG5gYGAifQ== -->

```r
sum(is.na(GDS5826)) # soma o numero de NAs
```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiV2FybmluZyBpbiBpcy5uYShHRFM1ODI2KSA6XG4gIGlzLm5hKCkgYXBwbGllZCB0byBub24tKGxpc3Qgb3IgdmVjdG9yKSBvZiB0eXBlICdTNCdcbiJ9 -->

```
Warning in is.na(GDS5826) :
  is.na() applied to non-(list or vector) of type 'S4'
```



<!-- rnb-output-end -->

<!-- rnb-output-begin eyJkYXRhIjoiWzFdIDBcbiJ9 -->

```
[1] 0
```



<!-- rnb-output-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuIyByIDwtIGRhdGEuZnJhbWUoc2FwcGx5KEdEUzU4MjYsIGZ1bmN0aW9uKHgpIHN1bShpcy5uYSh4KSkpKVxuIyBzdW0oISFyKSAjIENvbHVuYXMgY29tIE5Bc1xuIyBuw6NvIHRlbSBOQXNcblxuXG5kYXRhID0gR0RTMmVTZXQoR0RTNTgyNilcbmBgYCJ9 -->

```r
# r <- data.frame(sapply(GDS5826, function(x) sum(is.na(x))))
# sum(!!r) # Colunas com NAs
# não tem NAs


data = GDS2eSet(GDS5826)
```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiVXNpbmcgbG9jYWxseSBjYWNoZWQgdmVyc2lvbiBvZiBHUEw1NzAgZm91bmQgaGVyZTpcbi90bXAvUnRtcEhTMUpMTy9HUEw1NzAuYW5ub3QuZ3ogXG4ifQ== -->

```
Using locally cached version of GPL570 found here:
/tmp/RtmpHS1JLO/GPL570.annot.gz 
```



<!-- rnb-output-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuZGltKGRhdGEpXG5gYGAifQ== -->

```r
dim(data)
```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiRmVhdHVyZXMgIFNhbXBsZXMgXG4gICA1NDY3NSAgICAgICAxMiBcbiJ9 -->

```
Features  Samples 
   54675       12 
```



<!-- rnb-output-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuY2xhc3MoZGF0YSlcbmBgYCJ9 -->

```r
class(data)
```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiWzFdIFwiRXhwcmVzc2lvblNldFwiXG5hdHRyKCxcInBhY2thZ2VcIilcblsxXSBcIkJpb2Jhc2VcIlxuIn0= -->

```
[1] "ExpressionSet"
attr(,"package")
[1] "Biobase"
```



<!-- rnb-output-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuTWV0YShHRFM1ODI2KVxuYGBgIn0= -->

```r
Meta(GDS5826)
```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiJGNoYW5uZWxfY291bnRcblsxXSBcIjFcIlxuXG4kZGF0YXNldF9pZFxuWzFdIFwiR0RTNTgyNlwiIFwiR0RTNTgyNlwiIFwiR0RTNTgyNlwiIFwiR0RTNTgyNlwiIFwiR0RTNTgyNlwiIFwiR0RTNTgyNlwiXG5cbiRkZXNjcmlwdGlvblxuWzFdIFwiQW5hbHlzaXMgb2YgcHJvdGVhc29tZSBpbmhpYml0b3IgY2FyZmlsem9taWItcmVzaXN0YW50IG11bHRpcGxlIG15ZWxvbWEgKE1NKSBjZWxsIGxpbmVzIEtNUy0xMS9DZnogYW5kIEtNUy0zNC9DZnosIGFmdGVyIDEgd2VlayBvZiBncm93dGggaW4gdGhlIGFic2VuY2Ugb2YgY2FyZmlsem9taWIuIFJlc3VsdHMgcHJvdmlkZSBpbnNpZ2h0IGludG8gdGhlIG1vbGVjdWxhciBtZWNoYW5pc21zIHVuZGVybHlpbmcgdGhlIGFjcXVpc2l0aW9uIG9mIHByb3RlYXNvbWUgaW5oaWJpdG9yIHJlc2lzdGFuY2UgaW4gTU0uXCJcblsyXSBcIktNUy0xMS9DZnpcIiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXG5bM10gXCJLTVMtMzQvQ2Z6XCIgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFxuWzRdIFwiS01TLTExXCIgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcbls1XSBcIktNUy0zNFwiICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXG5bNl0gXCJjYXJmaWx6b21pYi1yZXNpc3RhbnQgTU1cIiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFxuWzddIFwicGFyZW50YWwgTU1cIiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcblxuJGVtYWlsXG5bMV0gXCJnZW9AbmNiaS5ubG0ubmloLmdvdlwiXG5cbiRmZWF0dXJlX2NvdW50XG5bMV0gXCI1NDY3NVwiXG5cbiRpbnN0aXR1dGVcblsxXSBcIk5DQkkgTkxNIE5JSFwiXG5cbiRuYW1lXG5bMV0gXCJHZW5lIEV4cHJlc3Npb24gT21uaWJ1cyAoR0VPKVwiXG5cbiRvcmRlclxuWzFdIFwibm9uZVwiXG5cbiRwbGF0Zm9ybVxuWzFdIFwiR1BMNTcwXCJcblxuJHBsYXRmb3JtX29yZ2FuaXNtXG5bMV0gXCJIb21vIHNhcGllbnNcIlxuXG4kcGxhdGZvcm1fdGVjaG5vbG9neV90eXBlXG5bMV0gXCJpbiBzaXR1IG9saWdvbnVjbGVvdGlkZVwiXG5cbiRwdWJtZWRfaWRcblsxXSBcIjI2MTA5NDMzXCJcblxuJHJlZlxuWzFdIFwiTnVjbGVpYyBBY2lkcyBSZXMuIDIwMDUgSmFuIDE7MzMgRGF0YWJhc2UgSXNzdWU6RDU2Mi02XCJcblxuJHJlZmVyZW5jZV9zZXJpZXNcblsxXSBcIkdTRTY5MDc4XCJcblxuJHNhbXBsZV9jb3VudFxuWzFdIFwiMTJcIlxuXG4kc2FtcGxlX2lkXG5bMV0gXCJHU00xNjkyNTg3LEdTTTE2OTI1ODgsR1NNMTY5MjU4OVwiICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXG5bMl0gXCJHU00xNjkyNTkwLEdTTTE2OTI1OTEsR1NNMTY5MjU5MlwiICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXG5bM10gXCJHU00xNjkyNTkzLEdTTTE2OTI1OTQsR1NNMTY5MjU5NVwiICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXG5bNF0gXCJHU00xNjkyNTk2LEdTTTE2OTI1OTcsR1NNMTY5MjU5OFwiICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXG5bNV0gXCJHU00xNjkyNTg3LEdTTTE2OTI1ODgsR1NNMTY5MjU4OSxHU00xNjkyNTkwLEdTTTE2OTI1OTEsR1NNMTY5MjU5MlwiXG5bNl0gXCJHU00xNjkyNTkzLEdTTTE2OTI1OTQsR1NNMTY5MjU5NSxHU00xNjkyNTk2LEdTTTE2OTI1OTcsR1NNMTY5MjU5OFwiXG5cbiRzYW1wbGVfb3JnYW5pc21cblsxXSBcIkhvbW8gc2FwaWVuc1wiXG5cbiRzYW1wbGVfdHlwZVxuWzFdIFwiUk5BXCJcblxuJHRpdGxlXG5bMV0gXCJNdWx0aXBsZSBteWVsb21hIGNlbGwgbGluZXMgd2l0aCBhY3F1aXJlZCByZXNpc3RhbmNlIHRvIGNoZW1vdGhlcmFwZXV0aWMgYWdlbnQgY2FyZmlsem9taWJcIlxuXG4kdHlwZVxuWzFdIFwiRXhwcmVzc2lvbiBwcm9maWxpbmcgYnkgYXJyYXlcIiBcImNlbGwgbGluZVwiICAgICAgICAgICAgICAgICAgICAgXCJjZWxsIGxpbmVcIiAgICAgICAgICAgICAgICAgICAgXG5bNF0gXCJjZWxsIGxpbmVcIiAgICAgICAgICAgICAgICAgICAgIFwiY2VsbCBsaW5lXCIgICAgICAgICAgICAgICAgICAgICBcImNlbGwgdHlwZVwiICAgICAgICAgICAgICAgICAgICBcbls3XSBcImNlbGwgdHlwZVwiICAgICAgICAgICAgICAgICAgICBcblxuJHVwZGF0ZV9kYXRlXG5bMV0gXCJKdWwgMTEgMjAxNlwiXG5cbiR2YWx1ZV90eXBlXG5bMV0gXCJ0cmFuc2Zvcm1lZCBjb3VudFwiXG5cbiR3ZWJfbGlua1xuWzFdIFwiaHR0cDovL3d3dy5uY2JpLm5sbS5uaWguZ292L2dlb1wiXG4ifQ== -->

```
$channel_count
[1] "1"

$dataset_id
[1] "GDS5826" "GDS5826" "GDS5826" "GDS5826" "GDS5826" "GDS5826"

$description
[1] "Analysis of proteasome inhibitor carfilzomib-resistant multiple myeloma (MM) cell lines KMS-11/Cfz and KMS-34/Cfz, after 1 week of growth in the absence of carfilzomib. Results provide insight into the molecular mechanisms underlying the acquisition of proteasome inhibitor resistance in MM."
[2] "KMS-11/Cfz"                                                                                                                                                                                                                                                                                         
[3] "KMS-34/Cfz"                                                                                                                                                                                                                                                                                         
[4] "KMS-11"                                                                                                                                                                                                                                                                                             
[5] "KMS-34"                                                                                                                                                                                                                                                                                             
[6] "carfilzomib-resistant MM"                                                                                                                                                                                                                                                                           
[7] "parental MM"                                                                                                                                                                                                                                                                                        

$email
[1] "geo@ncbi.nlm.nih.gov"

$feature_count
[1] "54675"

$institute
[1] "NCBI NLM NIH"

$name
[1] "Gene Expression Omnibus (GEO)"

$order
[1] "none"

$platform
[1] "GPL570"

$platform_organism
[1] "Homo sapiens"

$platform_technology_type
[1] "in situ oligonucleotide"

$pubmed_id
[1] "26109433"

$ref
[1] "Nucleic Acids Res. 2005 Jan 1;33 Database Issue:D562-6"

$reference_series
[1] "GSE69078"

$sample_count
[1] "12"

$sample_id
[1] "GSM1692587,GSM1692588,GSM1692589"                                 
[2] "GSM1692590,GSM1692591,GSM1692592"                                 
[3] "GSM1692593,GSM1692594,GSM1692595"                                 
[4] "GSM1692596,GSM1692597,GSM1692598"                                 
[5] "GSM1692587,GSM1692588,GSM1692589,GSM1692590,GSM1692591,GSM1692592"
[6] "GSM1692593,GSM1692594,GSM1692595,GSM1692596,GSM1692597,GSM1692598"

$sample_organism
[1] "Homo sapiens"

$sample_type
[1] "RNA"

$title
[1] "Multiple myeloma cell lines with acquired resistance to chemotherapeutic agent carfilzomib"

$type
[1] "Expression profiling by array" "cell line"                     "cell line"                    
[4] "cell line"                     "cell line"                     "cell type"                    
[7] "cell type"                    

$update_date
[1] "Jul 11 2016"

$value_type
[1] "transformed count"

$web_link
[1] "http://www.ncbi.nlm.nih.gov/geo"
```



<!-- rnb-output-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->

## Metadados


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuTWV0YShHRFM1ODI2KSRjaGFubmVsX2NvdW50XG5gYGAifQ== -->

```r
Meta(GDS5826)$channel_count
```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiWzFdIFwiMVwiXG4ifQ== -->

```
[1] "1"
```



<!-- rnb-output-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuTWV0YShHRFM1ODI2KSRkZXNjcmlwdGlvblxuYGBgIn0= -->

```r
Meta(GDS5826)$description
```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiWzFdIFwiQW5hbHlzaXMgb2YgcHJvdGVhc29tZSBpbmhpYml0b3IgY2FyZmlsem9taWItcmVzaXN0YW50IG11bHRpcGxlIG15ZWxvbWEgKE1NKSBjZWxsIGxpbmVzIEtNUy0xMS9DZnogYW5kIEtNUy0zNC9DZnosIGFmdGVyIDEgd2VlayBvZiBncm93dGggaW4gdGhlIGFic2VuY2Ugb2YgY2FyZmlsem9taWIuIFJlc3VsdHMgcHJvdmlkZSBpbnNpZ2h0IGludG8gdGhlIG1vbGVjdWxhciBtZWNoYW5pc21zIHVuZGVybHlpbmcgdGhlIGFjcXVpc2l0aW9uIG9mIHByb3RlYXNvbWUgaW5oaWJpdG9yIHJlc2lzdGFuY2UgaW4gTU0uXCJcblsyXSBcIktNUy0xMS9DZnpcIiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXG5bM10gXCJLTVMtMzQvQ2Z6XCIgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFxuWzRdIFwiS01TLTExXCIgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcbls1XSBcIktNUy0zNFwiICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXG5bNl0gXCJjYXJmaWx6b21pYi1yZXNpc3RhbnQgTU1cIiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFxuWzddIFwicGFyZW50YWwgTU1cIiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcbiJ9 -->

```
[1] "Analysis of proteasome inhibitor carfilzomib-resistant multiple myeloma (MM) cell lines KMS-11/Cfz and KMS-34/Cfz, after 1 week of growth in the absence of carfilzomib. Results provide insight into the molecular mechanisms underlying the acquisition of proteasome inhibitor resistance in MM."
[2] "KMS-11/Cfz"                                                                                                                                                                                                                                                                                         
[3] "KMS-34/Cfz"                                                                                                                                                                                                                                                                                         
[4] "KMS-11"                                                                                                                                                                                                                                                                                             
[5] "KMS-34"                                                                                                                                                                                                                                                                                             
[6] "carfilzomib-resistant MM"                                                                                                                                                                                                                                                                           
[7] "parental MM"                                                                                                                                                                                                                                                                                        
```



<!-- rnb-output-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuTWV0YShHRFM1ODI2KSRmZWF0dXJlX2NvdW50XG5gYGAifQ== -->

```r
Meta(GDS5826)$feature_count
```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiWzFdIFwiNTQ2NzVcIlxuIn0= -->

```
[1] "54675"
```



<!-- rnb-output-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuTWV0YShHRFM1ODI2KSRwbGF0Zm9ybVxuYGBgIn0= -->

```r
Meta(GDS5826)$platform
```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiWzFdIFwiR1BMNTcwXCJcbiJ9 -->

```
[1] "GPL570"
```



<!-- rnb-output-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuTWV0YShHRFM1ODI2KSRzYW1wbGVfaWRcbmBgYCJ9 -->

```r
Meta(GDS5826)$sample_id
```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiWzFdIFwiR1NNMTY5MjU4NyxHU00xNjkyNTg4LEdTTTE2OTI1ODlcIiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFxuWzJdIFwiR1NNMTY5MjU5MCxHU00xNjkyNTkxLEdTTTE2OTI1OTJcIiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFxuWzNdIFwiR1NNMTY5MjU5MyxHU00xNjkyNTk0LEdTTTE2OTI1OTVcIiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFxuWzRdIFwiR1NNMTY5MjU5NixHU00xNjkyNTk3LEdTTTE2OTI1OThcIiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFxuWzVdIFwiR1NNMTY5MjU4NyxHU00xNjkyNTg4LEdTTTE2OTI1ODksR1NNMTY5MjU5MCxHU00xNjkyNTkxLEdTTTE2OTI1OTJcIlxuWzZdIFwiR1NNMTY5MjU5MyxHU00xNjkyNTk0LEdTTTE2OTI1OTUsR1NNMTY5MjU5NixHU00xNjkyNTk3LEdTTTE2OTI1OThcIlxuIn0= -->

```
[1] "GSM1692587,GSM1692588,GSM1692589"                                 
[2] "GSM1692590,GSM1692591,GSM1692592"                                 
[3] "GSM1692593,GSM1692594,GSM1692595"                                 
[4] "GSM1692596,GSM1692597,GSM1692598"                                 
[5] "GSM1692587,GSM1692588,GSM1692589,GSM1692590,GSM1692591,GSM1692592"
[6] "GSM1692593,GSM1692594,GSM1692595,GSM1692596,GSM1692597,GSM1692598"
```



<!-- rnb-output-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuTWV0YShHRFM1ODI2KSRzYW1wbGVfY291bnRcbmBgYCJ9 -->

```r
Meta(GDS5826)$sample_count
```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiWzFdIFwiMTJcIlxuIn0= -->

```
[1] "12"
```



<!-- rnb-output-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuTWV0YShHRFM1ODI2KSRzYW1wbGVfb3JnYW5pc21cbmBgYCJ9 -->

```r
Meta(GDS5826)$sample_organism
```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiWzFdIFwiSG9tbyBzYXBpZW5zXCJcbiJ9 -->

```
[1] "Homo sapiens"
```



<!-- rnb-output-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuTWV0YShHRFM1ODI2KSRzYW1wbGVfdHlwZVxuYGBgIn0= -->

```r
Meta(GDS5826)$sample_type
```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiWzFdIFwiUk5BXCJcbiJ9 -->

```
[1] "RNA"
```



<!-- rnb-output-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuTWV0YShHRFM1ODI2KSR0aXRsZVxuYGBgIn0= -->

```r
Meta(GDS5826)$title
```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiWzFdIFwiTXVsdGlwbGUgbXllbG9tYSBjZWxsIGxpbmVzIHdpdGggYWNxdWlyZWQgcmVzaXN0YW5jZSB0byBjaGVtb3RoZXJhcGV1dGljIGFnZW50IGNhcmZpbHpvbWliXCJcbiJ9 -->

```
[1] "Multiple myeloma cell lines with acquired resistance to chemotherapeutic agent carfilzomib"
```



<!-- rnb-output-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuTWV0YShHRFM1ODI2KSR0eXBlXG5gYGAifQ== -->

```r
Meta(GDS5826)$type
```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiWzFdIFwiRXhwcmVzc2lvbiBwcm9maWxpbmcgYnkgYXJyYXlcIiBcImNlbGwgbGluZVwiICAgICAgICAgICAgICAgICAgICAgXCJjZWxsIGxpbmVcIiAgICAgICAgICAgICAgICAgICAgXG5bNF0gXCJjZWxsIGxpbmVcIiAgICAgICAgICAgICAgICAgICAgIFwiY2VsbCBsaW5lXCIgICAgICAgICAgICAgICAgICAgICBcImNlbGwgdHlwZVwiICAgICAgICAgICAgICAgICAgICBcbls3XSBcImNlbGwgdHlwZVwiICAgICAgICAgICAgICAgICAgICBcbiJ9 -->

```
[1] "Expression profiling by array" "cell line"                     "cell line"                    
[4] "cell line"                     "cell line"                     "cell type"                    
[7] "cell type"                    
```



<!-- rnb-output-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->





<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuXG5leHAgPSBleHBycyhkYXRhKVxuc2RzID0gcm93U2RzKGV4cClcbnN1bShpcy5uYShzZHMpKVxuYGBgIn0= -->

```r

exp = exprs(data)
sds = rowSds(exp)
sum(is.na(sds))
```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiWzFdIDYyXG4ifQ== -->

```
[1] 62
```



<!-- rnb-output-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuc2RzPW5hLmV4Y2x1ZGUoc2RzKVxubSA9IG1lZGlhbihzZHMsIG5hLnJtPVQpXG5tXG5gYGAifQ== -->

```r
sds=na.exclude(sds)
m = median(sds, na.rm=T)
m
```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiWzFdIDAuMTczNzE0MVxuIn0= -->

```
[1] 0.1737141
```



<!-- rnb-output-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuaGlzdChzZHMsIGJyZWFrcyA9IDUwLCBjb2wgPSBcIm1pc3R5cm9zZVwiKVxuYWJsaW5lKHY9bSwgY29sPVwiYmx1ZVwiLGx3ZCA9IDQsbHR5ID0gMilcbmBgYCJ9 -->

```r
hist(sds, breaks = 50, col = "mistyrose")
abline(v=m, col="blue",lwd = 4,lty = 2)
```

<!-- rnb-source-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuYWJsaW5lKHY9bSoyLCBjb2w9XCJyZWRcIixsd2QgPSA0LGx0eSA9IDIpXG5gYGAifQ== -->

```r
abline(v=m*2, col="red",lwd = 4,lty = 2)
```

<!-- rnb-source-end -->

<!-- rnb-plot-begin eyJoZWlnaHQiOjQzMi42MzI5LCJ3aWR0aCI6NzAwLCJzaXplX2JlaGF2aW9yIjowLCJjb25kaXRpb25zIjpbXX0= -->

<img src="data:image/png;base64,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tWuX9hC9dryUSIcPH9Zu8uabb0ZERLi4uPz888+meK/3SXkLtar1o9NNIrl+/XrdVhs2bDBqJ0BDwzUNgHzdu3dv2rRpQoiUlJS5c+e6uLj897//zc/PF0LovpGoO7n+n//8p1evXo0aNRo/fvzy5ctPnjz58ssvx8XFXb16VftlBFtb20mTJmkHjxw5ctmyZRcuXHjllVdee+21LVu2JCcnS4nk4uKifZCYmOju7n7lypU5c+Zolzz0AoU6cHR0lPIWDKv1o3v66adbtmx56dKlHTt2jBs3rn///jt27NiyZYtROwEaHAtehAk0HPrmadD9qH7otycqKysf+rULOzs73XccDh06pFuunachLy/vySefrLGJQqGYN2+e7nWzs7NrXIuge/rgPA3VA1+/fl13YF+rbdu22gcTJkzQjtF+e+LJJ5/UbfXgklatWtVYUp2Ut2D42xNSProVK1bUWKu7pkH77QkpOwEaFE5PAPJlY2OzcePG5cuXP/vss15eXvb29i1bthw4cOC///1v3dV5vXr1evvtt7WTM2q/raBSqX766ae5c+d27969SZMmLVu2jIyM3L9//+zZs3V79vb2Pn78+CuvvNK8efOWLVtOmTJF4i2aPD09d+/e3bNnz0aNGnXo0OF//ud/Dh48qD2bkJiYeO/ePZO8cSlvwTApH93EiROTk5OjoqI8PDxatWo1bty4zz//3NidAA2KQmOiI4oAHiPam206OTkFBQVplyxatGj69OlCiKysLH1fiADQwFEagIYoODj4zJkzNjY2X3/9db9+/VJTU2NjY3NycoKDg0+fPm3pdABkitIANERHjhyJjo6u8Q1DLy+v5OTk9u3bWyoVAJmjNAANVEFBwapVq3799dcbN274+Ph07tx5+PDhTk5Ols4FQL4oDQAAQBK+PQEAACShNAAAAEkoDQAAQBJKAwAAkITSAAAAJKE0AAAASSgNAABAEkoDAACQRGnpANZvxw6xevXvj6OixLBhxmyckyPi4n5/7Okpli0zbTYAAKSjNJjd2bNi48bfH/v5GblxYeH9jf39KQ0AAAvi9AQAAJCE0gAAACShNAAAAEkoDQAAQBJKg9l5eNx/7O5u5MZqtVD+cbGq0RsDAGBKlAazGzxYPP+8EEKEhIgxY4zcuFkzkZAg7OyEm5tYsMAM6QAAkEqh0WgsnaFBKCoSzs513bi0VCiVwtbWlIEAADASpQEAAEjC5E6mlJWVdfToUcNjnJ2d+/Xrp1Ao6icSAACmwpEGU4qPjz+YvC+wdWsDY3YdSDn9229+Rs8NCQCAhcnuSENZWVlBQYGtra1KpbJ0FqNpNJqRsX+dNGqkgTGtu/esqqqqt0gAAJiKXL49kZWVNWvWLH9/f0dHR3d3d7Va7ejoGBgY+Oabb2ZkZFg6HQAAkMeRhp9//rlXr16urq5RUVHt2rVTq9UajSY/P//cuXMbNmxITExMSUnp1KmTpWMCANCgyaI0TJkyJTQ0dPPmzU5OTjVWLV26dNiwYdOmTdu3b59FsgEAAC1ZnJ5ITU0dPXr0g41BCKFUKidOnHjixIn6TwUAAKqTRWkICgo6cOCAvrX79+9v27ZtfeYBAAAPksXpiZkzZ8bExGRkZAwePDg4OFilUikUCu01DVu3bt22bds333xj6YwAADR0sigNgwYNSkpKWrx48ejRo6svVygUYWFhSUlJffv2tVQ2AACgJYvSIISIjIyMjIy8fft2dnb2tWvXhBCenp4+Pj5qtdrS0R5VZaWYOk+15TvnsB4lie/lGb39xo1i1izh7i4++0w88YQZAgIAIIlcSoOWSqVycXHx9PR8TCd3eqhte50//LyxEGL1JpeuncqM27iwUIwYIYqKRHq6mDRJJCebJSIAABLI4kJIYdWTO2Vm3b87ZeYVI+9UmZMjiop+f/yYfw4AgMedLI40WPfkThoN96YCAFgDWZQGJncCAED+ZHF6gsmdAACQP1mUBiZ3AgBA/mRxeoLJnQAAkD9ZlAYmdwIAQP5kURqEVU/uBACAdZBLadCyysmdAACwDrK4EFJY9eROAABYB1kcabDuyZ0AALAOsigNTO4EAID8yaI0pKamfvHFFwYmd4qOjpaynyNHjkyaNEnf2qtXr27durVbt251DwoAQAMmi9Kgndxp0KBBD10rfXKnZ555JjExUd/akSNHVlVV1TEiAAANnixKg6kmd3JwcHjmmWf0rXVxcVEouHcUAAB1JIvSwOROAADInyxKg7DqyZ083Cp1j91djTw/olYLpVJUVAghhLu7SXMBAGAcuczToKVSqdq3bx8eHt67d29nZ+fy8nJLJzKBwZFFz4eWCCFCOpeNeeWecRs3ayYSEoSdnXBzEwsWmCUfAADSyKI0vP766wcPHtQ9XbRokUqlatOmjZeXV8uWLTdt2mTBbI/OyVGz56vcwvSsYzuuu6qMvxJzxgxx9664fl2EhZkhHQAAUsmiNKxaterMmTPax4mJidOnTx8wYMCWLVu2b98eERERGxu7e/duyyZ8dM5Omrpv7OAgbG1NlwUAgLqQyzUNOitWrBg7dqzum5PR0dH29vYJCQkRERGWDQYAQAMniyMN1V24cCEyMrL6koiIiLS0NEvlAQAAWrIrDW3btr106VL1JampqT4+PpbKAwAAtORSGubMmRMeHj5u3Dh3d/e5c+devHhRCKHRaL7++uvly5e/8MILlg4IAEBDJ4trGrZs2XLxD5mZmXfv3j1+/Li/v39qauqQIUO6des2f/58S2cEAKChk0VpeOmll6o/raqqqqioEEL4+vru27evd+/eNjZyOSICAECDJaMfxpWVlRkZGSUlJTY2Nvb29kIIDw+P8PDwsrKynJwcS6cDAKChk0VpqKiomDNnTuPGjf39/dVq9YwZMyor70+9/M0333h5eVkwHgAAEDI5PbF06dIFCxa88cYb3bt3P3z48JIlS27cuPHFF19YOhcAALhPFqVh1apV06dP117tOGjQoGeeeWbo0KEDBgzo37+/paMBAIDfyeL0xNWrV3v16qV7OmTIkKFDh8bHx5eUlFgwlalUVor4d1R+IT7D411LShVGb79xowgKEj17itOnzZAOAACpZFEannjiif3791dfsnjx4sLCwmnTplkqkglt2+v84eeNs7JtV29yWfV1I+M2LiwUI0aI9HRx5IiYNMk8AQEAkEQWpydee+21SZMmVVRUREVF/fnPf3ZwcPDw8Pjyyy/79+9/9+5db29vSwd8JJlZ9+81lXnFyPtO5eSIoqLfH2dkmC4UAABGk0VpiIuLu3PnzsKFCz/66KPz588HBAQIIfr167d169Zx48ZlZ2dbOuAj0WiMPyUBAID8yOL0hBDi7bffvnHjxoULF1q0aKFbGBUVdenSpf3793/yyScWzAYAAIRMjjRo2dvb+/v711ioVCp79+7du3dvi0QCAAA6cjnSAAAAZI7SAAAAJKE0AAAASSgNAABAEkoDAACQhNIAAAAkoTQAAABJKA0AAEASSgMAAJCE0gAAACShNAAAAEkoDQAAQBJKg9l5uFXqHru7Vhm3sVotlH/cVMzd3XShAAAwGqXB7AZHFj0fWiKECOlcNuaVe8Zt3KyZSEgQdnbCzU0sWGCWfAAASCOjW2NbKydHzZ6vcouKFc5OmrpsP2OGiI8XSqWwtTV1NAAAjEBpqCd1bAxaDg6mCwIAQB1xegIAAEhCaQAAAJJQGgAAgCSUBgAAIAmlAQAASEJpAAAAklAaAACAJJQGAAAgCaUBAABIYlUzQlZVVV2+fLmq6uE3hSorK6vnPAAAWBOrKg2HDh0aNWqUvrXZ2dk5OTn1mUerslJMnafa8p1zWI+SxPfyjN5+40Yxa5ZwdxeffSaeeMIMAQEAkMSqSkNoaOiFCxf0re3evbunp2d95tHattf5w88bCyFWb3Lp2snIox2FhWLECFFUJNLTxaRJIjnZLBEBAJCAaxrMLjPr/t0pM68YeafKnBxRVPT744wM04UCAMBolAaz02gUlo4AAIAJUBoAAIAklAYAACAJpQEAAEhCaQAAAJJQGgAAgCSUBgAAIAmlAQAASEJpAAAAklAaAACAJJQGAAAgCaUBAABIQmkAAACSUBoAAIAklAYAACAJpcHsPNwqdY/dXauM21itFkrlHxu7my4UAABGozSY3eDIoudDS4QQIZ3Lxrxyz7iNmzUTCQnCzk64uYkFC8ySDwAAaZS1D8GjcXLU7Pkqt6hY4eykqcv2M2aI+HihVApbW1NHAwDACJSGelLHxqDl4GC6IAAA1BGnJwAAgCSUBgAAIAmlAQAASEJpAAAAkugtDe+//35WVlZ9RgEAAHKmtzTMnz+/ZcuWoaGhn376aV5eXn1mAgAAMqS3NOTk5Hz77bfNmzefMmWKl5dX//79v/nmm6KiInMHKisru3nz5u3bt839QgAAwCh6S4ODg8OLL764bt263NzctWvXKpXK4cOHe3p6Dhs2bPfu3RUVFabNkZWVNWvWLH9/f0dHR3d3d7Va7ejoGBgY+Oabb2ZkZJj2tQAAQB3UPrmTs7Nz586dMzMzz58/n5aWtnXr1rVr13p6er777rvDhw83SYiff/65V69erq6uUVFR7dq1U6vVGo0mPz//3LlzGzZsSExMTElJ6dSpk0leCwAA1I2h0nDixImtW7d+++23p0+f9vLyGjBgwJIlS/7yl79cuXJl/vz5o0aNioiI8PT0fPQQU6ZMCQ0N3bx5s5OTU41VS5cuHTZs2LRp0/bt2/foLwQAAOpM7+kJX1/frl27rlmz5vnnnz906NDVq1dXrlwZHh6uVCpbtWq1ZMmSqqqqM2fOmCREamrq6NGjH2wMQgilUjlx4sQTJ06Y5IUsorJSxL+j8gvxGR7vWlKqMHr7jRtFUJDo2VOcPm2GdAAASKX3SMNrr702cODALl26PHStk5PT+fPnW7ZsaZIQQUFBBw4cGDRo0EPX7t+/v23btiZ5IYvYttf5w88bCyFWb3Lp2qnMuI0LC8WIEaKoSKSni0mTRHKyWSICACCB3tIwf/78rKysjz76aOzYsY6OjufOnduxY0dsbKyPj48QwtbWNiAgwFQhZs6cGRMTk5GRMXjw4ODgYJVKpVAotNc0bN26ddu2bd98842pXqv+ZWbdvztl5hUj71SZkyN031jhglAAgEXpPT1x6tSpJ554Yvr06eXl5UKIoqKi+fPnP/nkk8ePHzd5iEGDBiUlJZWUlIwePbpHjx7BwcHt2rV79tlnhw8fnp+fn5SUpO8gxGNBozH+lAQAAPKj90jDlClTOnfuvGHDhsaNGwshOnfufOXKlb/+9a//+Mc/UlJSTJ4jMjIyMjLy9u3b2dnZ165dE0J4enr6+Pio1WqTvxYAAKgDvaXh+PHjK1as8Pb21i1xdHT8+9//HhMTY740KpXKxcXF09PT1tZWpVKZ74UAAICx9J6ecHd3z8/Pr7EwIyPDw8PDHDmY3AkAAJnTe6QhNjb2rbfecnd3HzBggFKprKqq2rlz51tvvTV69GiTh2ByJwAA5E9vaXjnnXeys7NjY2NtbGzc3Nzy8vLKyspiYmISEhJMHoLJnQAAkD+9pcHW1nbVqlUzZsw4evTo5cuXvby8unbt2rFjR3OESE1N/eKLLwxM7hQdHW2O1wUAANLVcu+JwMDAwMBAc4ew7smdAACwDnpLQ15e3pw5c3755ZfKysoaq/7973+bNoR1T+4EAIB10Fsaxo0bt3Xr1r59+9bDb/nayZ0WL15c4ypLhUIRFhaWlJTUt29fc2cAAACG6S0Ne/funT9//vTp0+snB5M7AQAgcw8vDeXl5QUFBd26davnNEzuBACAbD18cic7O7s+ffp89tln9ZaDyZ0AAJA5vacnBg4c+Pbbbz/99NMRERFqtVqhuH/XpalTp5o2BJM7AQAgf3pLw4IFC1xcXG7cuLFmzZoaq0xeGkw1uVNOTs62bdv0rc3NzS0tLX3UrAAANFR6S0NWVla9hTDV5E63bt366aef9K0tLCwsKSmpe8q68nC7/51Vd9cq4zZWq4VSKSoqhBDC3d2kuQAAME4tkzsJIbKysu7cufPkk0+aL4SpJnd64oknEhMT9a09efJk06ZN6xjxEQyOLFq7xWXvQceQzmVjXrn3Sc0DNwY1ayYSEsTs2aJpU7FggbkiAgAggaHSsGfPnhEjRly/fl0IodFowsPDX3zxxbi4OJOHsO7JnZwcNXu+yi0qVjg7aeqy/YwZIj5eKJXC1tbU0QAAMILe0rB27dqRI0eOGjXqT3/602uvvSaE6NGjxxtvvOHg4DB27FjThmgIkzvVsTFoOTiYLggAAHVk6ELIuLi4JUuW3Lp1S7tk3rx5JSUly5YtM3lpEEzuBACA7OktDZmZmc8991yNhWFhYR9//LH50qhUKpVK1b59+6qqqoyMjPLycvO9FgAAMMrDJ3cSQrRr1+7YsWM1Fp44cSIgIMDkIV5//fWDBw/qni5atEilUrVp08bLy6tly5abNm0y+SsCAABj6T3SEBcXN3bsWKVS2adPHyHEjRs3kpKSEhIS3n33XZOHWLVqVZcuXUJDQ4UQiYmJ06dPHzZs2IABA5RK5Y4dO2JjY3fu3BkREWHy1wUAANLpLQ0jRoy4e/fu3LlzZ8+eLYTw8PBwcHCYMmVKfHy8WQOtWLFi7Nixum9ORkdH29vbJyQkUBoAALAsQ1+5jIuLGzVq1KlTpzIzM93d3Tt06ODm5mbuQBcuXEhISKi+JCIiYvXq1eZ+XQAAYFgtkzu5uLiEhISEhITUTxohRNu2bS9dulR9SWpqqo+PT70FAAAAD6W3NBiYuTkpKcnkOebMmbNx48aAgAB3d/e5c+dGRkb6+/trNJp169YtX7781VdfNfkrAgAAo+gtDTXOROTn5x85cuTmzZsTJ040eYgtW7Zc/ENmZubdu3ePHz/u7++fmpo6ZMiQbt26zZ8/3+QvCgAAjKK3NHz55Zc1lhQWFg4cOPDB72E+updeeqn606qqqoqKCiGEr6/vvn37evfubWOj96uhAACgfhjxw9jFxeWtt946fvz4jRs3zBdICGFjY2Nvby+E8PDwCA8Pf9wbQ2WliH9H5RfiMzzetaRUYfT2GzeKoCDRs6c4fdoM6QAAkKr2u1xWl5mZ6eTkVA/fobAm2/Y6f/h5YyHE6k0uXTuVGbdxYaEYMUIUFYn0dDFpkkhONktEAAAk0FsavvrqqxpLzp8/n5iY2K1bN4XC+F+XG7DMrPt3p8y8YuSdKnNyRFHR748zMkwXCgAAo+ktDWPGjKmxxMbGpmPHjitXrjR5iO3btx86dMjwmIULF5r8deuHRkPHAgBYA72lobi4uN5CODs7//DDD//5z38cHR1btGjx0DGPb2kAAMA6GHdNg5mEh4eHhYWFhYVpNJpaDzkAAACL0Fsa3N3dDW/p6Oh4+fJlU13fYGtrGxsbu379epPsDQAAmJze0rBo0aJRo0a1bNly8ODBzZs3v379+qZNm27evDlv3jwHBwftmJKSEicnJ1NFefHFF/39/U21NwAAYFp6S8Pu3bv79Omza9cupfL3MQkJCVFRUWfOnPn444/NEaVFixb6LmgAAAAWp3fepJSUlPHjx+sagxBCqVROmDBh+/bt9RIMAADIi97S4ODgcPny5RoLL1++XFVVZeZIAABAjvSWhhdffPGdd96pfkPLnTt3zp49+/nnn6+XYAAAQF70XtPw/vvvZ2Rk9O/fX61WN2/ePDs7Oy8vLyQk5IMPPqjPfAAAQCb0lgZHR8ekpKRjx44dOXLk8uXLHh4eTz31VERERH2GAwAA8lHL5E7dunVr3rz5nTt3nnzyyfoJBAAA5MnQXaf37Nnj7e3t5+fXoUMHIUR4ePiyZcvqKxgAAJAXvaVh7dq1UVFR/fv3X7NmjXZJjx493njjjU8//bS+sgEAABnRe3piwYIFcXFxS5YsuXXrlnbJvHnzSkpKli1bNnbs2PqKBwAA5ELvkYbMzMznnnuuxsKwsLCMjAwzR7I2Hm6VusfurkbOcqFWC938WrXdDQQAALPSWxratWt37NixGgtPnDgREBBg5kjWZnBk0fOhJUKIkM5lY165Z9zGzZqJhARhZyfc3MSCBWbJBwCANHpPT8TFxY0dO1apVPbp00cIcePGjaSkpISEhHfffbce41kDJ0fNnq9yi4oVzk6aumw/Y4aIjxdKpbC1NXU0AACMoLc0jBgx4u7du3Pnzp09e7YQwsPDw8HBYcqUKfHx8fUYz3rUsTFo/XFbUQAALOjhpaG8vDw7O/v1118fNWrUqVOnMjMz3d3dO3To4ObmVs/5AACATDy8NFRVVT311FOrVq0aOHBgSEhISEhIPccCAABy8/ALIR0cHEaNGrV69WqN5hEOqgMAACui95qGkJCQgwcPau834enpaWNzv1688cYb9ZINAADIiN7SoGsG//rXv/StAgAADUfN0rB///4OHTp4eHhcu3bNIoEAAIA81bymITw8/NChQ7qn//u//3vmzJn6jQQAAOTI0F0uhRBvv/32qVOn6icKAACQs1pKAwAAgBalAQAASEJpMLvKShH/jsovxGd4vMJ7k4gAABl8SURBVGtJqcLo7TduFEFBomdPcfq0GdIBACCV3q9cwlS27XX+8PPGQojVm1y6diozbuPCQjFihCgqEunpYtIkkZxslogAAEjwkCMNI0eOdP9Djae6hZAuM+v+3Skzrxh5p8qcHFFU9PvjjAzThQIAwGg1jzQ81jexTE1NNXDn7vT09Fu3btVnHi2NRvF/n2ru3Llz+/ZtA5s4Ojo6OTmZORcAAMapWRqWLl1qkRwmERAQEBsbW1VV9dC1J0+ebNq0aT1HelB2Tk7on/5UfVruGiqrKoOfaH/06NH6TAUAQK2s6pqGpk2bDho0SN/aJUuWKJWWf78ajfhl7+5Wvi30DTj522+vTf1HfUYCAEAKvj0BAAAkoTQAAABJKA0AAEASSgMAAJCE0gAAACShNAAAAEkoDQAAQBJKAwAAkITSAAAAJKE0AAAASSgNAABAEkoDAACQhNJgdh5ulbrH7q4PvwOnXmq10N1ky93ddKEAADAapcHsBkcWPR9aIoQI6Vw25pV7xm3crJlISBB2dsLNTSxYYJZ8AABIY/lbRVs9J0fNnq9yi4oVzk6aumw/Y4aIjxdKpbC1NXU0AACMQGmoJ3VsDFoODqYLAgBAHXF6AgAASEJpAAAAklAaAACAJJQGAAAgCaUBAABIQmkAAACSUBoAAIAklAYAACAJpQEAAEhCaQAAAJJQGsyuslLEv6PyC/EZHu9aUqowevuNG0VQkOjZU5w+bYZ0AABIxb0nzG7bXucPP28shFi9yaVrpzLjNi4sFCNGiKIikZ4uJk0SyclmiQgAgAQcaTC7zKz7d6fMvGLknSpzckRR0e+PMzJMFwoAAKNRGsxOozH+lAQAAPIju9MTZWVlBQUFtra2KpXK0lkAAMB9cjnSkJWVNWvWLH9/f0dHR3d3d7Va7ejoGBgY+Oabb2ZwWB4AABmQxZGGn3/+uVevXq6urlFRUe3atVOr1RqNJj8//9y5cxs2bEhMTExJSenUqZOlYwIA0KDJojRMmTIlNDR08+bNTk5ONVYtXbp02LBh06ZN27dvn0WyAQAALVmcnkhNTR09evSDjUEIoVQqJ06ceOLEifpPBQAAqpNFaQgKCjpw4IC+tfv372/btm195gEAAA+SxemJmTNnxsTEZGRkDB48ODg4WKVSKRQK7TUNW7du3bZt2zfffGPpjAAANHSyKA2DBg1KSkpavHjx6NGjqy9XKBRhYWFJSUl9+/a1VDYAAKAli9IghIiMjIyMjLx9+3Z2dva1a9eEEJ6enj4+Pmq12tLRAACAEPIpDVoqlcrFxcXT05PJnQAAkBtZXAgpmNwJAADZk8WRBiZ3AgBA/mRRGpjcCQAA+ZNFaUhNTf3iiy8MTO4UHR1d/6lqSEtLGzBggOExt27dmjd1cv3kAQCgnsmiNGgndxo0aNBD18pkcqfc3FwvtWrNRx8YGBM1fNSDCz3cKnWP3V2rjHtVtVoolaKiQggh3N2N2xYAAJOSRWl4XCZ3cnJ09PfzMzDAwd7uwYWDI4vWbnHZe9AxpHPZmFfuvf2+MS/ZrJlISBCzZ4umTcWCBUbmBQDAlGRRGqx7cicnR82er3KLihXOTpq6bD9jhoiPF0qlsLU1dTQAAIwgi9IgTDS5U0pKSu/evfWtVSgUOTk5JshaJ3VsDFoODqYLAgBAHcmlNGg94uROYWFhGo3en83du3f39PR8tIAAADRcTO4EAAAkkcWRBiZ3AgBA/mRRGpjcCQAA+ZPF6YnU1NTRo0cbmNzpxIkT9Z8KAABUJ4vSoJ3cSd9amUzuBABAAyeL0xOPy+ROAAA0ZLIoDdY9uZOxcm7cPHv2bEBAgOFhb731Vo2PCwAAs5JFaRAmmtzJOty5e7e1b4udq/+fgTEr/t/qzMzMegoEAIAQQj6lQUulUqlUqvbt21dVVWVkZJSXl1s6kWXY29kbvsmFqmnT0npLAwCAEEImF0K+/vrrBw8e1D1dtGiRSqVq06aNl5dXy5YtN23aZMFsj66yUsS/o/IL8Rke71pSqjB6+x07xZ/+IgYMFOfSzZAOAACpZFEaVq1adebMGe3jxMTE6dOnDxgwYMuWLdu3b4+IiIiNjd29e7dlEz6KbXudP/y8cVa27epNLqu+bmTcxkVFYvJUkZEhTvwk5rxjnoAAAEgir9MTQogVK1aMHTs2MTFR+zQ6Otre3j4hISEiIsKyweosM+v+3Skzrxh5p8qbN0Vx8e+PL2eZLhQAAEaTxZGG6i5cuBAZGVl9SURERFpamqXyPDqNxvhTEgAAyI/sSkPbtm0vXbpUfUlqaqqPj4+l8gAAAC25lIY5c+aEh4ePGzfO3d197ty5Fy9eFEJoNJqvv/56+fLlL7zwgqUDAgDQ0MnimoYtW7Zc/ENmZubdu3ePHz/u7++fmpo6ZMiQbt26zZ8/39IZAQBo6GRRGl566aXqT6uqqioqKoQQvr6++/bt6927t42NXI6IAADQYMmiNNRgY2Njb28vhPDw8AgPD7d0HAAAIIR8rmkAAAAyR2kAAACSUBoAAIAklAYAACAJpQEAAEgix29PoFa/nj37w/ET3333nYExNjY2mzZt8jN4i20AAKSjNDyWbublRf4ldPxrQw2MGRY/OScnh9IAADAVSsPjytvT45mOHQwMcHFyrrcwAICGgGsaAACAJJQGs/Nwq9Q9dnetMm7jZs2E8o+jQa6upgsFAIDRKA1mNziy6PnQEiFESOeyMa/cM27jJk3EP6YJpVKo1eLNGWbJBwCANFzTYHZOjpo9X+UWFSucnTR12X7iBDFmtLC1Fba2po4GAIARKA31pI6NQcve3nRBAACoI05PAAAASTjSYLUuXr780ksvOTg4GBjj7++/b9++eosEAHisURqsVklp6aqF73V8IljfgNybtwaOm1CfkQAAjzVKgzXz8fby1z8jpJOjY32GAQA87rimAQAASEJpAAAAklAaAACAJFzT0HCVlJbevn37ueeeMzxs5MiRr776av1EAgDIGaWh4bpz966zo+OM0SMNjNm5/8APP/xAaQAACEpDA+dgbx/+p14GBqRnZBxKO3Xx4kUDY+zt7X18fBQKhanTAQDkhdJgdpWVYuo81ZbvnMN6lCS+l2f09jt2ivcWCle1eP89ERRohoCGHE39eeeePccOHzYw5vqNG9uTkvr06VNvqQAAFkFpMLtte50//LyxEGL1JpeuncqM27ioSEyeKoqLRUaGmPOOWP+1WSLqV1FZ8eqAASsX/K+BMZHDRpSWltZbJACApVAazC4z6/7dKTOvGHmnyps3RXHx748vZ5kulCmdu5jxxhtvzJkzx8AYlUrFfNUA8LiTXWkoKysrKCiwtbVVqVSWzmIaGo2Vn+y/c/fu22/E9ezSRd+AisqK7v1fqvVrGpGRkfHx8aZOBwAwGbmUhqysrE8++WTdunWZmZkajUYI4eDg4OvrO3jw4LFjx7Zu3drSAWFIm1atnunYQd/akpJSGxsbw1/TWL99+9tvv71s2TIDYyoqKv7+978/9dRTBsY0bdq0a9eutQYGANSBLErDzz//3KtXL1dX16ioqHbt2qnVao1Gk5+ff+7cuQ0bNiQmJqakpHTq1MnSMVF3CoXC8Nc09v7ww5+6dl0x/38MjOn915f/tWrVXk8PA2MO/Hi4Y6dOtrZ6TwNVVVX5+voOHTrUwE7Ky8tbt27t6elpYIyzs7OXl5eBAUKIwsLCsrJarmKxmiNqABoCWZSGKVOmhIaGbt682cnJqcaqpUuXDhs2bNq0aVLOiN+5cyc5Obmqquqha2/fvl1RUfEoOXNu3Ny4Y6eBAfl3Cn45dbrGmP+e7iRED+3jcxczNFVVuw4ccHd11beTo6k/37lboN2Jy82b/f5YXlhUtOuPPefevHX6XLrhMBUVFft/PJyZdUXfgMwrV0pKSw3v5PLVbCdHR8NjSstKf/zP8cKiIn0DysrLNRqN4Z2cvXDx9p07P51MMzCmpLQsKrxP6LPPGhiz/8fDoc887aj/huCpv/66ffv27du3G9iJRH5+fkql3n9BFRUVV69eraysNLyT5s2be3t7Gxhw5coVd3d3A7c412g02dnZTz75pIGd3Lt3r6KiwvARu+zsbF9f3yZNmhgYk56eHhISYmBAcXHxnTt3goP13ltVCJGVldWkSZOmTZsafqEuXbrY2Oidsra8vDwzMzMw0ND3iW7evKlQKFz1/1vThmnfvr3hlnny5MmWLVsa2ElBQYGDg4PhEnnmzBm1Wm1nZ2dgTGZmpuEDaYWFhfn5+T4+PgbGXLt2rVmzZg/+X1pdRkaG4b8MZWVl165dM/yuc3NznZycGjdubGDMhQsX/P39DXwfu7KyUvsXz8BO7ty5o9FomjVrZmDM1atXvby8DPw5ajSay5cvG35HhYWFxcXFbm5uBsZcu3bNzc3N8J/jpUuXDL+QEKJTp05BQUGGx8iTLEpDamrqF1988dC/5UqlcuLEidHR0VL2c+HChW+++UbfWgcHB8N/vw0LCAgI7thx4/4UA2O8fX0zc2/UGHPufCNdaTh76XLbdu2ST6Qa+GFTWlrauJlKuxPPwkJdabhXXKzbs8rDI6fgruEwQW3bHvn1dOp5vVMsVFZWenp7G96JfaPGBYWFhscEtAn85cLF8zm5hsYEBBjeSanCpsLGxvCYlq1bn72SnWv4XQcFXcm/Y+CHjX3jJrX+c7169aq9vb27u7uBMadOnQoKCjLwf0dZWVnjxo2feOIJAzvJzc2tqKho3ry5gTHl5eV+fn4GfgZUVVUVFhbWesSioKDA8IDi4uL8/HzDLefq1auGJ+24e/duTk5Okf4GKYS4fPlyo0aN1Gq1gTFpaWm3b9828DOgpKTk0qVLly5dMrCTa9eu2djYGD5i9Ntvv6Wlpdnb2+sbUFFRcfbs2fbt2xvYyc2bN0tLSw3/LE9PT/fx8XF2dtY3QKPRpKWlnT9/3sBO8vPz8/PzW7VqZWBMRkaGq6ur4fJ38uTJjh07GhhQWFh47dq1Nm3aGBhz5coVJycnw53s119/DQ4ONvDnWFpampGR0a5dOwM7ycnJ0Wg0tXay1q1bG+jWlZWVZ86cMfzneOvWreLi4hYtWhgYc/78eW9vbxcXFwNjav14hRAFBQWPaWlQaC8gsKyuXbuGhISsWLHioWvnzp373XffHT16tJ5TmcrixWLatN8fT50qFi0yZuOLF0VAwO+P/f3FhQumzQYAgHSyONIwc+bMmJiYjIyMwYMHBwcHq1QqhUKhvaZh69at27ZtM3D8AAAA1A9ZlIZBgwYlJSUtXrx49OjR1ZcrFIqwsLCkpKS+fftaKhsAANCSRWkQQkRGRkZGRt6+fTs7O/vatWtCCE9PTx8fH8MnPgEAQL2RS2nQUqlUKpXK8LUqAADAIvReZA4AAFAdpQEAAEhCaQAAAJJQGgAAgCTyuhBSzq5evfrhhx8aniLtoY4e7SjE7zMfnzx58tNPjZilqsnNmy//8bigoGD9p58a++qPi8LCwtzcXO5MZiZ37twpKCgwPFkv6uzWrVulpaWGZ/ZEneXm5nbr1q3W2+SifshiRsjHwqZNmyZPntyvX7/ah/5flZV2FRW/TzZsa1uuVJZL31ah0Tj9cccjjUJRrH+y28ddRkbG+fPn+X/BTM6ePXv9+vXQ0FBLB7FOaWlp9+7d6969u6WDWKfU1NTOnTt/9tlnlg4CITjSIJ1KpWrbtm1iYqKlg1inzZs3r1u3jo/XTL744ovDhw/z8ZrJ0qVLs7KylixZYukg1umf//ynpSPgPq5pAAAAklAaAACAJJQGAAAgCaUBAABIQmkAAACSUBoAAIAklAYAACAJ8zRIZWdnp1TycZkLH69Z8fGaFR+vWdnZ2Vk6Au5jRkipNBrN7du31Wq1pYNYp4qKisLCwqZNm1o6iHUqKysrLS1t3LixpYNYp5KSksrKShcXF0sHsU7FxcVCCCcnJ0sHgRCUBgAAIBHXNAAAAEkoDQAAQBJKAwAAkITSAAAAJKE0AAAASSgNAABAEkoDAACQhNIAAAAkoTQAAABJKA0AAEASSgMAAJCE0iDJ5s2bQ0JCmjVr1rt3719++cXScaxTcnLy9u3bLZ3CCi1fvvzZZ59t3Lhxu3btFi1aVFFRYelE1uPevXuTJ0/29/dv1KhRly5dNm3aZOlEVmvXrl3ffvutpVOA0iDBzp07Y2JiunTpsmrVKgcHh169emVlZVk6lLWpqqqaNWvWoUOHLB3E2iQkJMTFxfXs2XP9+vUvvfTSm2++OW/ePEuHsh4TJkz4/PPPJ02atHHjxo4dO8bExCQnJ1s6lBU6e/ZsTEzMli1bLB0EQmhQm7CwsIiICO3joqIiX1/ft956y7KRrElWVtaKFSv+/Oc/CyGmTZtm6ThWpbS0tEmTJpMmTdItmTp1qpOTU0VFhQVTWY38/HyFQrF69Wrt06qqqrZt244YMcKyqaxPWVnZM888I4QYOnSopbNAw5GGWty+fTslJSUmJkb71MnJKTIyct26dZZNZU3S0tLWrVtXVVXl6Oho6SzW5sqVKwUFBVFRUbol3bt3Ly4uvnz5sgVTWY3c3NzQ0NAePXponyoUCk9Pz6KiIsumsj6zZ8+2tbXV9gZYHKWhFtnZ2UKI4OBg3ZLg4ODMzMyysjLLhbIqffv2PXTo0KFDh1q0aGHpLNbGx8fn/Pnz2qM4WocPH3ZycvL29rZgKqsRGBiYkpISEBAghNBoNN99993x48ejo6MtncuqHDx4cMWKFWvWrLGzs7N0FghBaajV9evXhRAqlUq3RK1WazSagoICy4UCJHFwcAgICHBwcNA+/eqrr5YtW/b3v/+dgzqm9dFHHzk7O/fr12/cuHFDhw61dBzrkZ+fP2zYsIULFwYFBVk6C36ntHQAudNoNEIIhUJRY4mtra3FMgFGunHjxtSpU9esWTNixIj58+dbOo61iY6O9vPzO3LkyAcffNCiRYupU6daOpGVGD9+fPv27cePH2/pILiP0lALT09PIUR+fr5uSX5+voODQ/VjD4Cc7dq1a+TIkS4uLlu3bn3xxRctHccKtW7dunXr1gMGDKioqFi5ciWlwSTWr1+fnJz866+/WjoI/g9OT9TCx8dHoVCcO3dOtyQ9PZ2z73hc7Nq1q3///oMHDz59+jSNwbQ2bdoUGRmpPfSo1b59+wsXLnAtpEkcO3bs1q1b3t7eCoVCoVAcPXp07dq1CoVi27Ztlo7WoFEaaqFWq8PCwrZu3ap9WlFRsXPnzsGDB1s2FSBFRUXF66+//uqrr65YsYLrGEzOxcVl165dJ06c0C3597//3aJFC2dnZwumshoTJkxIriY4OLhPnz7Jycm6r6vAIjg9Ubtp06ZFR0fPmzevT58+H3/88e3bt8eNG2fpUEDtfvjhh+zsbG9v73/961/Vl8fGxtIhHl14eHiPHj1efvnlOXPmeHl5JScnf/nllx9//LGlc1mJoKCg6tc/Nm3a1Nvbu0+fPhaMBEFpkKJv377r169fuHDhkiVLunTpcuDAgdatW1s6FFC79PR0IcT7779fY3lERASl4dHZ2dnt3r17ypQpCxYsuHr1art27dasWfPqq69aOhdgRorqJ+QAAAD04ZoGAAAgCaUBAABIQmkAAACSUBoAAIAklAYAACAJpQEAAEhCaQAAAJJQGgAAgCSUBgAAIAmlAQAASEJpAAAAklAaAACAJJQGAAAgCaUBAABIQmkAAACSUBoAAIAklAYAACAJpQEAAEhCaQAAAJJQGgAAgCSUBgAAIAmlAQAASEJpAAAAklAaAACAJJQGAAAgCaUBgOklJia6u7tbOgUAE6M0AAAASSgNAABAEkoDACOcPHkyIiJCpVK5ubkNHDgwKytLu7ywsHD8+PG+vr6+vr4TJ04sKSmpdRMAjx2FRqOxdAYAj4fi4uLWrVsHBATExcUVFBTMnTu3ffv2e/fuFUL07t37p59++uc//+nr67t8+fKffvrJ0dHxxo0bBjYB8NhRWjoAgMfG6dOnc3Jyvvrqqz59+gghPD09f/zxRyFESkpKSkrK9u3bo6OjhRDR0dEBAQGlpaUGNgHwOOJIAwCpbt261bp161atWk2ZMuWFF17w9vbWLn/vvffee++9vLw83cjJkyevXbv2xo0b+jYB8DjimgYAUrm6un7//fetWrWaMGFC8+bNn3766W+//VYIcf36dR8fn+oj/fz8DG8C4HFEaQBghKeffnr79u15eXl79+718PCIiYk5e/asj49PdnZ29WG3bt0yvEm9BwdgApQGAFJt3Lixbdu2hYWFTk5Ozz333MqVKysrKzMyMkJCQvLy8nbu3KkdVllZuXnzZsObWO5NAKg7LoQEIFXnzp0zMzNjY2MnTJhw9+7dL7/80s3NLSQkRK1Wh4WFDRkyJCEhwc/P75NPPiksLDS8iWXfCIC64UgDAKkCAwM3b958/fr12NjYuLg4pVKZnJysVquFENu3b4+NjX3//fe1szWsXLmy1k0APHb49gQAAJCEIw0AAEASSgMAAJCE0gAAACShNAAAAEkoDQAAQBJKAwAAkITSAAAAJKE0AAAASSgNAABAEkoDAACQhNIAAAAkoTQAAABJKA0AAEASSgMAAJCE0gAAACShNAAAAEkoDQAAQBJKAwAAkITSAAAAJKE0AAAASSgNAABAEkoDAACQhNIAAAAkoTQAAABJ/j/GwX8ZFkmhSAAAAABJRU5ErkJggg==" />

<!-- rnb-plot-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuZGF0YXIgPSBkYXRhW3NkcyA+PTMqbWVkaWFuKHNkcyksIF1cbmRhdGFyXG5gYGAifQ== -->

```r
datar = data[sds >=3*median(sds), ]
datar
```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiRXhwcmVzc2lvblNldCAoc3RvcmFnZU1vZGU6IGxvY2tlZEVudmlyb25tZW50KVxuYXNzYXlEYXRhOiA0OTAxIGZlYXR1cmVzLCAxMiBzYW1wbGVzIFxuICBlbGVtZW50IG5hbWVzOiBleHBycyBcbnByb3RvY29sRGF0YTogbm9uZVxucGhlbm9EYXRhXG4gIHNhbXBsZU5hbWVzOiBHU00xNjkyNTg3IEdTTTE2OTI1ODggLi4uIEdTTTE2OTI1OTggKDEyIHRvdGFsKVxuICB2YXJMYWJlbHM6IHNhbXBsZSBjZWxsLmxpbmUgY2VsbC50eXBlIGRlc2NyaXB0aW9uXG4gIHZhck1ldGFkYXRhOiBsYWJlbERlc2NyaXB0aW9uXG5mZWF0dXJlRGF0YVxuICBmZWF0dXJlTmFtZXM6IDE0ODdfYXQgMTU1MjI2NF9hX2F0IC4uLiBBRkZYLXIyLUJzLXRoci1NX3NfYXQgKDQ5MDEgdG90YWwpXG4gIGZ2YXJMYWJlbHM6IElEIEdlbmUgdGl0bGUgLi4uIEdPOkNvbXBvbmVudCBJRCAoMjEgdG90YWwpXG4gIGZ2YXJNZXRhZGF0YTogQ29sdW1uIGxhYmVsRGVzY3JpcHRpb25cbmV4cGVyaW1lbnREYXRhOiB1c2UgJ2V4cGVyaW1lbnREYXRhKG9iamVjdCknXG4gIHB1Yk1lZElkczogMjYxMDk0MzMgXG5Bbm5vdGF0aW9uOiAgXG4ifQ== -->

```
ExpressionSet (storageMode: lockedEnvironment)
assayData: 4901 features, 12 samples 
  element names: exprs 
protocolData: none
phenoData
  sampleNames: GSM1692587 GSM1692588 ... GSM1692598 (12 total)
  varLabels: sample cell.line cell.type description
  varMetadata: labelDescription
featureData
  featureNames: 1487_at 1552264_a_at ... AFFX-r2-Bs-thr-M_s_at (4901 total)
  fvarLabels: ID Gene title ... GO:Component ID (21 total)
  fvarMetadata: Column labelDescription
experimentData: use 'experimentData(object)'
  pubMedIds: 26109433 
Annotation:  
```



<!-- rnb-output-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->






<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuXG5tYXhpbW9zID0gYXBwbHkoZXhwLDEsbWF4KVxubWluaW1vcyA9IGFwcGx5KGV4cCwxLG1pbilcbnZsID0gbWF4aW1vcy9taW5pbW9zID4gMlxuZGF0YWY9ZGF0YXJbdmwsXVxuYGBgIn0= -->

```r

maximos = apply(exp,1,max)
minimos = apply(exp,1,min)
vl = maximos/minimos > 2
dataf=datar[vl,]
```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiRXJyb3IgaW4gb3JpZ1tbbm1dXVtpLCAsIC4uLiwgZHJvcCA9IGRyb3BdIDogXG4gIChzdWJzY3JpcHQpIGxvZ2ljYWwgc3Vic2NyaXB0IHRvbyBsb25nXG4ifQ== -->

```
Error in orig[[nm]][i, , ..., drop = drop] : 
  (subscript) logical subscript too long
```



<!-- rnb-output-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->





<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxubWV0YSA9IHBoZW5vRGF0YShkYXRhKVxubWV0YVxuYGBgIn0= -->

```r
meta = phenoData(data)
meta
```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiQW4gb2JqZWN0IG9mIGNsYXNzICdBbm5vdGF0ZWREYXRhRnJhbWUnXG4gIHNhbXBsZU5hbWVzOiBHU00xNjkyNTg3IEdTTTE2OTI1ODggLi4uIEdTTTE2OTI1OTggKDEyIHRvdGFsKVxuICB2YXJMYWJlbHM6IHNhbXBsZSBjZWxsLmxpbmUgY2VsbC50eXBlIGRlc2NyaXB0aW9uXG4gIHZhck1ldGFkYXRhOiBsYWJlbERlc2NyaXB0aW9uXG4ifQ== -->

```
An object of class 'AnnotatedDataFrame'
  sampleNames: GSM1692587 GSM1692588 ... GSM1692598 (12 total)
  varLabels: sample cell.line cell.type description
  varMetadata: labelDescription
```



<!-- rnb-output-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxudmFyTGFiZWxzKG1ldGEpXG5gYGAifQ== -->

```r
varLabels(meta)
```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiWzFdIFwic2FtcGxlXCIgICAgICBcImNlbGwubGluZVwiICAgXCJjZWxsLnR5cGVcIiAgIFwiZGVzY3JpcHRpb25cIlxuIn0= -->

```
[1] "sample"      "cell.line"   "cell.type"   "description"
```



<!-- rnb-output-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



Analise diferencial

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuXG5kYXRhciRjZWxsLnR5cGU9IGZhY3RvcihkYXRhciRjZWxsLnR5cGUpXG50YWJsZShkYXRhciRjZWxsLnR5cGUpXG5gYGAifQ== -->

```r

datar$cell.type= factor(datar$cell.type)
table(datar$cell.type)
```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiXG5jYXJmaWx6b21pYi1yZXNpc3RhbnQgTU0gICAgICAgICAgICAgIHBhcmVudGFsIE1NIFxuICAgICAgICAgICAgICAgICAgICAgICA2ICAgICAgICAgICAgICAgICAgICAgICAgNiBcbiJ9 -->

```

carfilzomib-resistant MM              parental MM 
                       6                        6 
```



<!-- rnb-output-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxubGlicmFyeShnZW5lZmlsdGVyKVxudHQgPSByb3d0dGVzdHMoZGF0YXIsIFwiY2VsbC50eXBlXCIpXG5uYW1lcyh0dClcbmBgYCJ9 -->

```r
library(genefilter)
tt = rowttests(datar, "cell.type")
names(tt)
```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiWzFdIFwic3RhdGlzdGljXCIgXCJkbVwiICAgICAgICBcInAudmFsdWVcIiAgXG4ifQ== -->

```
[1] "statistic" "dm"        "p.value"  
```



<!-- rnb-output-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxueD10dCRwLnZhbHVlXG5yYW5rID0gb3JkZXIodHQkcC52YWx1ZSlcbnAyMCA9IHJhbmtbMToyMF1cbnR0JHAudmFsdWVbcDIwXVxuYGBgIn0= -->

```r
x=tt$p.value
rank = order(tt$p.value)
p20 = rank[1:20]
tt$p.value[p20]
```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiIFsxXSAwLjAwMDEzMTQyMDkgMC4wMDAzNDQ5MzQ3IDAuMDAxNTkxOTQ4NSAwLjAwMjU0Nzk4OTkgMC4wMDMwMzk5NDE5IDAuMDA2NzM1MjcyNCAwLjAwODExNzcwMjAgMC4wMDg0NjM0NzYyIDAuMDEwMTMxMTMyNFxuWzEwXSAwLjAxMDU1NTgzNTAgMC4wMTE0Mjc4MjAzIDAuMDEyMjI5ODYzMyAwLjAxMjg0ODA5MjMgMC4wMTMxODIyMTUwIDAuMDE1NDg3MTcwNSAwLjAxNjI3NzAxODAgMC4wMTY2MDY0MjQ5IDAuMDE5NjQ1NDA3MFxuWzE5XSAwLjAyMjE1MjIyMDkgMC4wMjMwODM5MTAwXG4ifQ== -->

```
 [1] 0.0001314209 0.0003449347 0.0015919485 0.0025479899 0.0030399419 0.0067352724 0.0081177020 0.0084634762 0.0101311324
[10] 0.0105558350 0.0114278203 0.0122298633 0.0128480923 0.0131822150 0.0154871705 0.0162770180 0.0166064249 0.0196454070
[19] 0.0221522209 0.0230839100
```



<!-- rnb-output-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->

Lista de 20 genes


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuZyA9IGZlYXR1cmVOYW1lcyhkYXRhcltwMjBdKVxuZ1xuYGBgIn0= -->

```r
g = featureNames(datar[p20])
g
```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiIFsxXSBcIjIwODcxMl9hdFwiICAgIFwiMjE5MTU5X3NfYXRcIiAgXCIyMDU4MjJfc19hdFwiICBcIjIxMzQ3OF9hdFwiICAgIFwiMjEwOTk3X2F0XCIgICAgXCIyMjg0OTlfYXRcIiAgICBcIjIyNTA2MF9hdFwiICAgIFwiMjIyODM4X2F0XCIgICBcbiBbOV0gXCIyMDQ1ODhfc19hdFwiICBcIjIwNTgzOV9zX2F0XCIgIFwiMjM2MjgwX2F0XCIgICAgXCIyMjY1MTdfYXRcIiAgICBcIjIwODk4M19zX2F0XCIgIFwiMjA1NTQ5X2F0XCIgICAgXCIyMzQzMDZfc19hdFwiICBcIjIwODcxMV9zX2F0XCIgXG5bMTddIFwiMjA4NjgzX2F0XCIgICAgXCIxNTU4NDM4X2FfYXRcIiBcIjIyNTQxNV9hdFwiICAgIFwiMjE0MDIzX3hfYXRcIiBcbiJ9 -->

```
 [1] "208712_at"    "219159_s_at"  "205822_s_at"  "213478_at"    "210997_at"    "228499_at"    "225060_at"    "222838_at"   
 [9] "204588_s_at"  "205839_s_at"  "236280_at"    "226517_at"    "208983_s_at"  "205549_at"    "234306_s_at"  "208711_s_at" 
[17] "208683_at"    "1558438_a_at" "225415_at"    "214023_x_at" 
```



<!-- rnb-output-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuI2ZhbHRhIHBvZGVyIGNvcnJlc3BvbmRlciBvIGlkZW50aWZpY2Fkb3IgYW8gbm9tZSBkbyBnZW5lIFxuYGBgIn0= -->

```r
#falta poder corresponder o identificador ao nome do gene 
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


Heatmap



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuZGF0YTIwPWRhdGFyW3AyMCxdXG5vcmRlcl9jb2xzID0gb3JkZXIoZGF0YTIwJGNlbGwudHlwZSlcbm9yZGVyX2NvbHNcbmBgYCJ9 -->

```r
data20=datar[p20,]
order_cols = order(data20$cell.type)
order_cols
```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiIFsxXSAgMSAgMiAgMyAgNCAgNSAgNiAgNyAgOCAgOSAxMCAxMSAxMlxuIn0= -->

```
 [1]  1  2  3  4  5  6  7  8  9 10 11 12
```



<!-- rnb-output-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuZGF0YTIwID1kYXRhMjBbLG9yZGVyX2NvbHNdXG5oZWF0bWFwKGV4cHJzKGRhdGEyMCksIENvbHYgPSBOQSwgbGFiQ29sID0gZGF0YTIwJGNlbGwudHlwZSlcblxuYGBgIn0= -->

```r
data20 =data20[,order_cols]
heatmap(exprs(data20), Colv = NA, labCol = data20$cell.type)

```

<!-- rnb-source-end -->

<!-- rnb-plot-begin eyJoZWlnaHQiOjQzMi42MzI5LCJ3aWR0aCI6NzAwLCJzaXplX2JlaGF2aW9yIjowLCJjb25kaXRpb25zIjpbXX0= -->

<img src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAArwAAAGwCAIAAADE8iHyAAAACXBIWXMAAA7DAAAOwwHHb6hkAAAgAElEQVR4nOzde1hTZ7Y4/pdLCibhdlIuAaVDTIJHBAJm6GRmnDQwGEu4DAKhFA6D04eivZzWWO083tqjjtRWpM+MPVYcnwYxouCBkKCADhyBgmLBr4DESLkqXvAHAbFJSEzI74+MGY4EipotEtbnr+yVd6+99m6rqzvv3q+NwWBAAAAAAAA/x3auCwAAAADA/ABNAwAAAABmBZoGAAAAAMwKNA0AAAAAmBX7uS4AzFeaWj4WaSeGsciK7KnLsEhbsucBFmnJuP8Pi7RYYO19c65LsF4PfsQi68SI3PJJbWwdfnfS8mnByweaBvCM9P3FmKS9hUVWZOP0WyzSdhQPYpFW69CPRVosvL7FZ65LsFqG+5ewSDtx93vLJ7Wxs3xO8FKCnycAAAAAMCvQNAAAAABgVqBpAAAAgKFr166x2WxnZ2cajXbs2DFj8ObNm1wu19XVlc1md3Z2GoN8Pt9mktu3byOELl68GBYWRiAQQkJCampqjCPLy8sDAgIIBEJoaOiFCxeMwbq6urCwMGdnZzabfePGjRd9ngsDNA0AAACwolKpOBwOm83u7+/Pzc3NyspqaWlBCKWlpVGpVLlczmKxkpKSjIO7u7tFIlHfY15eXhqNJi4ujsvldnZ2Jicnx8fHP3jw4P79+6mpqVu3bu3r6+PxeAkJCSqV6tatWzExMQKBYGBgIDExMTIyUqfTYXFGLBart7cXi8zzgg28Rho8G9UxGyzSYjQREheGyUTIvasxmQhJmz8TIflN/znXJVit+TUREv8f5v+G/sc//pGenn7nzh3jZkJCQkhICJ/PDwwMHB4eJhKJOp2ORCJVV1czmUxXV9f29vYlS5aYdu/s7Fy2bNmDBw+cnJwePXrk5ORUX19/+/btvXv3Xr58GSGkVqtdXFyam5sbGhpKSkrOnz9v3DE4OPgvf/lLdHS0xc/V3d39+++/9/f3t3jmeQHuNAAAAMCKr6/vgQMHTJtDQ0NkMvn69esBAQFEIhEhZG9vHxoaKpPJhoeHx8bGMjMziURiYGDg2bNnEUIUCoVMJh85ckStVguFQicnp+XLl0dGRhq/RQjJ5XJbW1sqlTo+Pj4yMmI60KJFi378cdpnVnfu3Ont7e3k5BQdHT15ryeUlZUFBQXh8XgKhSIUChFCXC5XoVBEREQ0NDQ876WZn+CRS/CkPXv2aLVaAoEw9auIiAgmk/niSwIAzFN0Op1Opxs/5+fnt7a2RkVFSSQSEolkGkMikQYHB7u7u3E4XHp6+okTJ06fPr127dqOjo6lS5ceP348PDx88+bNExMTEonE+EcTgUDo7e1du3bt9evXT506hcfjV69evWXLFqlUGh4eXlRU1NTUFBERYbakysrKQ4cOXbhwwc3NLSkpKTc3d9euXVOHjY6O8vl8oVAYFRVVUVGRkZGRmppaVVXl7u5eXV29YO80QNMAnlRaWrp06VIKhTL1K/gxCwDwDMbGxgQCgVgsPnfuHJlMnpiYeGKATqcLCwvTaDTGzXfffbewsLC0tPStt95KTU09efIkj8erqanJzMysq6szdiGenp7Z2dlisVggELz++usBAQHHjx/fuHHjvXv3eDxeZGSkq6ur2WJ0Op1er7958yadTq+trZ2uZgKBIJfL/fz8tFotiUTSaDQqlcrFxcVCl2S+gqYBPMnDw+NPf/rTmjVr5roQAIA1kMlkMTExwcHBra2tPj4+CCFPT0+FQmEaMDIy4uXl9cRe/v7+d+/eLS8vDw4OTk5ORgjFxsYaO4kNGzbY2toSicQ1a9asWbNm2bJlVVVVf/zjH5OTk40jEUKrVq3y9PQ0W090dPTnn3++ffv2pKQkNpv91VdfLV++fOowOzs7kUgkFArxeHxgYKAlroQ1gDkNAAAAsDI+Ps7j8dLS0kpKSowdA0IoKChIJpOp1WqEkF6vv3LlSmBg4LFjx7Kyskw79vT0LF26dHx83MbmX3OuDQaDRqP58ssv33vvPVPQwcFBqVTK5fKsrCzjPYyHDx82NzdzOByzJfX29vJ4vJaWllu3blEoFIFAYHaYVCrNy8urra1ta2szPSkKoGkAAACAFYlEotPp1q9fP/iYUqmkUqlMJnPv3r1arTYnJ8fX15fJZK5YseLIkSNHjx4dHR3Nz89vaWlJTk7m8Xh1dXXFxcVKpVIqlZaXl8fGxvJ4vLKysgsXLqhUqqKiohs3bqxevXrx4sXFxcU5OTmDg4MffPABh8OZ/BTGZJWVlX/4wx9u3rxpfBuEVqs1O0yhUDg4ODg6OqrV6h07diCEjF0OQkipVGJ0uV5+0DQAAADASnt7+8DAgLe3t9dj+/fvRwgVFRVdunSJTCZXVFRIJBKEUGhoqFgsPnjw4JIlSw4fPlxVVUUikWg0WklJSXZ2toeHx7Zt20QiEYPBYLFYhw8f3rhxo4eHx759+8RiMZVKJRKJp06d+vvf/758+fJHjx6dOnVqupLeeecdJpMZEhLi6+srl8vz8vLMDktJSWEwGBQKJSQkhEajxcXFxcbGIoQSExM5HE59fT02F+xlB+9pAE968803P/roo5+d0wDvaUDwngZ4TwOWrOM9DcDKwERIAAAAVujq1as5OTlT4+np6ZGRkU87DBhB0wAAAMAKMRiMgoICSw0DRtA0LFAff/yxVCp1c3Ob+tW1a9d++9vfwiOXAAAAngBNwwK1aNGiyMjIzMzMqV+9//77ISEhL74kI4wmH2BkHk0+QAj9qHltrkuwTvNp8gFCtmQM/hOzgTn1CwU0DQvUK6+84u3tvXLlyqlfubm52drCHwEAAACeBH83AAAAAGBWoGkAAACAoWvXrrHZbGdnZxqNZnq1otkgn8+3meT27duT82i1Wr1ePzly586dnp4e0+bMuwOLgKYBAAAAVlQqFYfDYbPZ/f39ubm5WVlZLS0tZoMIoe7ubpFI1PfY5AUpVCpVUFCQWCyenHznzp2lpaWmzRl2nw0Wi9Xb2/t8p2v9YE4DAAAArDQ2NuJwOOPa09HR0cZlpkdGRqYGV65c2d3dvWrVKrOvfxYIBJ2dnabNwsJCkUh05swZ4/sljWbYfTa6urqme6U0MIE7DQAAALDi6+t74MAB0+bQ0BCZTDYbHB4eHhsby8zMJBKJgYGBZ8+eNQ2QSqUdHR0MBsMUIZFI8fHxk6dyz7D7VGVlZUFBQXg8nkKhCIVChBCXy1UoFBEREQ0NDZY4b6sFdxrAk3Q63b179yb/Umji5eWFx+NffEkAgHmKTqfT6XTj5/z8/NbW1qioKDKZPDXY3d2Nw+HS09NPnDhx+vTptWvXdnR0LF26dHBwcOPGjefPn+fz+aa0q1evRghVVlaaItPtPrWk0dFRPp8vFAqNdzgyMjJSU1Orqqrc3d2rq6v9/f0xvBzzHzQN4Ek//fTTjh07XnnllSfitra227Zty8jImIuiAADz2NjYmEAgEIvF586dI5PJZoNkMlmj0Ri/evfddwsLC0tLSz/55JN169Zt377dz89v5kOEhYWZ3X3qSAKBIJfL/fz8tFotiUTSaDQqlcrFxcVyp2vNoGkAT7p48eJclwAAsB4ymSwmJiY4OLi1tdXHx2eG4GT+/v5379799ttvcThcWlqaXq83GAx6vX5iYmI2L5Ix7m72Kzs7O5FIJBQK8Xh8YGDgc5zZQgRzGgAAAGBlfHycx+OlpaWVlJSYmgOzwWPHjmVlZZl27OnpWbp0aUNDg0QiweFw9vb2LS0tycnJCQkJZg9kdnezI6VSaV5eXm1tbVtbm+lpTzBL0DQAAADAikQi0el069evH3xMqVSaDa5YseLIkSNHjx4dHR3Nz883tggFBQWGx5hM5unTpyc/YzmZ2d3NjlQoFA4ODo6Ojmq1eseOHQghtVpt/EqpVGJ0HawG/DxhtZqbm6OjoxcvXmz2256enoCAgJ07d77gqgAAC0p7e/vAwIC3t7cp8vnnn+t0uqnBzz77TCwWf/bZZx9//HFgYGBVVRWJRJr9gUJDQ2e5e0pKytmzZykUCplM/vTTT+Pi4mJjYy9fvpyYmMjhcMrLy1etWvXM52v1bAwGw1zXADBx+fLlP/3pT/n5+Wa/PXjw4KuvvvrVV189c37VMZtn3ncGmKymg5nimMtzXcJTwGLBqm1X4yyec96BBauQja3j72stnxa8fOBOg9WysbHB4/Fml6RCCPn6+trZ2b3gkgAA4IW5evVqTk7O1Hh6enpkZOSLr8c6QNMAAADACjEYjIKCgrmuwtrAREgAAAAAzArcaQAvl/n0Oy5msJh8gBCiOfRbPCdGP+fPLzYev8IirR02aeEfGXgecKcBAAAAALMCTQMAAAAAZgWaBgAAAADMCjQNAAAAMHTt2jU2m+3s7Eyj0UyvbTYb1Gq1GzZs8PLyolKpUqnUGOTz+TaT3L5925T5zp07pvV4v/nmG5v/Kz4+/gWe5UIBTQMAAACsqFQqDofDZrP7+/tzc3OzsrJaWlrMBhFCGzZsGBgYuHz58qZNm95+++2hoSGEUHd3t0gk6nvMy8vLlHznzp2mt0qnp6f3TcJisSYvpT0bLBart7fXcqduneDpCQAAAFhpbGzE4XC7du1CCEVHR0dFRVVUVIyMjEwN0un0kydP3rhxY/HixRs2bGhra+vr63v11Ve7u7tXrVq1ZMmSyWkLCwtFItGZM2f2799vjDg5OTk5ORk/l5eXEwiElJSUpyq1q6tLq9Va4JytGtxpAAAAgBVfX98DBw6YNoeGhshkstlgY2Ojr6+vabmcQ4cOMZnM4eHhsbGxzMxMIpEYGBh49uxZ47ckEik+Pt7sG291Ot2mTZu+/vrrGaoqKysLCgrC4/EUCkUoFCKEuFyuQqGIiIhoaGh47pO2ZnCnYYG6cePG999/L5FInmqvjz/+OC0tDaOSAADWh06n0+l04+f8/PzW1taoqCgymTw1eP78eTc3t8zMzJKSEjc3t+3bt2dkZHR3d+NwuPT09BMnTpw+fXrt2rUdHR1Lly5dvXo1QqiysnLqEUUi0YoVKwICAqYraXR0lM/nC4VC4x2OjIyM1NTUqqoqd3f36upqf39/DC6D9YCmYYHKycnp6enB4/FPtdd069MDAMAMxsbGBAKBWCw+d+4cmUw2GxwdHb148WJiYmJfX98PP/ywZs0af39/Foul0WiM4999993CwsLS0tJPPvlkugMZDIbs7OzvvvtuhmIIBIJcLvfz89NqtSQSSaPRqFQqFxcXC56vFYOmYYHy8fHx8fGZ6yoAANZPJpPFxMQEBwe3traa/tiZGnRzc6PT6QKBACEUHh7+xhtvnDlzhsViTU7l7+9/9+7dGY7V2Nio1Wqf2OsJdnZ2IpFIKBTi8fjAwMDnO7kFB+Y0AAAAwMr4+DiPx0tLSyspKTF1DGaDVCr10aNHph3d3NwWLVp07NixrKwsU7Cnp2fm+50nTpxISkqauSSpVJqXl1dbW9vW1mZ62hPMEjQNAAAAsCKRSHQ63fr16wcfUyqVZoMsFotIJGZnZz98+PDMmTNnz56Nj49fsWLFkSNHjh49Ojo6mp+f39LSkpycPMPhqqqq2Gz2zCUpFAoHBwdHR0e1Wr1jxw6EkFqtNn6lVCotdeLWCpoGAAAAWGlvbx8YGPD29vZ6bP/+/WaDCKHy8vLq6mofH5+tW7eePHly+fLloaGhYrH44MGDS5YsOXz4cFVVFYlEmu5YAwMD3d3dYWFhM5eUkpLCYDAoFEpISAiNRouLi4uNjUUIJSYmcjic+vp6y14BK2NjMBjmugbw8+7evVtTU+Pp6Tn7Xa5fv/63v/2ts7MTo5JUx2wwyowFjFa5LI65jEXaebTKZZL0Z/6AXggwWuUSI5iscmlj6/j7WsunBS8fmAg5P3z33XcHDx6c4SGiqYaHh+/fv49dSQAA8DK7evVqTk7O1Hh6enpkZOSLr8c6QNMwPyxZsuT3v//9U83Z+eGHH95//33sSgIAgJcZg8EoKCiY6yqsDTQN4OWiv4VR2u+xSMv+JRZZEfrB8r8jYMR+eTomaV2X/Pygl8bDr9+c6xLmmo2d4+/nugbwQsBESAAAAADMCjQNAAAAAJgVaBoAAAAAMCvQNAAAAMDcnTt3enp6TJt8Pt9mktu3b08XvHjxYlhYGIFACAkJqampeSKtVqvV6/XGz+Xl5QEBAQQCITQ09MKFC8bgzLuDpwVNAwAAAMzt3LmztLTUtNnd3S0Sifoe8/LyMhvUaDRxcXFcLrezszM5OTk+Pv7BgwemJCqVKigoSCwWI4Tu37+fmpq6devWvr4+Ho+XkJCgUqlm3t2yWCxWb28vRslfHtA0AAAAwFBhYWF0dPTRo0cnB7u7u1etWvXaY3Z2dmaD/f39Q0NDW7Zs8fHx2bRpk0ajmfzCOoFAYNpsbGz09/dPTU11d3ffunXrw4cPu7q6Zt7dsrq6urRaLUbJXx7QNAAAAMAQiUSKj49fuXKlKTI8PDw2NpaZmUkkEgMDA8+ePTtdkEKhkMnkI0eOqNVqoVDo5OS0fPlyYxKpVNrR0cFgMIybkZGRxl0QQnK53NbWlkqlzrD7VGVlZUFBQXg8nkKhCIXCGc5o6kgul6tQKCIiIhoaGp7jUs0D8J6GBaqxsbGoqMi0sP0s8Xi8FStWYFQSAMAqrV69GiFUWVlpinR3d+NwuPT09BMnTpw+fXrt2rUdHR3Dw8NTg0uXLj1+/Hh4ePjmzZsnJiYkEgmBQEAIDQ4Obty48fz583w+35iTQCAQCITe3t61a9dev3791KlTeDweIWR296lGR0f5fL5QKIyKiqqoqMjIyEhNTcXhcLMcWVVV5e7uXl1d7e/vb/EL+FKBpmGBOnv2bH19/dO+S3VkZASjegAAC0dYWJhGozF+fvfddwsLC0tLSz/55JOpwbfeeis1NfXkyZM8Hq+mpiYzM7Ouro5Op69bt2779u1+fn5PZPb09MzOzhaLxQKB4PXXX9fpdGZ3n1oSgUCQy+V+fn5arZZEImk0GpVK5eLi8jwjrRI0DQuUl5fXr3/96y+++GKuCwEALHT+/v537941GywvLw8ODjYuhx0bG2vsJFxcXHA4XFpaml6vNxgMer1+YmLip59+srW1JRKJa9asWbNmzbJly6qqqtRq9dTdP/3006k12NnZiUQioVCIx+MDAwNnqHb2I60SzGkAAADwQh07diwrK8u02dPTs3TpUrPB8fFxG5t/LahrMBg0Gk1DQ4NEIsHhcPb29i0tLcnJyQkJCV9++eV7771nGung4KBUKs3ubrYkqVSal5dXW1vb1tY28yo/sx9plaBpAAAA8EKtWLHiyJEjR48eHR0dzc/PN/7FbzbI4/Hq6uqKi4uVSqVUKi0vL4+NjS0oKDA8xmQyT58+XVpayuPxysrKLly4oFKpioqKbty4sXr1arO7my1JoVA4ODg4Ojqq1eodO3YghNRq9dOOVCqV2Fywlwg0DQAAAF6o0NBQsVh88ODBJUuWHD58uKqqikQimQ3SaLSSkpLs7GwPD49t27aJRCLT4xJPYLFYhw8f3rhxo4eHx759+8RiMZVKnf3uKSkpDAaDQqGEhITQaLS4uLjp2ovpRiYmJnI4nPr6ektdpZeTjcFgmOsawM8rKCg4f/78MyyNffnyZbPfHjx48MaNG3/729+euSTVMZufH/T0MFrlEiMjlT8/5hnU/vAKJnkxkNJ9EIu0sMrlPGNj57RVN9dFgBcBJkICAABYKK5evZqTkzM1np6e/sTTZLMfuaBA0wAAAGChYDAYBQUFlh25oEDTYLVGRkbkcjmTyTT77a1bt9zd3V9wSQAAAOY1aBqsVnh4+KlTpzw8PMx+KxKJhoaGXnBJs2GHzW/ZtuTfYpG2drf5KSMvpx81r1k8p34Ikzf562SYPMlmvzwdi7S4MEz+7cKIjcevMEgKc+oXCmgarJa9vf2bb047P+vixYuPHj16kfUAAACY76A9BAAAAMCsQNMAAAAAQ9euXWOz2c7OzjQazfTcuNkgn8+3meT27dsIoZs3b3K5XFdXVzabbVrYWqvVbtiwwcvLi0qlSqVSY9DsSGBZ0DQAAADAikql4nA4bDa7v78/Nzc3KyurpaXFbBAh1N3dLRKJ+h7z8vJCCKWlpVGpVLlczmKxkpKSjGk3bNgwMDBw+fLlTZs2vf3228YZWmZHzh6Lxert7bXo2VshmNMAAAAAK42NjTgcbteuXQih6Oho43LSIyMjU4MrV67s7u5etWrVkiX/mg7d2dnZ1NR09uxZIpG4Z8+eQ4cONTc3+/v7nzx58saNG4sXL96wYUNbW1tfX59CoZg6crrHx8zq6urSarUWvwJWBpoGy7t06ZLZVdifR3d3t0qlsmxOAADAmq+v74EDB0ybQ0NDZDLZbHB4eHhsbCwzM/P777/38/Pbt29fVFTU9evXAwICiEQiQsje3j40NFQmkw0PD/v6+i5evNi4+6FDhxBCZWVlU0dO1zSUlZXt2LGjq6vLy8tr586dGRkZXC5XoVBEREScOnXqN7/5DXYXZL6DpsHCFArF7373u6CgIMumvXPnjr29Jf9hGVd7GxkZeaq9nJycLFsGAMC60el0Op1u/Jyfn9/a2hoVFUUmk6cGu7u7cThcenr6iRMnTp8+vXbt2o6Ojnv37pFIJFM2Eok0ODiIEHJzc8vMzCwpKXFzc9u+fXtGRsZ0I6caHR3l8/lCodB4hyMjIyM1NbWqqsrd3b26utrf3x+ra2EV4C8Ay3NycmpubrZsTuPaExZMeP/+fZFIdPr06afaa9++fZmZmRYsAwCwEIyNjQkEArFYfO7cOTKZbDZIJpNN61a/++67hYWFpaWlBALhiVQ6nW50dPTixYuJiYl9fX0//PDDmjVr/P39JyYmpo40WwyBQJDL5X5+flqtlkQiaTQalUrl4uJi0TO2WtA0LFC7d+/evXv3XFcBALB+MpksJiYmODi4tbXVx8dnhuBk/v7+d+/e/c1vfqNQKEzBkZERLy8ve3t7Op0uEAgQQuHh4W+88caZM2dCQ0OnjjRbj52dnUgkEgqFeDw+MDDQcie6IMDTEwAAALAyPj7O4/HS0tJKSkpMzYHZ4LFjx7Kyskw79vT0LF26NCgoSCaTqdVqhJBer79y5UpgYCCVSp38bjo3N7dFixaZHWm2JKlUmpeXV1tb29bW9lRLBwMETQMAAADsSCQSnU63fv36wceUSqXZ4IoVK44cOXL06NHR0dH8/PyWlpbk5GQqlcpkMvfu3avVanNycnx9fZlMJovFIhKJ2dnZDx8+PHPmzNmzZ+Pj482ONFuSQqFwcHBwdHRUq9U7duxACBlbDYSQUql8cZdmfoKmAQAAAFba29sHBga8vb29Htu/f7/ZYGhoqFgsPnjw4JIlSw4fPlxVVWWc2FhUVHTp0iUymVxRUSGRSIxpy8vLq6urfXx8tm7devLkyeXLl083cqqUlBQGg0GhUEJCQmg0WlxcXGxsLEIoMTGRw+HU19e/kAszX9kYDIa5rsGqKBQKGo02PDxs2bTGiZAv1Z001TGbuS7hKWC0YFVxzEJfsGrb1TiL50QIGe5fwiItRgtWYbS8FkYwWrDKIWif5dOClw9MhAQAAGCFrl69mpOTMzWenp4eGRn54uuxDtA0AAAAsEIMBqOgoGCuq7A2MKcBAAAAALMCdxrAMxo6MtcVPJ3vsUjKS8UiK/oJm/X52OhHi+fUlu+3eE6E0KMeLLIiHAWTfw2ctmEyAwMpn+6FrXNqPs1wAs8D7jQAAAAAYFagaQAAAADArEDTAAAAAIBZgaYBAAAAhq5du8Zms52dnWk0mullMzdv3uRyua6urmw2u7PzyVk8Wq1Wr9cbP/P5fJtJbt++jRC6ePFiWFgYgUAICQmpqakxjjQbBJYFTQMAAACsqFQqDofDZrP7+/tzc3OzsrJaWloQQmlpaVQqVS6Xs1ispKSkJ3YJCgoSi8XGze7ubpFI1PeYl5eXRqOJi4vjcrmdnZ3Jycnx8fEPHjwwG8TijFgsVm9vLxaZ5wVoGgAAAGClsbERh8Pt2rXLzc0tOjo6KiqqoqKis7Ozqalp3759Xl5ee/bs6evra25uNu0iEAgm33vo7u5etWrVa4/Z2dn19/cPDQ1t2bLFx8dn06ZNGo2ms7PTbBCLM+rq6tJqtVhknhegaQAAAIAVX1/fAwcOmDaHhobIZPL169cDAgKIRCJCyN7ePjQ0VCaTGQdIpdKOjg4Gg2HcHB4eHhsby8zMJBKJgYGBZ8+eRQhRKBQymXzkyBG1Wi0UCp2cnJYvX242OF1VO3fu9Pb2dnJyio6OHhmZ9tHWsrKyoKAgPB5PoVCEQiFCiMvlKhSKiIiIhoYGC1ydeQje07BA7dmzp7i42MPD46n2yszM5PP5GJUEALA+dDqdTqcbP+fn57e2tkZFRUkkEuNiVEYkEmlwcBAhNDg4uHHjxvPnz5v+nOnu7sbhcOnp6SdOnDh9+vTatWs7OjqWLl16/Pjx8PDwzZs3T0xMSCQSAoGAEDIbnKqysvLQoUMXLlxwc3NLSkrKzc3dtWvX1GGjo6N8Pl8oFBrvjmRkZKSmplZVVbm7u1dXV/v7+1v2Qs0X0DQsUFqtNiAg4E9/+tPsd7GxsTG1/wAAMHtjY2MCgUAsFp87d45MJk9MTDwxQKfTIYTWrVu3fft2Pz8/UzwsLEyj0Rg/v/vuu4WFhaWlpW+99VZqaurJkyd5PF5NTU1mZmZdXR0ej58aNPUrTxxLr9ffvHmTTqfX1tZOVzOBQJDL5X5+flqtlkQiaTQalUrl4uJigcsxn0HTsEAtWrToF7/4xe9///u5LgQAYOVkMllMTExwcHBra6uPjw9CyNPTU6FQmAaMjIx4eXl9++23OBwuLS1Nr9cbDAa9Xj8xMWFr+39+Q/f397979255eXlwcHBycjJCKDY21thJuLi4THaOOEYAACAASURBVA1++umnU+uJjo7+/PPPt2/fnpSUxGazv/rqK7M/ZNjZ2YlEIqFQiMfjAwMDLXpJ5jGY0wAAAAAr4+PjPB4vLS2tpKTE2DEghIKCgmQymVqtRgjp9forV64EBgY2NDRIJBIcDmdvb9/S0pKcnJyQkHDs2LGsrCxTtp6enqVLl46Pj9vY/OvF1QaDQaPRmA2aLam3t5fH47W0tNy6dYtCoQgEArPDpFJpXl5ebW1tW1ub6UlRAE0DAAAArEgkEp1Ot379+sHHlEollUplMpl79+7VarU5OTm+vr5MJrOgoMDwGJPJPH36dGlp6YoVK44cOXL06NHR0dH8/HxjM8Hj8erq6oqLi5VKpVQqLS8vj42NNRs0W1JlZeUf/vCHmzdvGl/8MN2jEAqFwsHBwdHRUa1W79ixAyFk7HIQQkqlEqPL9fKDpgEAAABW2tvbBwYGvL29vR7bv38/QqioqOjSpUtkMrmiokIikUy3e2hoqFgsPnjw4JIlSw4fPlxVVUUikWg0WklJSXZ2toeHx7Zt20QiEYPBMBs0m/Odd95hMpkhISG+vr5yuTwvL8/ssJSUFAaDQaFQQkJCaDRaXFycsQtJTEzkcDj19fWWuDzzj43BYJjrGqyKQqGg0WjDw8OWTVtQUHD+/HkL3iLLzs5++PDh3r17nznDzVWwrh0impllZQEYrXKJBbc1mKTFbJVLTNLCKpcI2SACd65rAC8CTIQEAABgha5evZqTkzM1np6eHhkZ+bTDgBE0DfODwWDQarUzvITkaSmVStOr3QEAwPowGIyCggJLDQNG0DTMD0NDQ2VlZefOnbNUQqVSuWzZsn379j1zBozuzM8v8+h3BIxg9DsCRrCqdj79jgDAc4GmYX4QCATTPRf0bIxzGiyYEAAAgNWDpycAAAAAMCvQNAAAAABgVqBpAAAAgKFr166x2WxnZ2cajfbEc+N37tzp6TEz00Sr1Zpmat+8eZPL5bq6urLZbNNq12ZzlpeXBwQEEAiE0NDQCxcuYHhKCxg0DQAAALCiUqk4HA6bze7v78/Nzc3KymppaTF9u3PnztLS0qm7BAUFicVi42ZaWhqVSpXL5SwWKykpabqc9+/fT01N3bp1a19fH4/HS0hIUKlUWJwRi8Xq7e3FIvO8AE0DAAAArDQ2NuJwuF27drm5uUVHRxuXmUYIFRYWRkdHHz16dOouAoHAdEehs7Ozqalp3759Xl5ee/bs6evra25uNpuzsbHR398/NTXV3d1969atDx8+7OrqwuKMurq6pnvz9EIATQMAAACs+Pr6HjhwwLQ5NDREJpMRQiQSKT4+fuXKlU+Ml0qlHR0dpjdAX79+PSAggEgkIoTs7e1DQ0NlMpnZnJGRkWfPnjVG5HK5ra0tlUo1W9Lp06f9/f2Ny1lt27YtISFhuuLLysqCgoLweDyFQhEKhQghLperUCgiIiIaGhqe+lpYBXjkcoEaGRm5du3adC9dn87vfve7ZcuWYVQSAMD60Ol0Ov2fL3XJz89vbW2NiopCCK1evRohVFlZOXnw4ODgxo0bz58/z+fzjZF79+6RSCTTABKJNDg4aDYngUAgEAi9vb1r1669fv36qVOn8Hi82ZISExNPnjy5e/fut95667vvvvt//+//mR02OjrK5/OFQqHxTkZGRkZqampVVZW7u3t1dbW/v/9zXZd5C5qGBYpGo12/fn3yj4uz4efnB00DAOBpjY2NCQQCsVh87tw5450Gs9atW7d9+3Y/Pz9TZGJi4okxOp1uhpyenp7Z2dlisVggELz++uteXl5mD/TNN98EBweLxeLc3FxPT0+zYwgEglwu9/Pz02q1JBJJo9GoVCoXF5fZn7VVgqZhgcrMzMzMzJzrKgAA1k8mk8XExAQHB7e2tvr4+Ew37Ntvv8XhcGlpaXq93mAw6PX6iYkJT09PhUJhGjMyMmLsA6bmHBsbs7W1JRKJa9asWbNmzbJly6qqqv74xz+aPZanp2dSUtLp06cTExOnq8fOzk4kEgmFQjweHxgY+Mynb2VgTgMAAACsjI+P83i8tLS0kpKSGToGhFBDQ4NEIsHhcPb29i0tLcnJyQkJCUFBQTKZTK1WI4T0ev2VK1cCAwPN5vzyyy/fe+89UzYHBwelUjndsW7cuPE///M/v/jFL8wuVWUklUrz8vJqa2vb2tosuMLwfAdNAwAAAKxIJBKdTrd+/frBx6b7u7ygoMDwGJPJPH36dGlpKZVKZTKZe/fu1Wq1OTk5vr6+TCbTbE4ej1dWVnbhwgWVSlVUVHTjxg3jtImp9Hp9RkbGrl278vPzs7Ozp3vIQqFQODg4ODo6qtXqHTt2IISMvQtCaIZ2xOpB0wAAAAAr7e3tAwMD3t7eXo/t37//qTIUFRVdunSJTCZXVFRIJJLpcrJYrMOHD2/cuNHDw2Pfvn1isXi6pydycnIcHBzeeecdOp2+adOmzMxMg8EwdVhKSgqDwaBQKCEhITQaLS4uLjY2FiGUmJjI4XDq6+uf/mJYAxuzFws8M4VCQaPRhoeH57oQzCnesZnrEuYerHIJi50ihP7trxVzXcKcs0EE7lzXAF4EmAgJAADACl29etXslIX09PTIyMinHQaM4E6DhcGdhgUF7jTAnQYEdxoQgjsNCwfMaQAAAADArMDPE+DlgqNgkvaRmYX0LMC3cqH/L+bNNW/OdQkvgf/E5CJgdB8Lk5tDtnb/dkSHQV7w0oE7DQAAAACYFWgaAAAAADAr0DQAAADA0LVr19hstrOzM41GM71asby8PCAggEAghIaGXrhwwRjUarUbNmzw8vKiUqlSqXSG3S9evBgWFkYgEEJCQmpqaozBmzdvcrlcV1dXNpttWlwbWBY0DQAAALCiUqk4HA6bze7v78/Nzc3Kymppabl//35qaurWrVv7+vp4PF5CQoJKpUIIbdiwYWBg4PLly5s2bXr77beHhobM7q7RaOLi4rhcbmdnZ3Jycnx8/IMHDxBCaWlpVCpVLpezWKykpKSnLZXFYvX29lr+ElgXeOTSwuCRy+c0vyZCwrN2MBESYfbcqXVMhPzHP/6Rnp5+584d42ZCQkJISMiKFSv27t17+fJlhJBarXZxcWlubvbz8/Py8rpx48bixYsRQhs2bHjnnXdGR0en7s7n85ctW/bgwQMnJ6dHjx45OTnV19e7uLgEBgYODw8TiUSdTkcikaqrq5lM5uxPwt3d/fvvv1+wa17PEtxpAAAAgBVfX98DBw6YNoeGhshkcmRk5NmzZ40RuVxua2tLpVIbGxt9fX2NHQNC6NChQ0wm0+zuFAqFTCYfOXJErVYLhUInJ6fly5dfv349ICCASCQihOzt7UNDQ2Uy2XRVlZWVBQUF4fF4CoUiFAoRQlwuV6FQRERENDQ0YHAZrAc8crng/PTTT48ePXqGHW1sbFxdXS1eDwDAitHpdDr9nzc38vPzW1tbo6KiCAQCgUDo7e1du3bt9evXT506hcfjBwcH3dzcMjMzS0pK3Nzctm/fnpGRYXZ3e3v748ePh4eHb968eWJiQiKREAiEe/fukUgk03FJJNLg4KDZkkZHR/l8vlAojIqKqqioyMjISE1Nraqqcnd3r66uhjsNM4OmYWGRy+VBQUHGZvwZ/OUvf9mwYYNlSwIAWL2xsTGBQCAWi8+dO0cmk41BT0/P7OxssVgsEAhef/310dHRixcvJiYm9vX1/fDDD2vWrPH392exWFN3HxgYSE1NPXnyJI/Hq6mpyczMrKurm5iYeOKgOp35X0wIBIJcLvfz89NqtSQSSaPRqFQqFxcXTK+A1YCmYWGxt7d/7bXXfvzxx7kuBACwUMhkspiYmODg4NbWVh8fH4TQ2NiYra0tkUhcs2bNmjVrli1bVlVV5ebmRqfTBQIBQig8PPyNN944c+YMi8Waunt5eXlwcHBycjJCKDY2trCwsLS0lEajKRQK00FHRka8vLzM1mNnZycSiYRCIR6PDwwMxP4CWBWY0wAAAAAr4+PjPB4vLS2tpKTE+Fc+QujLL7987733TGMcHByUSiWVSp38y6mbm9uiRYvM7j4+Pm5j86+J2AaDQaPRBAUFyWQytVqNENLr9VeuXJmuIZBKpXl5ebW1tW1tbaZnOMEsQdMAAAAAKxKJRKfTrV+/fvAxpVLJ4/HKysouXLigUqmKiopu3LixevVqFotFJBKzs7MfPnx45syZs2fPxsfHT7d7XV1dcXGxUqmUSqXl5eWxsbFUKpXJZO7du1er1ebk5Pj6+k736IRCoXBwcHB0dFSr1Tt27EAIGVsNhJBSqXxxl2Z+gqYBAAAAVtrb2wcGBry9vb0e279/P4vFOnz48MaNGz08PPbt2ycWi6lUKkKovLy8urrax8dn69atJ0+eXL58udndaTRaSUlJdna2h4fHtm3bRCIRg8FACBUVFV26dIlMJldUVEgkkulKSklJYTAYFAolJCSERqPFxcXFxsYihBITEzkcTn19/Qu7OPMRvKfBwl7y9zR0dXW9+eabFpnTAO9pQPCeBnhPA0II3tOAYMGqBQQmQgIAALBCV69ezcnJmRpPT0+PjIx88fVYB2gaAAAAWCEGg1FQUDDXVVgbaBqsR05ODpFInDypeKrBwcHR0dEXVhIAAABrAk2D9di1a1dsbCwej59hjEKhMM0Tfk6ObIukeZIt+bdYpMWFYZEVabqqMcmLDbtXLf9r9t+/p1k8J3ZoDv1YpE3I1GKRFqP/xDCByQQn8DKCpsF6ODk5ffHFF6ZHmc3q6uq6evXqCysJAACANYFHLgEAAAAwK3CnAd2/f/+zzz77xS9+YZFsSqXSuDA8AAAAYGXgTgNqa2uTSCQjFvLgwYOp66YAAMCCde3aNTab7ezsTKPRTK9tvnnzJpfLdXV1ZbPZnZ3/fCUFn8+3meT27dvTBU20Wq1erzd+rqurCwsLc3Z2ZrPZN27ceIGnuIDAnQbk5OS0ZMmSL774wiLZFArF8ePHLZIKAADmO5VKxeFwNmzYIBaLGxoakpKSAgICVq5cmZaWFhgYmJ+f//XXXyclJbW2tiKEuru7RSLRb37zG+O+xhWnzAZNyUNDQ//yl78kJCTcunUrJibm8OHDUVFR+fn5kZGRPT099vZP8Xcci8U6ceKEn5+fxU7eGsGdBgAAAFhpbGzE4XC7du1yc3OLjo6OioqqqKjo7Oxsamrat2+fl5fXnj17+vr6mpubEULd3d2rVq167TE7O7vpgkYCgcB0l6K8vDwsLOytt95ydnb+8MMP3dzcKisrn6rUrq4urRaTB2GsCTQNAAAAsOLr63vgwAHT5tDQEJlMvn79ekBAAJFIRAjZ29uHhobKZLLh4eGxsbHMzEwikRgYGHj27FmEkNmgkVQq7ejoMK46gRAaHx8fGRkxfbto0aIZ3pdfVlYWFBSEx+MpFIpQKEQIcblchUIRERHR0NBg4UtgXeDniYVFrVYPDg7y+fxn2z01NTUuLs6yJQEArBidTqfT//mCkPz8/NbW1qioKIlEQiKRTGNIJNLg4GB3dzcOh0tPTz9x4sTp06fXrl3b0dExPDw8Nbh06dLBwcGNGzeeP3/e9KfZ6tWrt2zZIpVKw8PDi4qKmpqaIiIizJY0OjrK5/OFQqHxtkdGRkZqampVVZW7u3t1dbW/vz/W12Reg6ZhYfn3f/9349Jwz7b7dEvNAgDADMbGxgQCgVgsPnfuHJlMnjpbXKfThYWFaTQa4+a7775bWFhYWlr6ySefmA2uW7du+/btk+cfBAQEHD9+fOPGjffu3ePxeJGRka6urmaLIRAIcrncz89Pq9WSSCSNRqNSqVxcXDA4bysETcPCYm9v//777891FQCABUQmk8XExAQHB7e2thrfPufp6alQKEwDRkZGJk9vNPL39797967Z4LfffovD4dLS0vR6vcFg0Ov1ExMTtra2ycnJycnJxpGrVq3y9PQ0W4+dnZ1IJBIKhXg8PjAw0FKnuUDAnAYAAABYGR8f5/F4aWlpJSUlpvfVBgUFyWQy4yvt9Xr9lStXAgMDjx07lpWVZdqxp6dn6dKlZoMNDQ0SiQSHw9nb27e0tCQnJyckJMjl8qysLOM9jIcPHzY3N3M4HLMlSaXSvLy82tratrY20yOgYJbgTsP8cOrUKa1W6+joOMOYhw8fwisiAAAvFYlEotPp1q9fPzg4aIwQiUQqlcpkMvfu3btjx46vv/7a19eXyWTa2tpmZGSEhYUlJCSUlZW1tLQUFhb29/dPDb733num5St/+ctf/vnPf05ISPjpp5+Ki4upVGp6evqWLVs4HM6SJUvMlqRQKBwcHBwdHdVq9e7duxFCarXa+POEUql8IVdlHoOmYX7Yu3evi4vL1Dt4k42Pj8PLKAEAL5X29vaBgQFvb29T5PPPP//ss8+KiorS09PJZHJQUJBEIkEIhYaGisXizz777OOPPw4MDKyqqiKRSCQSaWrQ7IGIROKpU6c++OCDL774gsvlnjp1arqSUlJSzp49S6FQyGTyp59+GhcXFxsbe/ny5cTERA6HU15evmrVKotfB6thYzAY5rqGOdbU1PTRRx9dunTJItkUCgWNRhseHrZINhMej/f+++9HRUXNMGbx4sVNTU0zL1hlQapjmCxsh9Eqlxix8fjVXJfwFLBY5XLX4q8snhM7WK1ymbfgH+63scP/h26uiwAvAtxpAAAAYIWuXr2ak5MzNZ6enh4ZGfni67EO0DQAAACwQgwGwzT1AVgKNA3g5TJx9/u5LuEp4H+9DYu0mq5qLNLqhzotnvPP+dO+dO956G9hkRXhwjD58at5PSb/0p7roWGRFgu2dnaf/8dcFwFeCHjkEgAAAACzAk0DAAAAAGYFmgYAAAAAzAo0DQAAADB07do1Npvt7OxMo9FMb2C8efMml8t1dXVls9mm5a0vXrwYFhZGIBBCQkJqamoQQt98843N/xUfH2/KfOfOnZ6enicOZzYILAWaBgAAAFhRqVQcDofNZvf39+fm5mZlZbW0tCCE0tLSqFSqXC5nsVhJSUkIIY1GExcXx+VyOzs7k5OT4+PjHzx4kJ6e3jcJi8WavEjvzp07S0tLnzii2aAFsVis3t5e7PK/5KBpAAAAgJXGxkYcDrdr1y43N7fo6GjjatSdnZ1NTU379u3z8vLas2dPX19fc3Nzf3//0NDQli1bfHx8Nm3apNFoOjs7nZycXnusvb2dQCCkpKQghAoLC6Ojo48ePTr5WGaDFtfV1aXVLtzXeUHTAAAAACu+vr4HDhwwbQ4NDZHJ5OvXrwcEBBCJRISQvb19aGioTCYzvtf5yJEjarVaKBQ6OTktX77ctKNOp9u0adPXX39t3CSRSPHx8StXrpx8LLNBs3bu3Ont7e3k5BQdHT0yMjLdsLKysqCgIDweT6FQhEIhQojL5SoUioiIiIaGhqe7ENYC3tNgPTQaTXl5+b/9279hdwgmkzl5AXsAAJgZnU6n0//5/vL8/PzW1taoqCiJRDJ5CQkSiTQ4OGhvb3/8+PHw8PDNmzdPTExIJBICgWAaIxKJVqxYERAQYNxcvXo1QqiysnLyscwGp6qsrDx06NCFCxfc3NySkpJyc3N37do1ddjo6CifzxcKhca7IxkZGampqVVVVe7u7tXV1f7+/s9yOeY/aBqsx4oVK86dO2dnZ4dRfhsbGycnJ2gaAABPa2xsTCAQiMXic+fOkcnkqevx6nS6gYGB1NTUkydP8ni8mpqazMzMuro6Y8NhMBiys7O/++47ixSj0+n0ev3NmzfpdHptbe10wwgEglwu9/Pz02q1JBJJo9GoVCrjYpgLGTQNc6Orq2tkZMTWdrY/D42MjOj1+pnH/O///u9z1wUAABYmk8liYmKCg4NbW1uNK+p5enoqFArTgJGRES8vr/Ly8uDg4OTkZIRQbGxsYWFhaWnpp59+ihBqbGzUarUsFssi9URHR3/++efbt29PSkpis9lfffXV5N9BTOzs7EQikVAoxOPxgYGBFjm0FYCmYW7ExsbqdDpnZ+dZju/o6Lhy5UpMTAymVQEAgGWNj4/zeLz09PT/+q//MgWDgoJkMplarV60aJFer79y5cq+ffu+//57G5t/rZ1rMBg0Go3x84kTJ4xPWFhEb28vj8f7z//8z5GRkZ07dwoEArO/aEil0ry8vIsXL/r4+Oj1+hMnTliqgHkNmoa5sXTp0qysrOjo6FmO5/F4v/zlLzEtCQAALE4ikeh0uvXr1w8ODhojRCKRSqUymcy9e/fu2LHj66+/9vX1ZTKZLi4u27dvLy4ujoqKqqmpKS8v//Of/2zcpaqq6q9//aulSqqsrPzv//7vM2fOODs729jYTPcohEKhcHBwcHR0VKvVu3fvRgip1WrjzxNKpdJSxcw78PQEAAAArLS3tw8MDHh7e3s9tn//foRQUVHRpUuXyGRyRUWFRCJBCNFotJKSkuzsbA8Pj23btolEIgaDgRAaGBjo7u4OCwuzVEnvvPMOk8kMCQnx9fWVy+V5eXlmh6WkpDAYDAqFEhISQqPR4uLiYmNjEUKJiYkcDqe+vt5S9cwvcKcBAAAAVnbv3m383/QnkMnk8+fPPxFcvXq18QmIyRYvXmwwGMwmLy4unmVwsldeeWU2cyodHR0np1q3bp3xw6FDhw4dOvSzu1sraBoAAABYoatXr+bk5EyNp6enR0ZGPu0wYARNAwAAACvEYDAKCgosNQwYwZwGAAAAAMwK3GkAz8iW/Fss0k7c/R6LtBhVqxu9hUVau1fpWKTVD3VaPOf4tK/GeS44CiZpVScx+bcLIzSH/rkuYbZsMHunHHjZwJ0GAAAAAMwKNA0AAAAAmBVoGgAAAAAwK9A0AAAAwNC1a9fYbLazszONRjt27JgxePPmTS6X6+rqymazOzv/Odvm4sWLYWFhBAIhJCSkpqYGIfTNN9/Y/F/x8fGmzHfu3Onp6Zn5QMCyoGkAAACAFZVKxeFw2Gx2f39/bm5uVlZWS0sLQigtLY1KpcrlchaLZVxXQqPRxMXFcbnczs7O5OTk+Pj4Bw8epKen903CYrH4fL4p+c6dO0tLS2c+0OyxWKze3l7Lnbp1gqYBAAAAVhobG3E43K5du9zc3KKjo6OioioqKjo7O5uamvbt2+fl5bVnz56+vr7m5ub+/v6hoaEtW7b4+Phs2rRJo9F0dnY6OTm99lh7ezuBQEhJSUEIFRYWRkdHHz16dOYDPVWpXV1d061DAUygaQAAAIAVX1/fAwcOmDaHhobIZPL169cDAgKIRCJCyN7ePjQ0VCaTUSgUMpl85MgRtVotFAqdnJwmr1it0+k2bdr09ddfGzdJJFJ8fPzKlStnPtB0VZWVlQUFBeHxeAqFIhQKEUJcLlehUERERDQ0NFjs5K0RvKdhfvjxxx95PN5cV4H++te/fvjhh3NdBQBg3qDT6XT6P987kp+f39raGhUVJZFISCSSaQyJRBocHLS3tz9+/Hh4ePjmzZsnJiYkEgmBQDCNEYlEK1asCAgIMG4al6iYvKS12QOZLWl0dJTP5wuFQuPdiIyMjNTU1KqqKnd39+rqan9/f4teAGsDTcP8YJooBAAA887Y2JhAIBCLxefOnSOTyRMTE08M0Ol0AwMDqampJ0+e5PF4NTU1mZmZdXV1xj7AYDBkZ2fPZpWpJw5kdgyBQJDL5X5+flqtlkQiaTQalUplXPMa/Cz4eQIAAACGZDJZSEiIQqFobW01rnDt6empUChMA0ZGRry8vMrLy4ODg5OTk4lEYmxsLIfDMU1ybGxs1Gq1LBbraQ9klp2dnUgkolKpTCbT+NsEmD1oGgAAAGBlfHycx+OlpaWVlJT4+PgYg0FBQTKZTK1WI4T0ev2VK1cCAwPHx8dtbGxMOxoMBo1GY/x84sQJ4xMWT3sgs6RSaV5eXm1tbVtbGzyZ+bSgaQAAAIAViUSi0+nWr18/+JhSqTT+X/7evXu1Wm1OTo6vry+TyeTxeHV1dcXFxUqlUiqVlpeXx8bGGpNUVVWx2exnOJDZkQqFwsHBwdHRUa1W79ixAyFkbF8QQtPtAkygaQAAAICV9vb2gYEBb29vr8f279+PECoqKrp06RKZTK6oqJBIJAghGo1WUlKSnZ3t4eGxbds2kUjEYDAQQgMDA93d3TP83DDzgaZKSUlhMBgUCiUkJIRGo8XFxRm7k8TERA6HU19fb+FLYF1sDAbDXNcwx5qamj766KNLly5ZJJtCoaDRaMPDwzMPi4mJycrKio6OtshB58T4+VVYpJ1fq1zaL0/HIi1GsFjlUvlX838uPyeMVrl81PPzY56B7AImaXtvv4JJXgzY2NmlKVVzXQV4EeDpCQAAAFbo6tWrOTk5U+Pp6emRkZEvvh7rAE0DAAAAK8RgMAoKCua6CmsDTQN4Ro8uY/I7AmYwqdY+EpOb80g5gkVWLH6e+AmbF4jcvIBJ2uVvzKe0vp2YvNKYSMcgqa0dBknBywgmQgIAAABgVqBpAAAAAMCsQNMAAAAAgFmBpgEAAACGrl27xmaznZ2daTTaE29gvHPnTk/Pvx6EvXnzJpfLdXV1ZbPZpgV36urqwsLCnJ2d2Wz2jRs3pgt+8803Nv9XfHz8izrFBQSaBgAAAFhRqVQcDofNZvf39+fm5mZlZbW0tJi+3blzp2mBCYRQWloalUqVy+UsFsv43uhbt27FxMQIBIKBgYHExMTIyEidTmc2mJ6e3jcJi8Xi8/lPVSqLxert7bXUiVsraBoAAABgpbGxEYfD7dq1y83NLTo62rgaNUKosLAwOjr66NGjppGdnZ1NTU379u3z8vLas2dPX19fc3NzeXl5WFjYW2+95ezs/OGHH7q5uVVWVpoNOjk5vfZYe3s7gUBISUl5qlK7urq0WkyeWLEm0DQAAADAiq+v74EDB0ybQ0NDxhWrSSRSfHz8ypUrTV9dv349ICCAEhcikgAAIABJREFUSCQihOzt7UNDQ2Uy2fj4+MjIv55AXrRo0Y8//mg2aNrU6XSbNm36+uuvZ6iqrKwsKCgIj8dTKBTjQpdcLlehUERERDQ0NDz/WVsxeE/DQlFfX//JJ5+89tprz5zBxsZm3bp1a9assWBVAADrRqfT6fR/vhoiPz+/tbU1KioKIbR69WqEUGVlpWnkvXv3SCSSaZNEIg0ODkZFRW3ZskUqlYaHhxcVFTU1NUVERLz99ttTg6YdRSLRihUrAgICpitpdHSUz+cLhULjbY+MjIzU1NSqqip3d/fq6mp/f3+LXwRrAk3DQnH//n2dTvezy8vOzPQfPwAAzN7Y2JhAIBCLxefOnTPeaZhqYmLiiYhOpwsICDh+/PjGjRvv3bvH4/EiIyNdXV3NBo27GAyG7Ozs7777boZiCASCXC738/PTarUkEkmj0ahUKhcXF4ucqdWDpmGhcHR0JJPJz9k0AADA05LJZDExMcHBwa2trT4+PtMN8/T0VCgUps2RkREvLy+EUHJycnJysjG4atUqT0/P6YIIocbGRq1Wy2KxZqjHzs5OJBIJhUI8Hh8YGPhc57bwwJwGAAAAWBkfH+fxeGlpaSUlJTN0DAihoKAgmUymVqsRQnq9/sqVK4GBgXK5PCsry3gT4uHDh83NzRwOx2zQmOTEiRM/+79GUqk0Ly+vtra2ra3tiUdAwc+CpgEAAABWJBKJTqdbv3794GNKpdLsSCqVymQy9+7dq9Vqc3JyfH19mUzm4sWLi4uLc3JyBgcHP/jgAw6Hs2TJErNBY5Kqqio2mz1zSQqFwsHBwdHRUa1W79ixAyFk7FQQQtPVBkygaQAAAICV9vb2gYEBb29vr8f27592mbeioqJLly6RyeSKigqJRIIQIhKJp06d+vvf/758+fJHjx6dOnVquiBCaGBgoLu7OywsbOaSUlJSGAwGhUIJCQmh0WhxcXGxsbEIocTERA6HU19fb7GTt0YwpwEAAABWdu/evXv37um+LS4unrxJJpPPnz//xJjIyEjTiyBnDi5evNhgMPxsSY6OjpOPu27dOuOHQ4cOHTp06Gd3X+CgaQAAAGCFrl69mpOTMzWenp4eGRn54uuxDtA0AAAAsEIMBqOgoGCuq7A20DRY2MTExE8//fSzbWxTU9OvfvWr6OhoSx2XSCT+7BQeeBAZAADA84CmwcJeffXVvLy8mZ8sQght3rx56dKlFjwuHo/v7e11d3efbsCZM2cs+3Od08cVFsxmovqfN7FIq7+FRVY0fv4TLNI6/nobFmn1rdNOQHvZLH8Dk7RnRK9gkZaXismCBR47f4tF2uKYyxbPaWNnl3bE4lnBywiaBsv74x//+LNjFi9ebHzFOgAAADBfwCOXAAAAAJgVaBoAAABg6Nq1a2w229nZmUajmd7AaDZ48eLFsLAwAoEQEhJSU1PzRB6tVqvX642f6+rqwsLCnJ2d2Wy26dlLs0FgWdA0AAAAwIpKpeJwOGw2u7+/Pzc3Nysrq6WlxWxQo9HExcVxudzOzs7k5OT4+PgHDx5MzhMUFCQWixFCt27diomJEQgEAwMDiYmJkZGROp3ObPCpSmWxWL29vRY+f6sDTQMAAACsNDY24nC4Xbt2ubm5RUdHG1ejNhvs7+8fGhrasmWLj4/Ppk2bNBpNZ2enKY9AIDBtlpeXh4WFvfXWW87Ozh9++KGbm1tlZaXZ4FOV2tXVpdViMqfVmkDTAAAAACu+vr4HDhwwbQ4NDZHJZLNBCoVCJpOPHDmiVquFQqGTk9Py5cuNA6RSaUdHB4PBMG6Oj4+PjIyYdl+0aNGPP/5oNjhdVWVlZUFBQXg8nkKhCIVChBCXy1UoFBEREQ0NDZY5cysFT08sFFqt9t69e0+8tPVp/frXv/7Zp0kBAMCETqfT6XTj5/z8/NbW1qioKDKZPDVob29//Pjx8PDwzZs3T0xMSCQSAoGAEBocHNy4ceP58+f5fL5xl9WrV2/ZskUqlYaHhxcVFTU1NUVERLz99ttTg2ZLGh0d5fP5QqHQeIcjIyMjNTW1qqrK3d29urra398f+6syj0HTsFA4OztPTEw8Z9PwyiuvQNMAAHhaY2NjAoFALBafO3eOTCabDQ4MDKSmpp48eZLH49XU1GRmZtbV1dHp9HXr1m3fvt3Pz8+ULSAg4Pjx4xs3brx37x6Px4uMjHR1dTUbNFsMgUCQy+V+fn5arZZEImk0GpVKBe++myVoGhaKiIiIK1euzHUVAIAFRyaTxcTEBAcHt7a2mv6vY2qwvLw8ODg4OTkZIRQbG1tYWFhaWuri4oLD4dLS0vR6vcFg0Ov1ExMTtra2ycnJxpEIoVWrVnl6eiKEzAansrOzE4lEQqEQj8cHBgZifPbWBuY0AAAAwMr4+DiPx0tLSyspKTF1DNMFbWxsTDsaDAaNRtPQ0CCRSHA4nL29fUtLS3JyckJCglwuz8rKmpiYQAg9fPiwubmZw+GYDZotSSqV5uXl1dbWtrW1mZ72BLMETQMAAACsSCQSnU63fv36wceUSqXZII/Hq6urKy4uViqVUqm0vLw8Nja2oKDA8BiTyTx9+nRpaenixYuLi4tzcnIGBwc/+OADDoezZMkSs0GzJSkUCgcHB0dHR7VavWPHDoSQWq02fvWzK/gAaBoAAABgpb29fWBgwNvb2+ux/fv3mw3SaLSSkpLs7GwPD49t27aJRCLT4xJPIBKJp06d+vvf/758+fJHjx6dOnVquqBZKSkpDAaDQqGEhITQaLS4uLjY2FiEUGJiIofDqa+vx+hSWAcbg8Ew1zXMsaampo8++ujSpUsv8qAxMTFZWVkWXOXSw8Ojo6NjhgWrLE/5dM9Az9L8WrAKF4bJkkIYLViFxbUdwmaZIuL/396dxzV15/vj/4RFkABKI0ugMldk8QFlCSJOfrY3IoNQwlIrghYupQ+HAn3YcYzdBpehTkfGXpF5TDvXite5UKUK+sUQIostClpwqfSKKESKZREXLA2bSQgE8vsjNZeR4KDmiCd5Pf/ogxzeefHmTKe+PTnn8/GkJJZeG1ZZraXVhlUyud5j4TmEGyEBAMAAXb58OScnZ/Lx5OTksLCwZ9+PYcDQAAAABiggIODgwYMz3YWhwT0NAAAAMC240gBPSNVPyW0Cs0LzqIg1m6v7PuqnRNEdGEN/pSR29Ef9Z9pF6D+TOqvzKLn54PuPqUgl7YX6v/mAELLARf8ngWFqqvdMeD7hSgMAAABMC4YGAAAAmBYMDQAAADAtGBoAAIBCV69e5fF4tra2Hh4e2mWbxWKxj48Pk8kMDAysqal56C0jIyNjY2OPqOzq6goPD587dy6Px2ttbdW+KyMjw8nJyd3dvays7Bn8akYIQwMAAFBFLpeHhITweLzOzs7c3Ny0tLSGhoZ79+4lJiZmZmZ2dHTw+fzVq1fL5fKJb/Hz8xMKhYSQqSqTkpLc3d0lEgmXy12zZo3mjRkZGd3d3RcvXty8efMbb7zR29tLxW/E5XLb29upSKYFDA0AAECV+vp6c3PzHTt22NnZRUVFRUZGVlRU1NfXe3l5JSYm2tvbZ2ZmDg0NtbW1ad8iEAi0Fw90Vra2tl64cGHXrl1OTk6ffPJJR0fHpUuXhoaGjhw5snfvXldX14yMjKSkpI6ODip+o7a2tpERSh7DoQUMDQAAQBVXV9c9e/ZoX/b29rLZ7LCwsPLycs0RiURiYmLi7u6ueVlWVnbt2jXtrhM6K1taWnx8fKytrQkhZmZmgYGBzc3N9fX1rq6uL774oqZ47969QUFBU3W1fft2Z2dnGxubqKiovr6+qcpKS0v9/PysrKzc3Nzy8/MJIeHh4VKpNDQ0tK6u7klPCb1haADdhEKhpaXlC/9s7969M90XANCJp6fn2rVrNV8XFBQ0NjZGRkYymcx58+a1t7dzOBwul1tUVGRlZUUI6enp2bRp05dffmn6YOEHnZV3795lsVjaH8FisTRbZdrZ2aWmprJYLHd3d82f8TpVVlbu3bv366+/vn79el9fX25urs6y/v7++Pj4P/zhD3fu3Nm5c2d6evro6GhVVdULL7xQXV29bNkyfZ0iesHiTqCbtbU1l8stKSnRHmEwGHPnzp3BlgCApgYHBwUCgVAoPHnyJJvN1hx0dHTMzs4WCoUCgWDp0qVOTk5vvfXW1q1bFyxY8NDbH6ocHx9/qEClUvX39587dy4uLq6jo+O7776LiIjw8vLicrmTm1GpVGNjY11dXZ6enrW1tVP1zGQyJRLJggULRkZGWCyWUqmUy+Vz5sx5ujNBe7jSAFMyMzOzmwATAwA8gebmZg6HI5VKGxsbg4ODCSGDg4P379+3srKKiIj44osvzM3Nq6qqNF8kJSWNjY2p1eqxsbHx8XGdlY6OjlKpVJvf19fn5ORkZ2fn6ekpEAhsbGxWrFixfPnyEydO6OwnKioqKytr69atdnZ2sbGx2vsnHmJqalpYWOju7h4UFPSI6xbGBkMDAABQZXh4mM/nJyUllZSUuLi4aA5++umn77zzjrbGwsJCJpPV1dWJRCJzc3MzM7OGhoaEhITVq1frrPTz82tublYoFISQsbGx77//3tfX193dfXR0VFtpZ2c3e/ZsnS21t7fz+fyGhoabN2+6ubkJBAKdZWVlZXl5ebW1tVeuXNE+KQoYGgAAgCoikUilUqWnp/c8IJPJ+Hx+aWlpTU2NXC4vLi6+fv36ypUrDx48qH4gKCjo2LFjx48f11mp+dv/zp07R0ZGcnJyXF1dg4KCuFyutbV1dnb20NDQiRMnysvLV61apbOlysrK1157rauri8FgMBiMqR6FkEqlFhYWlpaWCoVi27ZthBDNmEIIkclkFJ2u5x+GBgAAoEpTU1N3d7ezs7PTA7t37+Zyufv27du0aZODg8OuXbuEQqH26YmHTFVZXFx8/vx5NptdUVEhEok0xWKxuLq62sXFJTMz88iRI97e3joz169fHxQUxOFwXF1dJRJJXp7uTfLWrVsXEBDg5ubG4XA8PDxiY2NjYmIIIXFxcSEhIWfPntXD2aEhhlqtnukeZlhVVVV6evozXqwjOjo6LS0tKipKX4EODg7Xrl2zt7fXV+A333yza9eur7/+eqoC1a39+vpZzwC9drkco2QDUUp2uTR3038mdUwp+beAsl0ub82iIpaiXS6Xtar0HgvPITw9QZhMpuZ5XwAAMBiXL1/OycmZfDw5OTksLOxxy0ADQwMxNzdnMpkz3cXMu3//fklJibOzs+bld999N/H+ZAAAegkICDh48KC+ykADQwP8QiQSCQQCDoejeXnnzp2ffvrp2bdB0ecIqn5qrvjTCr0+SqACRZ/7UPQ5AkWo6JZhamqkSx0ZHwwN8AsHBwcOh6O9iUFzT8PMtgQAAM8VPD0BAAAA04KhAQAAAKYFQwMAAFDo6tWrPB7P1tbWw8NDu7Ti9A+eO3cuODiYyWRyOJxTp05pDorFYh8fHyaTGRgYWFNT84hK0C8MDQAAQBW5XB4SEsLj8To7O3Nzc9PS0hoaGqZ/UKlUxsbGhoeHt7a2JiQkrFq1amBg4N69e4mJiZmZmR0dHXw+f/Xq1XK5XGflY7XK5XKf8YI9dIQbIWfG8PBwW1tbQ0ODvgKnWgkVAGAG1dfXm5ub79ixgxASFRUVGRlZUVHR19c3zYM2Nja9vb0ffPCBjY3N5s2bs7KyWltbb9265eXllZiYSAjJzMzctWtXW1ubpaXl5MolS5ZMv9W2tjb8h/RfwtAwM8zMzPbv33/o0CF9BQ4NDWnXRQcAeE64urru2bNH+7K3t5fNZk//oJubG5vN3r9/f0ZGxqFDh2xsbLy9vb29vV9++WVNmUQiMTExcXd3nzVr1uTKqboqLS3dtm1bW1ubk5PT9u3bU1JSwsPDpVJpaGhoUVHRsmV4gHRKGBpmRkVFhX4DHRwcptrS7clIpdKWlpa0tLSJB9esWfOb3/xGjz8FAAybp6enp6en5uuCgoLGxsbIyEg2mz3Ng2ZmZocOHVqxYsX7778/Pj4uEok0a/Exmcz29vbXX3+9paWlqKjIysqKEKKzcrL+/v74+Pj8/HzNxYyUlJTExMSqqip7e/vq6movLy/KTwqdYWgA3V555ZWkpCQ3t39aD0i7sy0AwPQNDg4KBAKhUHjy5Ek2mz39g93d3YmJiUeOHOHz+adOnUpNTT1z5oxmtnB0dMzOzhYKhQKBYOnSpSqVaqrKhzCZTIlEsmDBgpGRERaLpVQq5XL5nDlzntnZoDUMDaAbm83+y1/+MtNdAADtNTc3R0dH+/v7NzY2av/iMc2DYrHY398/ISGBEBITE3P48OHjx49nZGSYmJhYW1tHREREREQsWrSoqqpKoVBMrvzwww8n92NqalpYWJifn29lZeXr6/tsToLBwNMTAABAleHhYT6fn5SUVFJSoh0OHusgg8HQpqnVaqVS+emnn77zzjvagxYWFjKZTGelzpbKysry8vJqa2uvXLmifbATpglXGoxXc3OzRCKxtbXVvMQOVQCgdyKRSKVSpaen9/T0aI5YW1ufOHFimgf5fP7WrVuPHj0aGRl56tQpsVj80UcfKRSKiIiImpqa4OBgsVh8/fr1lStXqtXqyZU6W5JKpRYWFpaWlgqF4k9/+hMhRKFQaD6ekMlklJ8RmsPQYLyysrL+93//99/+7d80L+/cuXP37t0Z7QgADE1TU1N3d7d2+1xCSFZWlkqlmubBP/7xjyUlJR999FFKSsrChQsLCwsDAgIIIfv27du0adMPP/zg5eUlFArd3d0JITorJ1u3bl15ebnmuYwPP/wwNjY2Jibm4sWLcXFxISEhYrH4lVdeoep00B9DrVbPdA8z7MKFCxs3bjx//vxMN/JUHBwcrl27Zm9vP/23bNmyhclkZmZmal5qdqjSblj1L6lu7X/sLqeBXrtcjlS/TUUsRZsxAkWEf6LTLpdUYJiaJsnkM90FPAu40gAAAAbo8uXLOTk5k48nJyeHhYU9+34MA4YGAAAwQAEBAQcPHpzpLgwNnp4AAACAacGVBni+UHTzAVW3SrBfpiKWkG+pidW/4dqZ7uA5sMCFkg0L2m9RcqsEFd0yTE31ngnPJ1xpAAAAgGnB0AAAAADTgqEBAAAodPXqVR6PZ2tr6+Hh8dAKjLdv3/7xxx+1L8+cORMcHGxra8vj8a5fv/6Iyvj4eMYEt27devTbQV8wNAAAAFXkcnlISAiPx+vs7MzNzU1LS2toaNB+d/v27cePH9d8ffPmzejoaIFA0N3dHRcXFxYWplKpdFYSQm7cuFFYWNjxgJOT06PfrkdcLre9vZ2KZFrA0AAAAFSpr683NzffsWOHnZ1dVFSUZjdqQsjhw4ejoqIOHDigrRSLxcHBwWvXrrW1tX333Xft7OwqKyt1VhJCbty48corr/zqAVNT06nerndtbW0jI5Tc+koLGBoAAIAqrq6ue/bs0b7s7e3V7ILNYrFWrVq1ePFi7beGh4f7+vq0L2fPnv3DDz/orPz5558HBwdTU1Otra19fX3Ly8sf8Xadtm/f7uzsbGNjExUVNfFdDyktLfXz87OysnJzc8vPzyeEhIeHS6XS0NDQurq6xzsRhgKPXMIvVCqVVCr95ptvpipgMBgBAQEsFutZdgUAtObp6enp6an5uqCgoLGxMTIykhCycuVKQsjEiwErV6784IMPysrKVqxYUVxcfOHChdDQUJ2VN27cMDc3T05O/uqrr44dO/b6669fu3ZtqrdPVllZuXfv3pqaGjs7uzVr1uTm5u7YsWNyWX9/f3x8fH5+vubqSEpKSmJiYlVVlb29fXV1tZeXl97OEa1gaIBfWFpaDg4O7tq16xE1qamp8fHxz6wlADAMg4ODAoFAKBSePHlSc6VhMh8fn0OHDm3atOnu3bt8Pj8sLGzu3Lk6K4ODg7XbXr/99tuHDx8+fvz4e++9N823q1SqsbGxrq4uT0/P2topVxphMpkSiWTBggUjIyMsFkupVMrlcs1mmMYMQ4OxGBkZGRwcNJ2wBotcLre0tNS+XL58+SOu5gEAPJnm5ubo6Gh/f//GxkYXF5dHVCYkJCQkJGi+fuWVVxwdHaeT7+XldefOnem/PSoqKisra+vWrWvWrOHxeP/5n//p7e09uczU1LSwsDA/P9/KysrX13c6nRgD3NNgLN58801XV9eFE/zXf/3XyZMnZ7ovADBkw8PDfD4/KSmppKTk0RODRCJJS0sbHx8nhAwNDV26dCkkJERn5ZdffpmWlqZ9+eOPPy5cuHD6b29vb+fz+Q0NDTdv3nRzcxMIBDrLysrK8vLyamtrr1y58tCTosYMQ4Ox8PLy+uCDD6QTvPfee6+++upM9wUAhkwkEqlUqvT09J4HZDKZzsoXX3zx6NGjOTk5PT09GzZsCAkJmT9f9+rvL7300v79+w8cONDf319QUNDQ0JCQkDD9t1dWVr722mtdXV2aNR6mehRCKpVaWFhYWloqFIpt27YRQhQKheZbU/0KxgBDAwAAUKWpqam7u9vZ2dnpgd27d+ustLa2Lioq+u///m9vb+/R0dGioqKpMgMDA4VC4eeffz5//vx9+/ZVVVWxWKzpv339+vVBQUEcDsfV1VUikeTl5eksW7duXUBAgJubG4fD8fDwiI2NjYmJIYTExcWFhIScPXv2Mc+EgWCo1eqZ7mGGXbhwYePGjefPn5/pRp6Kg4PDtWvX7O3tpyrIysrS/lNjy5YtTCYzMzPzyX6i6tb+J3vjjKBow6rh+j9TETt+BxtW0UlzDSWx9NqwalkrJSspwfMGN0ICAIABunz5ck5OzuTjycnJYWFhj1sGGhgaAADAAAUEBBw8eFBfZaCBexoAAABgWnClAZ4vY72tM93CYxi9SMnNB+bBL1MRSwVzN0rOgCklt6AQEzYlJ/bkP3qoiKUIz5mC5Vjw10+jgf+pAQAAYFowNAAAAMC0YGgAAACAacHQAAAAFLp69SqPx7O1tfXw8NCux9zV1RUeHj537lwej9fa+sudTOfOnQsODmYymRwO59SpU5qDZ86cCQ4OtrW15fF4169ff0SlWCz28fFhMpmBgYE1NTXP9Jc0GhgaAACAKnK5PCQkhMfjdXZ25ubmpqWlNTQ0EEKSkpLc3d0lEgmXy12zZg0hRKlUxsbGhoeHt7a2JiQkrFq1amBg4ObNm9HR0QKBoLu7Oy4uLiwsTKVS6ay8d+9eYmJiZmZmR0cHn89fvXq1XC6n4jficrnt7e1UJNMChgYAAKBKfX29ubn5jh077OzsoqKiIiMjKyoqWltbL1y4sGvXLicnp08++aSjo+PSpUudnZ29vb0ffPCBi4vL5s2blUpla2urWCwODg5eu3atra3tu+++a2dnV1lZqbOyvr7ey8srMTHR3t4+MzNzaGiora2Nit+ora1tqu0qjAGGBgAAoIqrq+uePXu0L3t7e9lsdktLi4+Pj7W1NSHEzMwsMDCwubnZzc2NzWbv379foVDk5+fb2Nh4e3sPDw/39fVp3z579uwffvhBZ2VYWFh5ebmmTCKRmJiYuLu7T9XV9u3bnZ2dbWxsoqKiJuY/pLS01M/Pz8rKys3NLT8/nxASHh4ulUpDQ0Pr6uqe+tzQEoYG49XS0vLJJ5+88Dim2tkFAEAnT0/PtWvXar4uKChobGyMjIy8e/cui8XS1rBYrJ6eHjMzs0OHDm3evNna2vrtt9/+xz/+wWQyV65c2djYWFZWJpPJ/ud//ufChQu9vb06K5lM5rx589rb2zkcDpfLLSoqsrKy0tlSZWXl3r17v/766+vXr/f19eXm5uos6+/vj4+P/8Mf/nDnzp2dO3emp6ePjo5WVVW98MIL1dXVy5Yt0/u5ogUs7mS8vvrqq/v375uamk7/LTY2NtT1AwCGanBwUCAQCIXCkydPstns8fHxhwpUKlV3d3diYuKRI0f4fP6pU6dSU1PPnDnj4+Nz6NChTZs23b17l8/nh4WFzZ07V2elp6cnIcTR0TE7O1soFAoEgqVLlzo5OU1uRqVSjY2NdXV1eXp61tZOueUak8mUSCQLFiwYGRlhsVhKpVIul8+ZM0e/Z4Z2cKXBeFlaWs6bN8/ucZiZYcoEgMfT3NzM4XCkUmljY2NwcDAhxNHRUSqVagv6+vqcnJzEYrG/v39CQoK1tXVMTExISMjx48cJIQkJCW1tbffv3y8qKlIoFI6OjjorBwcH79+/b2VlFRER8cUXX5ibm1dVVensJyoqKisra+vWrXZ2drGxsdpnNx5iampaWFjo7u4eFBSk+WwCCIYGAACgzvDwMJ/PT0pKKikpcXFx0Rz08/Nrbm5WKBSEkLGxse+//97X13d4eJjBYGjfqFarlUqlRCJJS0vTXJkYGhq6dOlSSEiIzspPP/30nXfe0R60sLCQyWQ6W2pvb+fz+Q0NDTdv3nRzcxMIBDrLysrK8vLyamtrr1y5on1SFDA0AAAAVUQikUqlSk9P73lAJpNp/vq+c+fOkZGRnJwcV1fXoKAgPp9/5syZo0ePymSysrIysVgcExPz4osvHj16NCcnp6enZ8OGDSEhIfPnz9dZyefzS0tLa2pq5HJ5cXHx9evXV65cqbOlysrK1157rauri8FgMBiMqR6FkEqlFhYWlpaWCoVi27ZthBDNlEMImWocMQa42mwghoeHP/vsMyaTOVVBVVWVq6vrs2wJAKCpqam7u9vZ2Vl7JCsr649//GNxcXFycjKbzfbz8xOJRIQQDw+PkpKSjz76KCUlZeHChYWFhQEBAYSQoqKiDRs2/OUvfwkPDy8qKnpE5b59+zZt2vTDDz94eXkJhcKpnp5Yv379xYsXORzO6Ojor3/966nu7163bl15ebnmSY0PP/wwNjY2Jibm4sWLcXFxISEhYrH4lVde0fvpev4x1Gr1TPcwwy5cuLBx48bz58/PdCNP5c3rkdq0AAAbl0lEQVQ333Rycpp4ye4hJ0+eXLhw4dGjR/X1E1W39usraiKKdrk0nedJRawi/20qYmm0yyVF+3zSa5fLnSvptMvlb1+mYpdLU9dalf5j4fmDKw0GoqCg4NEFlpaWz6YTAIDnweXLl3NyciYfT05ODgsLe9wy0MDQAAAABiggIODgwYP6KgMNDA3whMzmUnIFmaJYilB0FZ2ia/5UGP2RTrGEUHJif0vNp0ldt2kTyzAluGHKSODpCQAAAJgWDA0AAAAwLRgaAAAAYFowNAAAAIWuXr3K4/FsbW09PDweWlrx9u3bP/74fzewiMViHx8fJpMZGBhYU1PzUM7IyMjY2Bgh5O9//zvjn61atYoQcu7cueDgYCaTyeFwTp06RfkvZpQwNAAAAFXkcnlISAiPx+vs7MzNzU1LS2toaNB+d/v27ZoNJggh9+7dS0xMzMzM7Ojo4PP5q1evlsvlE3P8/PyEQiEhJDk5uWMCLpcbHx+vVCpjY2PDw8NbW1sTEhJWrVo1MDDwWK1yudz29nZ9/NKGDEMDAABQpb6+3tzcfMeOHXZ2dlFRUZGRkRUVFYSQw4cPR0VFHThwYGKll5dXYmKivb19Zmbm0NBQW1ub9rsCgUC7s5SNjc2vHmhqamIymevWrevs7Ozt7f3ggw9cXFw2b96sVCqn2olqKm1tbVMtKQ1aGBoAAIAqrq6ue/bs0b7s7e1ls9mEEBaLtWrVqsWLF2u/FRYWVl5ervlaIpGYmJho14EuKyu7du2aZq3oiVQq1ebNm//6178SQjTrPe/fv1+hUOTn59vY2Hh7e0/VVWlpqZ+fn5WVlZubm2YHy/DwcKlUGhoaWldXp5df3FBhnQbQTaVSxcbG+vr6TjzI5/ONc7l1AHgynp6enp6/LOJeUFDQ2NgYGRlJCNHsJlVZWamtZDKZTCazvb399ddfb2lpKSoqsrKyIoT09PRs2rTp66+/jo+Pfyi8sLDwpZde8vHxIYSYmZkdOnRoxYoV77///vj4uEgkmmovnv7+/vj4+Pz8fM1lj5SUlMTExKqqKnt7++rqai8vLwpOg+HAlQbQbWho6PTp03aTzHRfAEA/g4ODv/3tbzdv3nzy5EnNlYapODo6Zmdnp6SkCASCu3fvEkLeeuutrVu3Lliw4KFKtVqdnZ393nvvaV52d3cnJiYeOXJkYGCgtLQ0NTV1qo8nmEymRCJZt27d7NmzWSyWUqmcePMEPBquNIBuDAbDwsLiww8/nOlGAIDempubo6Oj/f39GxsbXVxcpiobHBw0MTGxtraOiIiIiIhYtGhRVVWVQqEwNzdPSkoaGxtTq9VjY2Pj4+MmJiaEkPr6+pGRES6Xq3m7WCz29/dPSEgghMTExBw+fPj48eM6/wtmampaWFiYn59vZWX10MVU+JdwpQEAAKgyPDzM5/OTkpJKSkoeMTEQQj799NN33nlH+9LCwkImk9XV1YlEInNzczMzs4aGhoSEhNWrV2sKvvrqqzVr1kz8QRO3+VWr1UqlUucPKisry8vLq62tvXLlykOPgMK/hKEBAACoIhKJVCpVenp6zwMymUxnJZ/PLy0trampkcvlxcXF169fX7ly5cGDB9UPBAUFHTt2TPuIZlVVFY/Hm/j2M2fOHD16VCaTlZWVicXimJgYnT9IKpVaWFhYWloqFIpt27YRQhQKheZbU/UGWhgaAACAKk1NTd3d3c7Ozk4P7N69W2cll8vdt2/fpk2bHBwcdu3aJRQKtU9PTNbd3X3jxo3g4GDtEQ8Pj5KSkuzsbAcHhy1bthQWFk5+2kJj3bp1AQEBbm5uHA7Hw8MjNjZWM17ExcWFhIScPXv26X5jA8dQq9Uz3cMMu3DhwsaNG8+fPz/TjVArKytL+8/p6O/vX7BgQV9f35QVssopv2U05P/vVSpix25SkUoJyrajpJP7j7ccwHRRtMslFRimpstaVTPdBTwLuBESAAAM0OXLl3NyciYfT05ODgsLe/b9GAYMDQAAYIACAgIOHjw4010YGgwNxkIikXz99dfT/7/Q+Pj44OAgpS0BAAC9YGgwFgcOHLh9+7apqek06wcHB5cvX05lR7qp+in5PN/MxY+KWBP2y1TEjt38lopYGqHoLgEAeEoYGowFk8n08PCYfn1/f//Eh54BAADwyCUAAABMC4YGAACg0NWrV3k8nq2trYeHh3YFxq6urvDw8Llz5/J4PO0mEfHx8YwJbt26RQg5d+5ccHAwk8nkcDinTp16RKbOg6BfGBoAAIAqcrk8JCSEx+N1dnbm5uampaU1NDQQQpKSktzd3SUSCZfL1a4GfePGjcLCwo4HnJyclEplbGxseHh4a2trQkLCqlWrBgYGdGZO9YP0jsvltre3U5FMC1jcyVgWd3pcM7W4E71uhBz++j0qYkcv0uZGSIoWd8KNkMRQFnf65ptvkpOTb9/+5ZdZvXo1h8OJj4/39fX9+eefra2tVSoVi8Wqrq4OCgqaO3duU1PT/PnztW9vbW1dtGjRwMCAjY3N6OiojY3N2bNnBwYGJmf++te/nnxw69atev9l7e3tv/32W6PdQRtXGgAAgCqurq579uzRvuzt7WWz2S0tLT4+PtbW1oQQMzOzwMDA5ubmn3/+eXBwMDU11dra2tfXt7y8nBDi5ubGZrP379+vUCjy8/NtbGy8vb11Zuo8OFVX27dvd3Z2trGxiYqKesRfjUpLS/38/KysrNzc3PLz8wkh4eHhUqk0NDS0rq7uqc4LbeHpCZjS2NjYQ9f3Fi5cOHfu3JnqBwBox9PT09PTU/N1QUFBY2NjZGSkSCRisVjaGhaL1dPTc+PGDXNz8+Tk5K+++urYsWOvv/76tWvXFi5ceOjQoRUrVrz//vvj4+MikYjJZOrMZLPZkw/qbKmysnLv3r01NTV2dnZr1qzJzc3dsWPH5LL+/v74+Pj8/PzIyMiKioqUlJTExMSqqip7e/vq6mqjvdKAoQF0s7CwmD17dlpa2sSDGzdu/I//+I+ZagkAaGpwcFAgEAiFwpMnT7LZ7PHx8YcKVCpVcHCwdjPrt99++/Dhw8ePH1+7dm1iYuKRI0f4fP6pU6dSU1PPnDmjGQ4eytT5g3Q2o1KpxsbGurq6PD09a2trp+qZyWRKJJIFCxaMjIywWCylUimXy+fMmaOH00Fn+HgCdJs9e3ZPT8+lf4aJAQAeV3NzM4fDkUqljY2Nmn0pHR0dpVKptqCvr8/Jyemhd3l5ed25c0csFvv7+yckJFhbW8fExISEhGi2xp6cOdXByaKiorKysrZu3WpnZxcbG6t9duMhpqamhYWF7u7uQUFBms8mgGBoAAAA6gwPD/P5/KSkpJKSEhcXF81BPz+/5uZmhUJBCBkbG/v+++99fX2//PLLiZc2f/zxx4ULFw4PD09cZU6tViuVSp2ZOg/q1N7ezufzGxoabt686ebmJhAIdJaVlZXl5eXV1tZeuXIFD3BqYWgAAACqiEQilUqVnp7e84BMJtP89X3nzp0jIyM5OTmurq5BQUEvvfTS/v37Dxw40N/fX1BQ0NDQkJCQwOfzz5w5c/ToUZlMVlZWJhaLY2JidGbqPKizpcrKytdee62rq0uzGsTIyIjOMqlUamFhYWlpqVAotm3bRgjRTDmEkKmSjQGGBgAAoEpTU1N3d7ezs7PTA7t37yaEFBcXnz9/ns1mV1RUiEQiQkhgYKBQKPz888/nz5+/b9++qqoqFovl4eFRUlKSnZ3t4OCwZcuWwsLCgIAAnZlT/aDJ1q9fHxQUxOFwXF1dJRJJXl6ezrJ169YFBAS4ublxOBwPD4/Y2NiYmBhCSFxcXEhIyNmzZyk7Z881rNOAdRqeFNZpwDoNWKeBSoaxTgMYGDw9AQAABujy5cs5OTmTjycnJ4eFhT1uGWhgaAAAAAMUEBBw8OBBfZWBBu5pAAAAgGnBlQZ4QtLfvTrTLcy8F/5WQUWs5f9HRSoluiLwrwGx9qQk1puaWErgr59GA/9TAwAAwLRgaAAAAIBpwdAAAAAA04KhAQAAKHT16lUej2dra+vh4aFdj7mrqys8PHzu3Lk8Hk+7+8O5c+eCg4OZTCaHwzl16pTm4MjISEZGhpOTk7u7e1lZ2SMqxWKxj48Pk8kMDAysqal5pr+k0cDQAAAAVJHL5SEhITwer7OzMzc3Ny0traGhgRCSlJTk7u4ukUi4XO6aNWsIIUqlMjY2Njw8vLW1NSEhYdWqVQMDA4SQjIyM7u7uixcvbt68+Y033ujt7dVZee/evcTExMzMzI6ODj6fv3r1arlc/litcrnc9vZ2Kk6CIcGKkFgR8glJ1zP+dZGho+jpCRrB0xOEsqcn6MTE9IX9uleE/Oabb5KTk2/f/mV5y9WrV3M4nPj4eF9f359//tna2lqlUrFYrOrqaltb20WLFg0MDNjY2IyOjtrY2Jw9e3bRokVOTk7Xr19/8cUXCSEZGRnr16/XWXnr1q2dO3devHiREKJQKObMmXPp0iU/v8dYCtbe3v7bb7/18vJ66tNhyHClAQAAqOLq6rpnzx7ty97eXjab3dLS4uPjY21tTQgxMzMLDAxsbm52c3Njs9n79+9XKBT5+fk2Njbe3t719fWurq6aiYEQsnfv3qCgIJ2VYWFh5eXlmjKJRGJiYuLu7j5VV6WlpX5+flZWVm5ubpptr8PDw6VSaWhoaF1dHXVnwwBgnQbS19d35cqVhQsXznQjzzsTE5MtW7akpKTMdCMAQBuenp6enr9ciikoKGhsbIyMjBSJRCwWS1vDYrF6enrMzMwOHTq0YsWK999/f3x8XCQSMZnMnp4eOzu71NTUkpISOzu7rVu3pqSk6KwkhDCZzPb29tdff72lpaWoqMjKykpnS/39/fHx8fn5+ZGRkRUVFSkpKYmJiVVVVfb29tXV1bjS8GgYGkhERMT58+c1My88mpOT00y3AAD0Mzg4KBAIhELhyZMn2Wz2+Pj4QwUqlaq7uzsxMfHIkSN8Pv/UqVOpqalnzpzp7+8/d+5cXFxcR0fHd999FxER4eXlNX/+/MmVmtHE0dExOztbKBQKBIKlS5fq/E8Wk8mUSCQLFiwYGRlhsVhKpVIul8+ZM+dZnAj6w9BACCGP9bkXAABMX3Nzc3R0tL+/f2Njo4uLCyHE0dFRKpVqC/r6+pycnMRisb+/f0JCAiEkJibm8OHDx48fd3Z29vT0FAgEhJAVK1YsX778xIkTL7744uTKjIwMExMTa2vriIiIiIiIRYsWVVVVvfnmm5P7MTU1LSwszM/Pt7Ky8vX1fSbnwHDgngYAAKDK8PAwn89PSkoqKSnRTAyEED8/v+bmZoVCQQgZGxv7/vvvfX19h4eHGYz/u71arVYrlUp3d/fR0VHtQTs7u9mzZ+us/PTTT9955x3tQQsLC5lMprOlsrKyvLy82traK1euaB8BhWnC0AAAAFQRiUQqlSo9Pb3nAZlM5u7uHhQUtHPnzpGRkZycHFdX16CgID6ff+bMmaNHj8pksrKyMrFYHBMTw+Vyra2ts7Ozh4aGTpw4UV5evmrVKp2VfD6/tLS0pqZGLpcXFxdfv3595cqVOluSSqUWFhaWlpYKhWLbtm2EEM34QgiZas4ALQwNAABAlaampu7ubmdnZ6cHdu/eTQgpLi4+f/48m82uqKgQiUSEEA8Pj5KSkuzsbAcHhy1bthQWFgYEBBBCxGJxdXW1i4tLZmbmkSNHvL29dVZyudx9+/Zt2rTJwcFh165dQqFwqqcn1q1bFxAQ4ObmxuFwPDw8YmNjY2JiCCFxcXEhISFnz559hqeHfrBOAzwhrNNAsE4D1mkghGCdBvKodRrAwOBGSAAAMECXL1/OycmZfDw5OTksLOzZ92MYMDQAAIABCggIOHjw4Ex3YWjw8QQAAABMC26EBAAAgGnB0AAAAADTgqEBAAAApgVDAwAAAEwLhgYAAACYFgwNAAAAMC0YGgAAAGBasLgT6N+BAwceXbB+/XrEGnwsjVpFLKWxYEiwuBPoH4PBYDAYc+bMmaqgr68PsQYfS6NWEUtpLBgUNYC+vfrqq+bm5i+99NLHH3/c3NyMWOOMpVGriKU0FgwJhgaghFQqzc/P5/P5s2bN0uN/gxBLr1gatYpYSmPBYGBoAGr19/d/+eWX0dHRFhYWPj4+WVlZ165dQ6yxxdKoVcRSGgt0h6EBnpGBgYGCggJXV1f9fiiGWHrF0qhVxFIaCzSFpyeAcmq1uq6urri4+NixYzKZLCkpCbFGGEujVhFLaSzQ28zMKmAExsfH6+vrf//737u4uNjY2LzxxhtCoXB4eBixRhVLo1YRS2ksGAYMDaB/Fy5c2Lx58/z585lM5tq1a0tKShQKBWKNLZZGrSKW0lgwJFinAfRP87T3smXLIiMjraysJhds3LgRsQYfS6NWEUtpLBgSDA2gf/Pnz390wc2bNxFr8LE0ahWxlMaCIcHQAAAAANOCpydA/+i1MD5iKYqlUauIpTQWDAmuNID+0WthfMRSFEujVhFLaSwYlJm9DxMMEr0WxkcsRbE0ahWxlMaCIcHQAJSg18L4iKUolkatIpbSWDAYGBqAWvRaGB+xFMXSqFXEUhoLdIehAZ4Rei2Mj1iKYmnUKmIpjQWawtMTQDk1rRbGRyxFsTRqFbGUxgK9zcysAkaAXgvjI5aiWBq1ilhKY8EwYGgA/aPXwviIpSiWRq0iltJYMCRYpwH0j14L4yOWolgatYpYSmPBkGBoAP2j18L4iKUolkatIpbSWDAkGBoAAABgWkxmugEAAACgBwwNAAAAMC0YGgAAAGBaMDQAAADAtGBFSKDQkiVL6urqZs2apT1y7969+Pj4mpoaxBp8bHV19aMLQkNDHzcTsXSMBUOCoQH0r6+v729/+xsh5NKlS3/605/MzP7vX7MbN260tLQg1hhiIyIiHl0wOjqKWGOIBUOCoQH0T6FQfPvtt5qv6+rqTE1Ntd8yMTHJzc1FrDHEUvQHDGJpFwuGBOs0AIUmX+tGrHHGqtXqwcFB7cvW1tbf/OY3AwMDiDW2WKA7XGkACn333Xc3btxoa2t76Hh4eDhijSf29OnT8fHxvb29Ew/yeLwnDkQsTWPBEMzUphdgDD7//HMGg2Fubm79zxBrVLFLliyJiYk5ffo0m80uLy+vqanx9vZubW19ylYRS7tYMAAYGoBCLi4u6enpg4ODiDXmWCsrq5qaGrVanZiYWFpaqlarhUJhdHQ0Yo0tFgwA1mkACkml0g0bNtjY2CDWmGNtbW3lcjkhxN3dXfMghq+vr/a+S8QaTywYAAwNQKFFixY1NDQg1shjX3755Y8//rilpcXf37+4uPinn34SiUQ6d15GrGHHggHAjZBAoe3bt2dkZEgkkqCgIAsLC+1xPp+PWOOJzcnJiY2NLSsr+/3vf//nP/+ZzWaPjY3t3r37afpELB1jwQDgkUug0FRXuYeGhhBrVLFacrm8trZ23rx5S5Ys0UsgYmkaC3Q10zdVAICBCwoKUiqVE4/09PTweDzEGlssGAB8PAHUUtNq5RnE6jGWRiteI5bSWDAkGBqAQvRaeQax+o2l0YrXiKU0FgzKTF7mAENHr5VnEEtR7ORr3XqBWNrFggHAjZBAISaTWV5ezuPxkpKS4uPjY2JiSktLDxw4IBKJEGtUsVSseI1YOsYC3eHjCaDQQ0vExMTEULHyDGKf89i///3v7777rpmZ2cRnOMlTP5GBWNrFgiGY6UsdYMji4uKWLl3a3NxcUlISGBh479693NxcFxcXxBpVLI1WvEYspbFgADA0AIU6OzsDAgJ27dqlVCoXL16sua9q9+7diDWq2NmzZ1+9evUpG0OsAcSCAcDQAM+ITCYrLy+/ePEiYo0tlsPhFBQU6KUlxNI6FgwAboQECi1ZsqSurm7WrFnaI/fu3YuPj6+pqUGs8cQKhcKMjIy33npLvyteI5Z2sWAAMDSA/mmXiMnKytq6detDS8RUVVX19PQg1hhiNei14jViqYsFA4CnJ0D/6LXyDGKpi9Wg6E8axNIuFgzBTH8+AoaMXivPIJa62I6Ojg8//DApKWlkZOTEiRNSqRSxxhkLdIehASg0Pj5+5MiRa9euqdXqzz77jMPhpKWlyWQyxBpV7NWrV52dnRcvXkwIGR4eXr58OZvNbmpqespWEUu7WDAAGBqAQpmZmQwGo7KysrW11cTE5He/+523t/fvfvc7xBpV7IoVKzZs2KBWqzV/Ag0ODsbExISHhz9lq4ilXSwYAAwNQCFnZ+cvvvhCrVbn5OS8/PLLarX62LFj8+fPR6xRxTKZTM1zm5o/gdRq9enTp21tbZ+yVcTSLhYMgMmzvocCjElfX5+fnx8hpLa2duXKlYSQhQsX/vTTT4g1qth58+b19/dPPPLzzz+zWKynyUQsHWPBAGBoAAp5e3t/9dVXTU1NlZWVmj+BKioqfvWrXyHWqGJjY2O3b9+ufWjzu+++27hx42uvvfaUrSKWdrFgCGb6UgcYsoqKCs3KMP/+7/8+Pj7+8ccfMxiMvXv3ItaoYu/fvx8ZGWliYkIIsbW1ZTAYb7zxhuaiN2KNKhYMABZ3Amrdvn37hx9+WLp0qaWlZW1trZmZ2bJlyxBrhLGXL1++du2ara3tSy+9tGDBgqcPRCxNY4HWMDSA/tXV1bHZbDc3t7q6Op0FT/aHEGJpF0sIkcvloaGhhw4dWrhw4ZMlINYwYsFAzPSlDjBAhBDt81p6/LcOsbSL1YiOjs7Ly3uaBMQaRiwYAFxpAP0bGxtjMBiaD0QRa8yxGufOnRMIBCEhIYsXL7a2ttYeDw8PR6xRxYIhmOmpBQwZjVYtRCx1sZZTeMpWEUu7WDAAGBqAQjRatRCxlMYCgGHA0AAUotGqhYilNHZ8fLx/gosXL+pleUHE0i4W6A5bYwOFaLRqIWKpiz19+nR8fHxvb+/Egzwe72kyEUvHWDAEMz21gCFbvHjxhg0brly5MmvWrPPnz6vV6p07d3p5eSHWqGKXLFkSExNz+vRpNptdXl5eU1Pj7e3d2tr6lK0ilnaxYAAwNACFaLRqIWKpi7WysqqpqVGr1YmJiaWlpWq1WigURkdHP2WriKVdLBgAfDwBFIqIiPjxxx81ywsyGAwej3f27NmnX14QsfSKtbW1lcvlhBB3d/eWlpaYmBhfX99vv/32KVtFLO1iwQBgwyqg1ujoaEVFRWpq6ujoqEwm8/b2Rqyxxb788ssff/xxS0uLv79/cXHxTz/9JBKJrKysEGtssWAIZvpSBxiyq1evOjs7L168mBAyPDy8fPlyNpvd1NSEWKOK7ezsDAgI2LVrl1KpXLx4sampKSFk9+7dT9kqYmkXCwYAK0IChUJDQ729vT/77DMGgzE8PDwyMpKUlKRUKisrKxFrVLFacrm8trZ23rx5S5Ys0UsgYmkaCzT1/wMs8nncRZfsfAAAAABJRU5ErkJggg==" />

<!-- rnb-plot-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->








####### clustering de genes e/ou amostras #######

# Clustering e outros métodos não supervisionados

## Redução de dimensionalidade

### PCA

...


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuXG5jZWxsX2xpbmVfdHlwZSA9IHJiaW5kKGRhdGEyMCRjZWxsLmxpbmUsIGRhdGEyMCRjZWxsLnR5cGUpXG5cbnBjYSA9IHByY29tcChjZWxsX2xpbmVfdHlwZSlcbm1pbih3aGljaChzdW1tYXJ5KHBjYSkkaW1wb3J0YW5jZVszLF0+MC44KSlcblxuY2VsbF9saW5lX3R5cGUkbGluZWFnZSA9IDBcbmNlbGxfbGluZV90eXBlJGxpbmVhZ2VbMTpsZW5ndGgoZGF0YTIwJGNlbGwubGluZVtdKV0gPSBcInNvbGlkb1wiXG5jZWxsX2xpbmVfdHlwZSRsaW5lYWdlW2xlbmd0aChkYXRhMjAkY2VsbC5saW5lW10pKzE6bGVuZ3RoKGRhdGEyMCRjZWxsLnR5cGVbXSldID0gXCJsaXF1aWRvXCJcbmNlbGxfbGluZV90eXBlJGxpbmVhZ2UgPSBhcy5mYWN0b3IoY2VsbF9saW5lX3R5cGUkbGluZWFnZSlcblxuIyBUZW50YXIgY3JpYXIgZ3J1cG9zIGRlIGNvcmVzIHBlbGEgY2xhc3NpZmljYcOnw6NvIGxpcXVpZG8gLyBzb2xpZG9cblxucGFpcnMocGNhJHhbXSwgcGNoPTE4LCBtYWluID0gXCJHcsOhZmljb3MgZGUgZGlzcGVyc8OjbyBtw7psdGlwbG9zXCIsIGNvbD1jZWxsX2xpbmVfdHlwZSRsaW5lYWdlKVxuXG5Jc3RvIMOpIHVtYSBleHBlcmnDqm5jaWFcbmBgYCJ9 -->

```r

cell_line_type = rbind(data20$cell.line, data20$cell.type)

pca = prcomp(cell_line_type)
min(which(summary(pca)$importance[3,]>0.8))

cell_line_type$lineage = 0
cell_line_type$lineage[1:length(data20$cell.line[])] = "solido"
cell_line_type$lineage[length(data20$cell.line[])+1:length(data20$cell.type[])] = "liquido"
cell_line_type$lineage = as.factor(cell_line_type$lineage)

# Tentar criar grupos de cores pela classificação liquido / solido

pairs(pca$x[], pch=18, main = "Gráficos de dispersão múltiplos", col=cell_line_type$lineage)

Isto é uma experiência
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


####### Fim de clustering de genes e/ou amostras #######


<!-- rnb-text-end -->

