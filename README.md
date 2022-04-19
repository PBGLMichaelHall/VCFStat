# VCFStat is a package I made to analyze info fields in a VCF files

---
title: "VCFstat"
author: "Michael Hall"
date: "4/14/2022"
output: pdf_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, message = FALSE, comment = FALSE,warning=FALSE,results = FALSE)
```

```{r Best}
  

devtools::install_github("PBGLMichaelHall/VCFStat",force = TRUE)

library(VCFstat)
library(vcfR)
library(data.table)
library(QTLseqr)
library(ggplot2)

```


```{r Plots}
setwd("/home/michael/Desktop/QTLseqr/extdata/")
chromlist <- c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10")


VCFstat::ChromQual(vcf = "freebayes_D2.filtered.vcf.gz", chromlist = chromlist,windowSize = 1e+05, binwidth = 50,Maximum=5000)
![1](https://user-images.githubusercontent.com/93121277/163991880-7eff336b-6eb2-45cf-ada4-0ae3787b489f.png)


VCFstat::ChromDP(vcf = "freebayes_D2.filtered.vcf.gz", chromlist = chromlist,windowSize = 1e+05, binwidth = 10)
![2](https://user-images.githubusercontent.com/93121277/163991898-568d96de-9773-4da1-a2e2-e22edd470ced.png)


VCFstat::ChromRO(vcf = "freebayes_D2.filtered.vcf.gz", chromlist = chromlist,windowSize = 1e+05, binwidth = 10)
![3](https://user-images.githubusercontent.com/93121277/163991973-cc2743e1-49df-4a59-a6b0-ef66bafbd07d.png)


VCFstat::ChromAO(vcf = "freebayes_D2.filtered.vcf.gz", chromlist = chromlist,windowSize = 1e+05, binwidth = 10)
![4](https://user-images.githubusercontent.com/93121277/163992003-8ca07920-803d-4ea0-91a6-27cb9baf5318.png)



VCFstat::ChromMQM(vcf = "freebayes_D2.filtered.vcf.gz", chromlist = chromlist,windowSize = 1e+05, binwidth = 10)

![5](https://user-images.githubusercontent.com/93121277/163992014-1d49cd02-6eb4-41da-97cc-3a5b4edec069.png)


VCFstat::ChromAC(vcf = "freebayes_D2.filtered.vcf.gz", chromlist = chromlist,windowSize = 1e+05, binwidth = 10)
![6](https://user-images.githubusercontent.com/93121277/163992029-85972022-0ed0-4613-b8dc-f8c521399447.png)



VCFstat::ChromAN(vcf = "freebayes_D2.filtered.vcf.gz", chromlist = chromlist,windowSize = 1e+05, binwidth = 10)

![7](https://user-images.githubusercontent.com/93121277/163992038-a8fd2a00-8e8b-4592-9aab-4836a6024763.png)


VCFstat::ChromnSNPs(vcf = "freebayes_D2.filtered.vcf.gz", chromlist = chromlist,windowSize = 1e+05, binwidth = 10)
![8](https://user-images.githubusercontent.com/93121277/163992055-ab66d757-7869-48f4-9f1d-83cdb75dd816.png)



VCFstat::FacetChromnSNPs(vcf = "freebayes_D2.filtered.vcf.gz", chromlist = chromlist,windowSize = 1e+05, ncol=10)
![9](https://user-images.githubusercontent.com/93121277/163992080-709ee512-a0c9-4deb-bc78-c8bb0eb091b9.png)



VCFstat::FacetChromQual(vcf = "freebayes_D2.filtered.vcf.gz", chromlist = chromlist,windowSize = 1e+05, ncol=10)
![10](https://user-images.githubusercontent.com/93121277/163992096-a6665a60-250d-4902-9ee1-9f164b81b658.png)



VCFstat::FacetChromDP(vcf = "freebayes_D2.filtered.vcf.gz", chromlist = chromlist,windowSize = 1e+05, ncol=10)
![11](https://user-images.githubusercontent.com/93121277/163992109-ae192335-ca3d-47cd-b1dd-573db78dd29c.png)



VCFstat::FacetChromAO(vcf = "freebayes_D2.filtered.vcf.gz", chromlist = chromlist,windowSize = 1e+05, ncol=10)
![12](https://user-images.githubusercontent.com/93121277/163992119-8272a52e-e0f2-4162-b94c-ec3c63fd8983.png)

![13](https://user-images.githubusercontent.com/93121277/163992127-da7dfea3-66e1-4a8b-9f75-e99cbe323773.png)



```


```{r MorePlots}

setwd("/home/michael/Desktop/Variants/Decompressed/")

chromlist <- c("NC_016131.3","NC_016132.3","NC_016133.3","NC_016134.3","NC_016135.3")


VCFstat::ChromQual(vcf = "freebayes~bwa~GCF_000005505.3_Brachypodium_distachyon_v3.0~all_samples~filtered-strict.vcf", chromlist = chromlist,windowSize = 1e+05, binwidth = 50,Maximum=5000)

![14](https://user-images.githubusercontent.com/93121277/163992145-1203a87a-1264-4e65-868f-ff31bae02c15.png)


VCFstat::ChromDP(vcf = "freebayes~bwa~GCF_000005505.3_Brachypodium_distachyon_v3.0~all_samples~filtered-strict.vcf", chromlist = chromlist,windowSize = 1e+05, binwidth = 10)
![15](https://user-images.githubusercontent.com/93121277/163992161-c74b85c0-6c94-4ffd-b8d9-65bceee90dc4.png)



VCFstat::ChromRO(vcf = "freebayes~bwa~GCF_000005505.3_Brachypodium_distachyon_v3.0~all_samples~filtered-strict.vcf", chromlist = chromlist,windowSize = 1e+05, binwidth = 10)

![16](https://user-images.githubusercontent.com/93121277/163992166-88f61376-3c6a-456f-bead-4f4a52af9780.png)


VCFstat::ChromAO(vcf = "freebayes~bwa~GCF_000005505.3_Brachypodium_distachyon_v3.0~all_samples~filtered-strict.vcf", chromlist = chromlist,windowSize = 1e+05, binwidth = 10)

![17](https://user-images.githubusercontent.com/93121277/163992184-ac8ad14b-1143-40e7-866e-296b1a6977d3.png)


VCFstat::ChromMQM(vcf = "freebayes~bwa~GCF_000005505.3_Brachypodium_distachyon_v3.0~all_samples~filtered-strict.vcf", chromlist = chromlist,windowSize = 1e+05, binwidth = 10)
![18](https://user-images.githubusercontent.com/93121277/163992194-db12fecd-816e-4178-b943-2b3b8ef57bbe.png)



VCFstat::ChromAC(vcf = "freebayes~bwa~GCF_000005505.3_Brachypodium_distachyon_v3.0~all_samples~filtered-strict.vcf", chromlist = chromlist,windowSize = 1e+05, binwidth = 10)
![19](https://user-images.githubusercontent.com/93121277/163992229-ecd45187-5527-4605-9fbc-6d36576c7631.png)



VCFstat::ChromAN(vcf = "freebayes~bwa~GCF_000005505.3_Brachypodium_distachyon_v3.0~all_samples~filtered-strict.vcf", chromlist = chromlist,windowSize = 1e+05, binwidth = 10)
![20](https://user-images.githubusercontent.com/93121277/163992243-46ef5458-dca7-4997-8946-f9218e477201.png)




VCFstat::ChromnSNPs(vcf = "freebayes~bwa~GCF_000005505.3_Brachypodium_distachyon_v3.0~all_samples~filtered-strict.vcf", chromlist = chromlist,windowSize = 1e+05, binwidth = 10)

![21](https://user-images.githubusercontent.com/93121277/163992252-21f96161-65eb-4af5-9cb7-1f62e60347c4.png)



VCFstat::FacetChromnSNPs(vcf = "freebayes~bwa~GCF_000005505.3_Brachypodium_distachyon_v3.0~all_samples~filtered-strict.vcf", chromlist = chromlist,windowSize = 1e+05, ncol=10)
![22](https://user-images.githubusercontent.com/93121277/163992268-8dc9d71c-ede5-487b-88be-804a240e7672.png)



VCFstat::FacetChromQual(vcf = "freebayes~bwa~GCF_000005505.3_Brachypodium_distachyon_v3.0~all_samples~filtered-strict.vcf", chromlist = chromlist,windowSize = 1e+05, ncol=10)

![23](https://user-images.githubusercontent.com/93121277/163992280-f49b0119-ac4d-424d-b6da-4b3a003cb261.png)



VCFstat::FacetChromDP(vcf = "freebayes~bwa~GCF_000005505.3_Brachypodium_distachyon_v3.0~all_samples~filtered-strict.vcf", chromlist = chromlist,windowSize = 1e+05, ncol=10)

![24](https://user-images.githubusercontent.com/93121277/163992294-353674eb-d71e-4058-bde1-ebd5b83eb4e8.png)



VCFstat::FacetChromAO(vcf = "freebayes~bwa~GCF_000005505.3_Brachypodium_distachyon_v3.0~all_samples~filtered-strict.vcf", chromlist = chromlist,windowSize = 1e+05, ncol=10)

![25](https://user-images.githubusercontent.com/93121277/163992323-54ffd134-b431-40ca-9a0b-8da38dbd901c.png)

![26](https://user-images.githubusercontent.com/93121277/163992348-46a0722c-2bb8-4e4c-9701-1f09c44b96a0.png)


```
