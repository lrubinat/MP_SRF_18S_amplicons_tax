---
title: "MP_SRF_18S_amplicons_tax"
author: "lrubinat"
date: "05/03/2016"
output:
  html_document:
    theme: united
    toc: yes
  pdf_document:
    highlight: zenburn
    toc: yes
---

<!--- INITIALIZATION
```{r, echo=FALSE}
#error hook to kill knitr in case of errors
library(knitr)
knit_hooks$set(error = function(x, options) stop(x))
opts_chunk$set(cache=TRUE, autodep=TRUE)
```
--->

# 1) 18S amplicons

## 1) Data overview

Let's read the dataset and remove the samples containing less than 8522 reads:

``` {r load_data, echo=FALSE, message=FALSE}
setwd("/home/laura/Documents/TFM/genwork/data_analysis/MP_18S_SRF_amplicons/MP_SRF_18S_amplicons_tax/")

#read data 
tb18_tax <- read.table(file="/home/laura/Documents/TFM/home/data/MALASPINA/Malaspina_18S_Surface/table_with_BLAST/MP_18S_SRF_MAS_BM_SILVA_classif_tab_filtered.txt", head=TRUE, fill=TRUE)

#table dimensions and format before setting column names
dim(tb18_tax) # 47144   132
tb18_tax[1:5,1:5]

#row names = OTU name (option A)
row.names(tb18_tax)<-tb18_tax[,1]

#row names = row number (option B)
#rownames(otu_tb18) <- 1:nrow(otu_tb18)

tb18_tax<-tb18_tax[,-1]
tb18_tax[is.na(tb18_tax)]<-0

dim(tb18_tax)
tb18_tax[1:5,1:5]

#table with taxonomi classification alone
tb18_class <- tb18_tax[,125:131]
dim(tb18_class)
tb18_class[1:5,1:7]

#table with occurence data alone
tb18_tax_occur <- tb18_tax[, 1:124]
dim(tb18_tax_occur) # 47144   124
tb18_tax_occur[1:5,1:5]

amplicons_per_sample_tb18<-colSums(tb18_tax_occur)
amplicons_per_sample_tb18[which(colSums(tb18_tax_occur)<8522)]
# sample st054 has less than 8522 reads.

#remove samples with less than 8522 reads
tb18_tax_occur_min8522 <- tb18_tax_occur[,colSums(tb18_tax_occur) >= 8522]
dim(tb18_tax_occur_min8522)

#remove samples with omitted in MP_SRF_18S_amplicons_tax dataset (so that we can compare the relative abundance of 16S and 18S OTUs considering the same samples)
tb18_tax_occur_min8522<-subset(tb18_tax_occur_min8522, select=-c(st122, st124, st137, st144))
dim(tb18_tax_occur_min8522) #47144   119
```

Table dimensions and content outline:

```{r starting_dataset, echo=FALSE}
dim(tb18_tax_occur_min8522)
tb18_tax_occur_min8522[1:5,1:5]
```

Minimum number of reads per station:

```{r reads_per_sample_overview1, echo=1}
min(colSums(tb18_tax_occur_min8522)) 
#8522
```

Maximum number of reads per station:

```{r reads_per_sample_overview2, echo=1}
max(colSums(tb18_tax_occur_min8522)) 
# max: 936570
```

Identification of station with higher number of reads:

```{r reads_per_sample_overview3, echo=TRUE}
amplicons_per_sample<-colSums(tb18_tax_occur_min8522)
amplicons_per_sample[which(colSums(tb18_tax_occur_min8522)>900000)]
```

Overall reads per sample:

``` {r reads_per_sample_overview4, echo=FALSE}
plot(sort(colSums(tb18_tax_occur_min8522)), pch=19, xlab="sample", ylab="reads per sample", cex=0.9)
```


## 2) Normalization

Let's normalize the original dataset by randomly subsampling 8522 reads in each station:

``` {r species_richness_rarefaction1, echo=TRUE}
library(vegan)
tb18_tax_occur_min8522_t<-t(tb18_tax_occur_min8522)
tb18_tax_occur_ss8522<-rrarefy(tb18_tax_occur_min8522_t, 8522)
```

The normalized table shows the following dimensions and format:

```{r species_richness_rarefaction2, echo=FALSE}
dim(tb18_tax_occur_ss8522)
tb18_tax_occur_ss8522[1:5,1:5]
```

Its content fits with the expected normalization values (8522 reads per station):

``` {r species_richness_rarefaction3, echo=TRUE}
rowSums(tb18_tax_occur_ss8522)
```

Let's check out how many OTUs don't appear in the new table:

```{r species_richness_rarefaction4, echo=1:5}
length(which(colSums(tb18_tax_occur_ss8522)==0)) 
```

There are 20237 OTUs that don't show any occurrence in the normalized data. Let's remove them from the table and take a look at its final dimensions:

```{r species_richness_rarefaction5, echo=1:3}
tb18_tax_occur_ss8522_no_cero<-tb18_tax_occur_ss8522[,-(which(colSums(tb18_tax_occur_ss8522)==0))]
dim(tb18_tax_occur_ss8522_no_cero)

#the final dimensions of the normalized table are 119 26907.
#26907 + 20237 = 47144
```

Datasets summary:
dim(tb18_tax) --> 47144   131
dim(tb18_tax_occur) --> 47144   124
dim(tb18_tax_occur_ss8522_no_cero) --> 119 26907



## 3) General community analysis

### 3.1) Richness and evenness (Shannon index)

```{r shannon_index1, echo=FALSE}
tb18_tax_occur_ss8522_div <- diversity(tb18_tax_occur_ss8522_no_cero, index="shannon")
```

Most of the samples take Shannon Index values between 5.5 and 6:

```{r shannon_index2, echo=FALSE}
boxplot(tb18_tax_occur_ss8522_div, pch=19, main="Shannon's index of diversity")
plot(sort(tb18_tax_occur_ss8522_div), pch=19, main="Shannon's index of diversity")
```

### 3.2) Richness: OTU number

```{r richness_otu_no1, echo=FALSE}
OTUs_per_sample_18S_tax_occur_ss8522<-specnumber(tb18_tax_occur_ss8522_no_cero)
```

Lowest number of OTUs per sample:

```{r richness_otu_no2, echo=FALSE}
min(OTUs_per_sample_18S_tax_occur_ss8522)
```

Maximum number of OTUs per sample:

```{r richness_otu_no3, echo=FALSE}
max(OTUs_per_sample_18S_tax_occur_ss8522)
```

In most of the samples, we can identify between 1000 and 1300 OTUs:

```{r richness_otu_no4, echo=TRUE}
plot(sort(OTUs_per_sample_18S_tax_occur_ss8522), pch=19)
boxplot(OTUs_per_sample_18S_tax_occur_ss8522, pch=19)
```

### 3.3) Index of evenness

#### 3.3.1) Pielou's index

```{r pielou_index_of_evenness1, echo=TRUE}
pielou_evenness_18S_tax_occur_ss8522 <- tb18_tax_occur_ss8522_div/log(OTUs_per_sample_18S_tax_occur_ss8522)
```

The Pielou index (constrained between 0 and 1) takes values closer to 1 as the variation of species proportion in a sample decreases. Most of the samples get values between 0.8 and 0.9, meaning that the numerical composition of different OTUs in a sample is highly similar:

```{r pielou_index_of_evenness2, echo=TRUE}
plot(sort(pielou_evenness_18S_tax_occur_ss8522), pch=19)
boxplot(pielou_evenness_18S_tax_occur_ss8522, pch=19)
```

The OTU_2, with 27292 reads, is the most abundant in the overall dataset:

```{r OTUs_overall_abundance, echo=TRUE}
head(sort(colSums(tb18_tax_occur_ss8522_no_cero), decreasing=T), n=10L)
```

Most of the OTUs show very few occurrences; the plot suggests that we will probably be able to identify a significant ammount of rare otus:

```{r OTUs_overall_abundance2, echo=TRUE}
plot(log(sort(colSums(tb18_tax_occur_ss8522_no_cero), decreasing=T)), pch=19)
```

<!---
#### 3.3.2) Sads
library(sads)
?sads
--->

### 3.4) Abundance Models
#### 3.4.1) Rank-Abundance or Dominance/Diversity Model ("radfit")

The OTUs abundance distribution fits relativelly close to log-normal model. 

```{r radfit, echo=FALSE}
#?radfit
#otu_tb18_t[1:5,1:5]

tb18_tax_occur_min8522_radfit<-radfit(colSums(tb18_tax_occur_min8522_t))
plot(tb18_tax_occur_min8522_radfit)
```

#### 3.4.2) Preston's Lognormal Model

According to Preston's lognormal model fit into species frequencies groups, we're missing ~3226 species:

```{r preston_model1, echo=T}
tb18_tax_occur_ss8522_prestonfit<-prestonfit(colSums(tb18_tax_occur_min8522_t))
plot(tb18_tax_occur_ss8522_prestonfit, main="Pooled species")

veiledspec(tb18_tax_occur_ss8522_prestonfit)
```

### COMA ### When computing Prestons lognormal model fit without pooling data into groups, we seem to miss ~3055 species:

```{r preston_model2, echo=4}
tb18_tax_occur_ss8522_dist_all<-prestondistr(colSums(tb18_tax_occur_min8522_t))
plot(tb18_tax_occur_ss8522_prestonfit, main="All malaspina")
lines(tb18_tax_occur_ss8522_dist_all, line.col="blue3")

veiledspec(tb18_tax_occur_ss8522_dist_all)
```

<!---
### 3.5) Rarefaction curve

(To be computed)

```{r rarefraction_curve, echo=TRUE}
#?rarecurve

#str(colSums(otu_tb18_t))

#otus_tb18_colsums<-colSums(otu_tb18_t)

#str(otus_tb18_colsums)
#otu_tb18_colsums<-as.matrix(otu_tb18_colsums)

#otu_tb18_colsums<-t(otu_tb18_colsums)

#otu_tb18_colsums[,1:3]

#rarecurve(otu_tb18_colsums, step = 1, 50000, xlab = "Sample Size", ylab = "OTUs", label = TRUE)
#rarecurve(colSums(otu_tb18_t), step = 1, 50000, xlab = "Sample Size", ylab = "OTUs", label = TRUE)
```
--->

### 3.6) Beta diversity

#### 3.6.1) Dissimilarity matrix using Bray-Curtis index:

The Bray-Curtis dissimilarity, constrained between 0 (minimum distance) and 1 (highest dissimilarity) allows us to quantify the differences between samples according to the composition and relative abundance of their OTUs. In our dataset, most of the samples pairs take dissimilarity values between between 0.7 and 0.8, meaning that their composition is substantially different.

```{r beta_div1, echo=FALSE}
#?vegdist
tb18_tax_occur_ss8522_no_cero.bray<-vegdist(tb18_tax_occur_ss8522_no_cero, method="bray")
boxplot(tb18_tax_occur_ss8522_no_cero.bray, main="Bray-Curtis dissimilarity matrix")
```

#### 3.6.2) Hierarchical clustering

The only relatively evident cluster we can distinguish in the dendogram stands out in the very left side of the plot. 

### COMA after SAMPLES' ### (To be done: assign Longhurst provinces information to each station and check if any of the central clusters is meaningful regarding to the samples geographical ubication)

```{r beta_div2, echo=FALSE}
#UPGMA
tb18_tax_occur_ss8522_no_cero.upgma<-hclust(tb18_tax_occur_ss8522_no_cero.bray, "average")
plot(tb18_tax_occur_ss8522_no_cero.upgma, cex=.35, main="Samples Hierarchical Clustering")
```

#### 3.6.3) Non-metric multidimensional scaling

We can identify a prominent group in the central part of the NMDS plot and a few outliers in the middle-left edge of the plot. The stress parameter takes a value below 0.3, suggesting that the plot is acceptable. 

```{r monoNMDS, echo=F}
#NMDS
tb18_tax_occur_ss8522_no_cero.nmds<-monoMDS(tb18_tax_occur_ss8522_no_cero.bray)
tb18_tax_occur_ss8522_no_cero.nmds
plot(tb18_tax_occur_ss8522_no_cero.nmds, main="monoMDs method")
```

When implementing a most robut function for computing NMDS plots, the result is quiet the same:

```{r metaNMDS, echo=F}
tb18_tax_occur_ss8522_no_cero.meta_nmds<-metaMDS(tb18_tax_occur_ss8522_no_cero.bray)
plot(tb18_tax_occur_ss8522_no_cero.meta_nmds, main="metaMDS method")
```

## 4) Geographical analysis

```{r load_geo_data, echo=F, results="hide", message=F}
#load geographical ubication of stations and sort according to otu_tb18 stations sequence.
MP_geo_18S<-read.table(file="/home/laura/Documents/TFM/home/data/MALASPINA/mp_surface_ubication.txt", sep="\t", header=T)

row.names(MP_geo_18S)<-MP_geo_18S[,1]

MP_geo_18S_sorted<-MP_geo_18S[row.names(tb18_tax_occur_ss8522_no_cero),]
dim(MP_geo_18S_sorted)
MP_geo_18S_sorted[1:5,1:5]
tb18_tax_occur_ss8522_no_cero[1:5,1:5]

#read lat-long in decimal degrees and translate into distance in km.
library(fossil)

#select only columns containing info about station, latitude and longitude.
MP_geo_18S_sorted_v2<-create.lats(MP_geo_18S_sorted, loc="sample", long="long", lat="lat")
head(MP_geo_18S_sorted)

#create a distance matrix (lower triangle) between a list of points.
geo_distances_MP_18S<-earth.dist(MP_geo_18S_sorted_v2, dist = TRUE)
head(geo_distances_MP_18S)
dim(geo_distances_MP_18S)

geo_distances_MP_18S<-as.matrix(geo_distances_MP_18S)
dim(geo_distances_MP_18S)

#geo distances dataset ready to use "geo_distances_MP_18S"
```

Working datasets:

1) Community matrix: tb18_tax_occur_ss8522_no_cero

```{r working_datasets1, echo=T}
dim(tb18_tax_occur_ss8522_no_cero)
tb18_tax_occur_ss8522_no_cero[1:5, 1:5]
```

2) Community Bray-Curtis: tb18_tax_occur_ss8522_no_cero.bray

```{r working_datasets2, echo=2}
dim(tb18_tax_occur_ss8522_no_cero.bray)
tb18_tax_occur_ss8522_no_cero.bray<-as.matrix(tb18_tax_occur_ss8522_no_cero.bray)
```

3) Stations distances in km: geo_distances_MP_18S

```{r working_datasets3, echo=T}
dim(geo_distances_MP_18S)
```

Communities quickly change their composition across geographical distances:

```{r working_datasets4, echo=T}
plot(geo_distances_MP_18S, tb18_tax_occur_ss8522_no_cero.bray, pch=19, cex=0.4, xlab="Geopgraphical distances", ylab="Bray-Curtis dissimilarities")
```

### 4.1) Mantel correlograms

Mantel statistic is -significantlly- so low, meaning that the correlation between samples dissimilarity and geographical distances is weak.

```{r mantel_correlogram1, echo=T}
mantel(geo_distances_MP_18S, tb18_tax_occur_ss8522_no_cero.bray)
```

Maximum distance between samples:

```{r mantel_correlogram2, echo=F}
max(geo_distances_MP_18S)
```

Minimum distance between samples:

```{r mantel_correlogram3, echo=F}
min(geo_distances_MP_18S)
```

Correlograms:

```{r mantel_correlogram4, echo=T}
MP_18s_ss8522_mantel_correl_by_1000km<-mantel.correlog(tb18_tax_occur_ss8522_no_cero.bray, D.geo=geo_distances_MP_18S, break.pts=seq(0,20000, by=1000))
plot(MP_18s_ss8522_mantel_correl_by_1000km)

MP_18s_ss8522_mantel_correl_by_100km<-mantel.correlog(tb18_tax_occur_ss8522_no_cero.bray, D.geo=geo_distances_MP_18S, break.pts=seq(0,20000, by=100))
plot(MP_18s_ss8522_mantel_correl_by_100km)
```

## 5) Abundance vs. occurence

```{r OTUs_mean_relative_abund, echo=F, results="hide"}
tb18_tax_occur_ss8522_no_cero[1:5,1:5]
tb18_tax_occur_ss8522_no_cero_t<-t(tb18_tax_occur_ss8522_no_cero)

colSums(tb18_tax_occur_ss8522_no_cero_t)

#local abundance percentage
tb18_tax_occur_ss8522_no_cero_t.rabund<-tb18_tax_occur_ss8522_no_cero_t/8522

colSums(tb18_tax_occur_ss8522_no_cero_t.rabund)
tb18_tax_occur_ss8522_no_cero_t.rabund[1:5,1:5]

#OTUs mean relative abundance
tb18_txa_occur_ss8522_no_cero_t.rabund_means<-rowMeans(tb18_tax_occur_ss8522_no_cero_t.rabund) 
tb18_tax_occur_ss8522_no_cero_t.rabund_means<-as.data.frame(tb18_tax_occur_ss8522_no_cero_t.rabund_means)

head(tb18_tax_occur_ss8522_no_cero_t.rabund_means)
```

```{r OTUs_occurence, echo=F, results='hide'}
tb18_tax_occur_ss8522_no_cero_t.rabund.occur<-tb18_tax_occur_ss8522_no_cero_t.rabund
tb18_tax_occur_ss8522_no_cero_t.rabund.occur[tb18_tax_occur_ss8522_no_cero_t.rabund.occur>0]<-1
tb18_tax_occur_ss8522_no_cero_t.rabund.occur[1:5,1:5] ### presence - absence table

#percentage of occurence in overall stations
tb18_tax_occur_ss8522_no_cero_t.rabund_means.occurence_perc<-as.data.frame(100*(rowSums(tb18_tax_occur_ss8522_no_cero_t.rabund.occur)/119))

str(tb18_tax_occur_ss8522_no_cero_t.rabund_means.occurence_perc)
```

```{r merge_rabund_peroccur, echo=F, results='hide'}
otu_tb18_ss8522_rabund_percoccur<-merge(tb18_tax_occur_ss8522_no_cero_t.rabund_means,tb18_tax_occur_ss8522_no_cero_t.rabund_means.occurence_perc, by="row.names")

colnames(otu_tb18_ss8522_rabund_percoccur)<-c("OTUs","mean_rabund","perc_occur")
otu_tb18_ss8522_rabund_percoccur[1:5,]

row.names(otu_tb18_ss8522_rabund_percoccur)<-otu_tb18_ss8522_rabund_percoccur[,1]
otu_tb18_ss8522_rabund_percoccur<-otu_tb18_ss8522_rabund_percoccur[,-1]
otu_tb18_ss8522_rabund_percoccur[1:5,]
```

In the following plot, we can appreciate the OTUs distribution according to their percentage of occurence and relative abundance. The red line keeps up OTUs that occur in more than 80% of the samples, the green line limits regionally rare OTUs (< 0.001%), and the blue one restricts regionally abundant OTUs (> 0.1%).

```{r abund_vs_occurence_table, echo=F}
plot(otu_tb18_ss8522_rabund_percoccur$mean_rabund,otu_tb18_ss8522_rabund_percoccur$perc_occur, log="x", pch=19, cex=0.8, xlab="Mean relative abundance", ylab="Percentage of occurence")
abline(h=80, col="red") #occurence higher than 80%
abline(v=0.00001, col="green") #rare OTUs
abline(v=0.001, col="blue") #cosmopolitan OTUs

#Conventional limits:
#Regionally rare     = 0.00001
#Regionally abundant = 0.001
```

Regionally abundant OTUs (relative abundance over 0.1%):

```{r abundant_OTUs, echo=7}
#regionally abundant
tb18_ss8522_abundant<-otu_tb18_ss8522_rabund_percoccur[otu_tb18_ss8522_rabund_percoccur$mean_rabund > 0.001,]

tb18_ss8522_abundant_sorted<-tb18_ss8522_abundant[order(tb18_ss8522_abundant$mean_rabund, tb18_ss8522_abundant$perc_occur, decreasing = T), c(1,2)]

tb18_ss8522_abundant_sorted
dim(tb18_ss8522_abundant_sorted)
```

Proportion of regionally abundant OTUs (%):

```{r abundant_OTUs2, echo=F}
#there are 83 regionally abundant OTUs.
(162/25191)*100 # = 7.29% of the OTUs are regionally abundant

#length(row.names(otu_tb18_ss8522_rabund_percoccur[otu_tb18_ss8522_rabund_percoccur$mean_rabund > 0.001,])) # 83 OTUs
#row.names(otu_tb18_ss8522_rabund_percoccur[otu_tb18_ss8522_rabund_percoccur$mean_rabund > 0.001,])
```

Cosmopolitan OTUs (relative abundance over 0.1% and occurence in more than 80% of samples):

```{r select_cosmopolitan, echo=6}
otu_tb18_ss8522_rabund_cosm<-otu_tb18_ss8522_rabund_percoccur[otu_tb18_ss8522_rabund_percoccur$mean_rabund > 0.001,]
otu_tb18_ss8522_rabund_poccur_cosm<-otu_tb18_ss8522_rabund_cosm[otu_tb18_ss8522_rabund_cosm$perc_occur > 80,]
otu_tb18_ss8522_cosmop_sorted<-otu_tb18_ss8522_rabund_poccur_cosm[order(otu_tb18_ss8522_rabund_poccur_cosm$perc_occur, otu_tb18_ss8522_rabund_poccur_cosm$mean_rabund, decreasing = T), c(1,2)]

otu_tb18_ss8522_cosmop_sorted
dim(otu_tb18_ss8522_cosmop_sorted)
```

Proportion of cosmopolitan OTUs (%):

```{r percentage_cosmopolitan, echo=F}
(64/25191)*100
```

Number and proportion (%) of rare OTUs:

```{r rare_OTUs, echo=1}
dim(otu_tb18_ss8522_rabund_percoccur[otu_tb18_ss8522_rabund_percoccur$mean_rabund < 0.00001 & otu_tb18_ss8522_rabund_percoccur$mean_rabund >0,])
 
(18551/25191)*100 # = 28.56% of the OTUs are regionally rare
```

<!---
```{r otu_col_chech, echo = T}
dim(otu_tb18_ss8522_rabund_percoccur)
dim(otu_tb18_min5000_v2)
otu_tb18_min5000_v2[1:5,1:5]

#rare OTUs:
# length(row.names(otu_tb16_ss8522_rabund_percoccur[otu_tb16_ss8522_rabund_percoccur$mean_rabund < 0.00001 & otu_tb16_ss8522_rabund_percoccur$mean_rabund >0 ,])) # 325 OTUs
```
--->


## 6) Taxonomic composition analysis

### COMA Let's ### Let's add the taxonomic classification to the left OTUs by merging "tb18_tax_occur_ss8522_no_cero" with "tb18_tax":

```{r merge_tables, echo=FALSE}
tb18_tax_occur_ss8522_no_cero_t<-t(tb18_tax_occur_ss8522_no_cero)
tb18_ss8522_tax<-merge(tb18_tax_occur_ss8522_no_cero_t,tb18_class, by="row.names")

dim(tb18_ss8522_tax)
tb18_ss8522_tax[1:5,1:5]

#fix OTU_no as new row
rownames(tb18_ss8522_tax)=tb18_ss8522_tax$Row.names

#add OTU_no as rowname
rownames.tb18_ss8522_tax<-tb18_ss8522_tax[,1]
tb18_ss8522_tax<-tb18_ss8522_tax[,-1]
#colnames(tb18_ss8522_tax, do.NULL=F)

dim(tb18_ss8522_tax)
tb18_ss8522_tax[1:5, 1:5]

#sort by OTU_no (split rowname, introduce no. into column "OTU_no" and sort)
tb18_ss8522_tax["OTU_no"] <- NA
tb18_ss8522_tax$OTU_no <- sapply(strsplit(rownames(tb18_ss8522_tax),split= "\\_"),'[',2)
tb18_ss8522_tax$OTU_no <- as.numeric(as.character(tb18_ss8522_tax$OTU_no))
tb18_ss8522_tax_sorted<-tb18_ss8522_tax[order(tb18_ss8522_tax$OTU_no, decreasing = FALSE), ]

dim(tb18_ss8522_tax_sorted)
tb18_ss8522_tax_sorted[1:5,1:5]
```

```{r select_phototrophs, echo=FALSE}
tb18_phototrophs <- tb18_ss8522_tax_sorted[which(tb18_ss8522_tax_sorted$MAS_plus_BM_plus_SILVA_class != "NA"),]
dim(tb18_phototrophs)
tb18_phototrophs[1:5,123]
```

```{r aggregate, echo=FALSE}

class_summary_reads_per_class<-aggregate(rowSums(tb18_phototrophs[1:123]), list(tb18_phototrophs$MAS_plus_BM_plus_SILVA_class), sum)
# count the different groups

class_summary_otus_per_class<-aggregate(rowSums(tb18_phototrophs[1:123]), list(tb18_phototrophs$MAS_plus_BM_plus_SILVA_class), length)

attach(class_summary_reads_per_class)
class_summary_reads_per_class_order<-class_summary_reads_per_class[order(-x),]
detach(class_summary_reads_per_class)
class_summary_reads_per_class_order

attach(class_summary_otus_per_class)
class_summary_otus_per_class_order<-class_summary_otus_per_class[order(-x),]
detach(class_summary_otus_per_class)
class_summary_otus_per_class_order


#class_summary_reads<-aggregate(sum~class, data=otutab_full_wTax, FUN="sum") 
# sum reads different groups
```


# 2) 16S amplicons

<!--- INITIALIZATION
```{r, echo=FALSE}
#error hook to kill knitr in case of errors
library(knitr)
knit_hooks$set(error = function(x, options) stop(x))
opts_chunk$set(cache=TRUE, autodep=TRUE)
```
--->


# 1) Data overview

``` {r load_data_2, echo=FALSE, message=FALSE}
setwd("/home/laura/Documents/TFM/genwork/data_analysis/MP_16S_SRF_amplicons/MP_16S_amplicons_tax/")

#read data 
tb16_tax <- read.table(file="/home/laura/Documents/TFM/home/data/MALASPINA/Malaspina_16S_Surface/table_with_BLAST/MP_16S_SRF_tax_classif_tab_filtered.txt", head=TRUE, fill=TRUE)

#table dimensions and format before setting column names
dim(tb16_tax) # 54979   132
tb16_tax[1:5,1:5]

#row names = OTU name (option A)
row.names(tb16_tax)<-tb16_tax[,1]

#row names = row number (option B)
#rownames(otu_tb16) <- 1:nrow(otu_tb16)

tb16_tax<-tb16_tax[,-1]
tb16_tax[is.na(tb16_tax)]<-0

dim(tb16_tax)
tb16_tax[1:5,1:5]

#table with taxonomi classification alone
tb16_class <- tb16_tax[,125:129]
dim(tb16_class)
tb16_class[1:5,1:5]

#table with occurence data alone
tb16_tax_occur <- tb16_tax[, 1:124]
dim(tb16_tax_occur) # 993 124
tb16_tax_occur[1:5,1:5]

amplicons_per_sample_tb16<-colSums(tb16_tax_occur)
amplicons_per_sample_tb16[which(colSums(tb16_tax_occur)<5844)]
# samples st122, st124, st137 and st144 have less than 5844 reads.

#remove samples with less than 5844 reads
tb16_tax_occur_min5844 <- tb16_tax_occur[,colSums(tb16_tax_occur) >= 5844]
dim(tb16_tax_occur_min5844)


#remove samples with omitted in MP_SRF_18S_amplicons_tax dataset (so that we can compare the relative abundance of 16S and 16S OTUs considering the same samples)
tb16_tax_occur_min5844<-subset(tb16_tax_occur_min5844, select=-c(st054))
dim(tb16_tax_occur_min5844)
```

Table dimensions and content outline:

```{r starting_dataset_2, echo=FALSE}
dim(tb16_tax_occur_min5844)
tb16_tax_occur_min5844[1:5,1:5]
```

Minimum number of reads per station:

```{r reads_per_sample_overview1_2, echo=1}
min(colSums(tb16_tax_occur_min5844)) 
#5844
```

Maximum number of reads per station:

```{r reads_per_sample_overview2_2, echo=1}
max(colSums(tb16_tax_occur_min5844)) 
# max: 156413
```

Identification of station with higher number of reads:

```{r reads_per_sample_overview3_2, echo=TRUE}
amplicons_per_sample<-colSums(tb16_tax_occur_min5844)
amplicons_per_sample[which(colSums(tb16_tax_occur_min5844)>150000)]
```

Overall reads per sample:

``` {r reads_per_sample_overview4_2, echo=FALSE}
plot(sort(colSums(tb16_tax_occur_min5844)), pch=19, xlab="sample", ylab="reads per sample", cex=0.9)
```


## 2) Normalization

Let's normalize the original dataset by randomly subsampling 5844 reads in each station:

``` {r species_richness_rarefaction1_2, echo=TRUE}
library(vegan)
tb16_tax_occur_min5844_t<-t(tb16_tax_occur_min5844)
tb16_tax_occur_ss5844<-rrarefy(tb16_tax_occur_min5844_t, 5844)
```

The normalized table shows the following dimensions and format:

```{r species_richness_rarefaction2_2, echo=FALSE}
dim(tb16_tax_occur_ss5844)
tb16_tax_occur_ss5844[1:5,1:5]
```

Its content fits with the expected normalization values (5844 reads per station):

``` {r species_richness_rarefaction3_2, echo=TRUE}
rowSums(tb16_tax_occur_ss5844)
```

Let's check out how many OTUs don't appear in the new table:

```{r species_richness_rarefaction4_2, echo=1:5}
length(which(colSums(tb16_tax_occur_ss5844)==0)) 
```

There are 54 OTUs that don't show any occurrence in the normalized data. Let's remove them from the table and take a look at its final dimensions:

```{r species_richness_rarefaction5_2, echo=1:3}
tb16_tax_occur_ss5844_no_cero<-tb16_tax_occur_ss5844[,-(which(colSums(tb16_tax_occur_ss5844)==0))]
dim(tb16_tax_occur_ss5844_no_cero)

#the final dimensions of the normalized table are 123 27067.
#939 + 54 = 993
```

Datasets summary:
dim(tb16_tax) --> 993 129
dim(tb16_tax_occur) --> 993 124
dim(tb16_tax_occur_ss5844_no_cero) --> 119 939

Let's add the taxonomic classification to the left OTUs by merging "tb16_tax_occur_ss5844_no_cero" with "tb16_tax".

```{r merge_tables_2, echo=FALSE}
tb16_tax_occur_ss5844_no_cero_t<-t(tb16_tax_occur_ss5844_no_cero)
tb16_ss5844_tax<-merge(tb16_tax_occur_ss5844_no_cero_t,tb16_class, by="row.names")

dim(tb16_ss5844_tax)
tb16_ss5844_tax[1:5,1:5]
colSums(tb16_ss5844_tax[,2:120])

#fix OTU_no as new row
rownames(tb16_ss5844_tax)=tb16_ss5844_tax$Row.names

#add OTU_no as rowname
rownames.tb16_ss5844_tax<-tb16_ss5844_tax[,1]
tb16_ss5844_tax<-tb16_ss5844_tax[,-1]
#colnames(tb16_ss5844_tax, do.NULL=F)

dim(tb16_ss5844_tax)
tb16_ss5844_tax[1:5, 1:5]

#sort by OTU_no (split rowname, introduce no. into column "OTU_no" and sort)
tb16_ss5844_tax["OTU_no"] <- NA
tb16_ss5844_tax$OTU_no <- sapply(strsplit(rownames(tb16_ss5844_tax),split= "\\_"),'[',2)
tb16_ss5844_tax$OTU_no <- as.numeric(as.character(tb16_ss5844_tax$OTU_no))
tb16_ss5844_tax_sorted<-tb16_ss5844_tax[order(tb16_ss5844_tax$OTU_no, decreasing = FALSE), ]

dim(tb16_ss5844_tax_sorted)
tb16_ss5844_tax_sorted[1:5,124:125]
```

```{r select_phototrophs_2, echo=FALSE}
tb16_bacteria <- tb16_ss5844_tax_sorted[which(tb16_ss5844_tax_sorted$taxonomic_classes != "Bolidophyceae" & tb16_ss5844_tax_sorted$taxonomic_classes != "Chlorarachniophyceae" & tb16_ss5844_tax_sorted$taxonomic_classes != "Cryptophyceae" & tb16_ss5844_tax_sorted$taxonomic_classes != "Pelagophyceae" & tb16_ss5844_tax_sorted$taxonomic_classes != "Rappemonad" & tb16_ss5844_tax_sorted$taxonomic_classes != "Chrysophyceae" & tb16_ss5844_tax_sorted$taxonomic_classes != "Bacillariophyta" & tb16_ss5844_tax_sorted$taxonomic_classes != "Dictyochophyceae" & tb16_ss5844_tax_sorted$taxonomic_classes != "Mamiellophyceae" & tb16_ss5844_tax_sorted$taxonomic_classes != "Prasinophyceae" & tb16_ss5844_tax_sorted$taxonomic_classes != "Prymnesiophyceae" & tb16_ss5844_tax_sorted$taxonomic_classes != "other_plastids"),]

tb16_protists <- tb16_ss5844_tax_sorted[which(tb16_ss5844_tax_sorted$taxonomic_classes != "bacteria" & tb16_ss5844_tax_sorted$taxonomic_classes != "other_cyanob" & tb16_ss5844_tax_sorted$taxonomic_classes != "Prochlorococcus" & tb16_ss5844_tax_sorted$taxonomic_classes != "Synechococcus"),]

dim(tb16_protists)
dim(tb16_bacteria)
tb16_protists[1:5,119:125]
tb16_bacteria[1:5,124:125]
```


```{r aggregate_2, echo=TRUE}

class_summary_reads_per_class_16S<-aggregate(rowSums(tb16_protists[1:119]), list(tb16_protists$taxonomic_classes), sum)
# count the different groups

class_summary_otus_per_class_16S<-aggregate(rowSums(tb16_protists[1:119]), list(tb16_protists$taxonomic_classes), length)

attach(class_summary_reads_per_class_16S)
class_summary_reads_per_class_16S_order<-class_summary_reads_per_class_16S[order(-x),]
detach(class_summary_reads_per_class_16S)
class_summary_reads_per_class_16S_order

attach(class_summary_otus_per_class_16S)
class_summary_otus_per_class_16S_order<-class_summary_otus_per_class_16S[order(-x),]
detach(class_summary_otus_per_class_16S)
class_summary_otus_per_class_16S_order


#class_summary_reads<-aggregate(sum~class, data=otutab_full_wTax, FUN="sum") 
# sum reads different groups
```


