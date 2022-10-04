1_HealthyHumanFlagellome
================
Andrea Borbon

## Load libraries

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──

    ## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
    ## ✔ tibble  3.1.7     ✔ dplyr   1.0.9
    ## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
    ## ✔ readr   2.1.2     ✔ forcats 0.5.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

``` r
library(ggplot2)
library(scales)
```

    ## 
    ## Attaching package: 'scales'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     discard

    ## The following object is masked from 'package:readr':
    ## 
    ##     col_factor

``` r
library(tidytext)
library(taxonomizr)
library(conflicted)
library(ggsci)
```

## set working directory

``` r
knitr::opts_knit$set(root.dir = "/ebio/abt3_projects/small_projects/aborbon/TLR5/SilentFlagellins_code/")
```

#Declare conflict preferences

``` r
conflict_prefer("filter","dplyr")
```

    ## [conflicted] Will prefer dplyr::filter over any other package

## Metagenomes mapped to flagellins for candidates selection

This notebook describes the characterization of top abundant flagellins
present in metagenomic samples in the dataset of Lloyd-Price et al.,
2019. This dataset comprised 270 samples from healthy adults. (indent 1
tab)

## 1. Import metadata and mapping stats

``` r
#Import sequencing qc stats (from multiFastQC)
lloyd_qc=tbl_df(read.delim("data/StatsLloyd_fastqc.txt",header=T))
```

    ## Warning: `tbl_df()` was deprecated in dplyr 1.0.0.
    ## Please use `tibble::as_tibble()` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.

``` r
names(lloyd_qc)[names(lloyd_qc)=="site_sub_coll"]<-"Sample"

#Import metadata from Lloyd-Price samples
metadataLloyd=tbl_df(read.delim("data/LloydPrice_metadata.txt",sep="\t",header=T))%>%
  mutate(Sample=str_remove_all(Sample,"_MGX"))

#Import statistics of mapped reads (this is only the mapped reads to flagellins per sample)
mappedReadsStats=tbl_df(read.delim("data/mapping_stats.txt", sep="", header=T))%>%
  rename(Sample="site_sub_coll")%>%
  mutate(Sample=str_remove_all(Sample,"_MGX"))

#Import diamond blastx results - reads mapped to flagellin accessions  
flagellin_diamond_lloyd=tbl_df(read.delim("data/mergedFlagellinLloyd.txt",header=F,sep="")%>%
  mutate(Study="LloydPrice_2019"))
names(flagellin_diamond_lloyd)=c("Sample","qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","Study")

head(flagellin_diamond_lloyd)
```

    ## # A tibble: 6 × 14
    ##   Sample  qseqid sseqid pident length mismatch gapopen qstart  qend sstart  send
    ##   <chr>   <chr>  <chr>   <dbl>  <int>    <int>   <int>  <int> <int>  <int> <int>
    ## 1 C3001C… SRR59… WP_05…  100       33        0       0     99     1    413   445
    ## 2 C3001C… SRR59… DAA95…  100       33        0       0      1    99    148   180
    ## 3 C3001C… SRR59… SCJ35…   75.8     33        8       0      1    99    218   250
    ## 4 C3001C… SRR59… CDC20…  100       33        0       0      3   101     46    78
    ## 5 C3001C… SRR59… CDD73…  100       33        0       0      3   101    354   386
    ## 6 C3001C… SRR59… WP_02…  100       33        0       0     99     1    103   135
    ## # … with 3 more variables: evalue <dbl>, bitscore <dbl>, Study <chr>

## 2. Join metadata, mapped reads and sequencing QC statistics

###### Read counts table for all samples (controls + IBD)

``` r
flagellin_lloyd_total=inner_join(lloyd_qc,metadataLloyd,by="Sample")%>%
  inner_join(mappedReadsStats,by="Sample")%>%
  full_join(flagellin_diamond_lloyd)%>%
  group_by(sseqid)%>%
  summarise(mappedReads=n())%>%
  rename(Accession=sseqid)
```

    ## Joining, by = "Sample"

###### Subset samples from healthy subjects

``` r
#Merge all mapping results and metadata. Keep only samples from healthy subjects (Diagnosis==Control)

flagellin_lloyd_metadata_control=inner_join(lloyd_qc,metadataLloyd,by="Sample")%>%
  inner_join(mappedReadsStats,by="Sample")%>%
  full_join(flagellin_diamond_lloyd)%>%
  filter(Condition=="Control")
```

    ## Joining, by = "Sample"

## 3. Transform mapping hits results into read counts per accession

Only for the subset of samples from healthy subjects

``` r
lloyd.healthy<-flagellin_lloyd_metadata_control%>%
  group_by(sseqid)%>%
  summarise(mappedReads=n())%>%
  rename(Accession=sseqid)

head(lloyd.healthy)
```

    ## # A tibble: 6 × 2
    ##   Accession  mappedReads
    ##   <chr>            <int>
    ## 1 AAA23797.1           6
    ## 2 AAA23798.1           8
    ## 3 AAA23799.1           2
    ## 4 AAA26411.1         419
    ## 5 AAB17947.1          57
    ## 6 AAB82613.1           1

## 4. Descriptive stats on mapping results

Descriptive stats on the read counts.

``` r
#Read counts table - subset healthy
read.count= lloyd.healthy

stats.mapping=data.frame(dataset="LloydPrice",
                          "mean"=mean(read.count$mappedReads),
                          "median"=median(read.count$mappedReads),
                          "min"=min(read.count$mappedReads),
                          "max"=max(read.count$mappedReads))
stats.mapping %>% print()
```

    ##      dataset     mean median min   max
    ## 1 LloydPrice 239.0225      9   1 16608

### Summary statistics on read count filters

#### Option 1: Filtering using the median of mapped reads

``` r
summary_sseqid_healthy.median<-lloyd.healthy%>%
   filter(mappedReads>median(mappedReads))

  #Summary statistics
  nrow(summary_sseqid_healthy.median) 
```

    ## [1] 1263

``` r
  min(summary_sseqid_healthy.median$mappedReads) 
```

    ## [1] 10

``` r
  max(summary_sseqid_healthy.median$mappedReads) 
```

    ## [1] 16608

``` r
#Plot mapped reads distribution using the median of mapped reads as cut-off
ggplot(data=summary_sseqid_healthy.median,aes(x=mappedReads))+
  geom_histogram(bins=50)+
  geom_vline(aes(xintercept = median(mappedReads)),col='red',size=1)
```

![](1_HealthyHumanFlagellome_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

#### Option 2: Filtering using the mean of mapped reads

``` r
summary_sseqid_healthy.mean<-lloyd.healthy%>%
   filter(mappedReads>mean(mappedReads))

  #Summary statistics
  nrow(summary_sseqid_healthy.mean) 
```

    ## [1] 461

``` r
  min(summary_sseqid_healthy.mean$mappedReads) 
```

    ## [1] 242

``` r
  max(summary_sseqid_healthy.mean$mappedReads)
```

    ## [1] 16608

``` r
#Plot mapped reads distribution using the mean of mapped reads as cut-off
ggplot(data=summary_sseqid_healthy.mean,aes(x=mappedReads))+
  geom_histogram(bins=50)+
  geom_vline(aes(xintercept = median(mappedReads)),col='red',size=1)
```

![](1_HealthyHumanFlagellome_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

## 5. Filtering flagellins

#### Filtering flagellins

Flagellins with a read count below the median of mapped reads/flagellin
will be discarded.

``` r
#Summarize flagellins and filter using the mean of mappedReads

subsetHealthy<-flagellin_lloyd_metadata_control%>%
  group_by(sseqid)%>%
  summarise(mappedReads=n())%>%
  rename(Accession=sseqid)%>%
  filter(mappedReads>median(mappedReads))

head(subsetHealthy)
```

    ## # A tibble: 6 × 2
    ##   Accession  mappedReads
    ##   <chr>            <int>
    ## 1 AAA26411.1         419
    ## 2 AAB17947.1          57
    ## 3 AAD28520.2          70
    ## 4 AAF71896.1          16
    ## 5 AAL30165.1          15
    ## 6 AAN77105.1          25

## 6. Taxonomic annotation of mapped flagellins

#### 6.1. Import GTDB metadata and taxonomy files

``` r
#Import metadata of GTDB version 202
gtdb_bac_metadata_v202=read.delim("data/bac120_metadata_r202.tsv")

#Import taxonomy of GTDB version 202
gtdb_bac.tax_v202=read.delim("data/bac120_taxonomy_r202.tsv",header=F)%>%
  rename(accession=V1)%>%
  separate(V2,into=c("Domain","Phylum","Class","Order","Family","Genus","Species"),sep=";")

#Merge taxonomy and metadata files
taxid.to.gtdb_taxonomy.v202=inner_join(gtdb_bac.tax_v202,gtdb_bac_metadata_v202)%>%
  select(Domain,Phylum,Class,Order,Family,Genus,Species,ncbi_taxid,accession)
```

    ## Joining, by = "accession"

#### 6.2. Import esearch-annotated accessions

``` r
#This was produced by using esearch to find the assembly associated to each protein accession in the IPG database, and subsequently searching said assembly ID into the Assembly database.
metadata_esearch=read_table("data/finalMetadataEsearchOnlyids.txt")%>%
  rename(ncbi_genbank_assembly_accession="Assembly")
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   Accession = col_character(),
    ##   Biosample = col_character(),
    ##   Assembly = col_character()
    ## )

``` r
head(metadata_esearch)
```

    ## # A tibble: 6 × 3
    ##   Accession      Biosample    ncbi_genbank_assembly_accession
    ##   <chr>          <chr>        <chr>                          
    ## 1 WP_108001953.1 SAMN08769590 GCF_003046475.1                
    ## 2 WP_075123400.1 SAMN05942569 GCF_001921205.1                
    ## 3 CUH95360.1     SAMEA3539092 GCA_001298735.1                
    ## 4 OWG58734.1     SAMN07137396 GCA_002193165.1                
    ## 5 OYW63100.1     SAMN06622388 GCA_002279885.1                
    ## 6 WP_092563040.1 SAMN05661086 GCF_900112775.1

#### 6.4. Prepares NCBI database to obtain taxids from accessions

``` r
#This step will take long to process
#prepareDatabase(sqlFile='accessionTaxa.sql',tmpDir=".",getAccessions=T)
```

#### 6.5. Find taxids of flagellin accessions with taxonomizr

``` r
#Get taxids of flagellins in subset of healthy subjects above abundance cutoff
taxids.subsethealthy=subsetHealthy %>%
  mutate(ncbi_taxid=
  accessionToTaxa((as.character(subsetHealthy$Accession)),"/ebio/abt3_projects/databases_no-backup/NCBI_accession2taxid/accessionTaxa.sql",version="version"))
```

#### 6.6. Assign taxonomy to flagellins found in samples from healthy subjects

``` r
#Filtered list (by read counts) 
taxonomy.lloyd.healthy=left_join(taxids.subsethealthy,taxid.to.gtdb_taxonomy.v202,by="ncbi_taxid")%>%
  distinct(Accession,.keep_all=T)%>%
  select(contains(c("Accession","Domain","Phylum","Class","Order","Family","Genus","Species","mappedReads")))%>%
  rename(gtdb_assembly="accession")%>%
  filter(!is.na(Phylum))
  
no.taxonomy.lloyd.healthy=anti_join(taxids.subsethealthy,taxid.to.gtdb_taxonomy.v202,by="ncbi_taxid")%>%
  inner_join(metadata_esearch,by="Accession")%>%
  inner_join(taxid.to.gtdb_taxonomy.v202%>%mutate(ncbi_genbank_assembly_accession=str_remove_all(accession,"RS_")),by="ncbi_genbank_assembly_accession")%>%
  select(contains(c("Accession","Domain","Phylum","Class","Order","Family","Genus","Species","mappedReads")))%>%
  select(-ncbi_genbank_assembly_accession)%>%
  rename(gtdb_assembly="accession")

taxonomy.lloyd.healthy = rbind(taxonomy.lloyd.healthy,no.taxonomy.lloyd.healthy)

#Export top abundant flagellins with their corresponding taxonomy
write_tsv(taxonomy.lloyd.healthy,"results/taxonomy.top.metagenomes.1114.tsv")
```

# 7. Plots

#### Plot of the ranked flagellins found in healthy subjects

``` r
#This ranks the flagellins by read counts and plots the top 50 of the ranking.
pdf(file="graphics/topabundant_flagellins.pdf",height=4,width=6)
taxonomy.plot.healthy=ggplot(data=taxonomy.lloyd.healthy %>%
unite(species_accn,c(Species,Accession))%>%
  mutate(species_accn=str_remove_all(species_accn,"s__"))%>%
  mutate(Family=str_remove_all(Family,"f__"))%>%
  top_n(50), 
  aes(x=mappedReads,y=reorder(species_accn,mappedReads),color=Family))+
  geom_point(stat="identity")+
  geom_segment(aes(x=0,xend=mappedReads,yend=species_accn,y=species_accn),size=1)+
  labs(x="Total mapped reads",y="Accession")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
          axis.text.y = element_text(size = 8))+
  scale_color_npg()
```

    ## Selecting by mappedReads

``` r
  theme(legend.position = "none")#
```

    ## List of 1
    ##  $ legend.position: chr "none"
    ##  - attr(*, "class")= chr [1:2] "theme" "gg"
    ##  - attr(*, "complete")= logi FALSE
    ##  - attr(*, "validate")= logi TRUE

``` r
taxonomy.plot.healthy
```

#### Summarize genus of flagellins found in metagenomes of healthy subjects

``` r
#Plot the mapped reads for each genus in the healthy subset.Here I'm using the top 100 in the ranking.
pdf(file="graphics/genus_count_topabundant.pdf",height=4,width=6)
plot.genus.final.healthy=ggplot(data=top_n(taxonomy.lloyd.healthy,100,mappedReads)%>%count(Genus),aes(x=n,y=reorder(Genus,n)))+
  geom_bar(stat="identity",color="black", fill="#8491B4FF")+
  labs(x="Genus count",y="GTDB_Genus")+
  theme_bw()
dev.off()
```

    ## png 
    ##   2

``` r
plot.genus.final.healthy
```

![](1_HealthyHumanFlagellome_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

# 8. Final 44 candidates

``` r
#Import list of accessions of final random candidates among top abundant flagellins
accn.nonepitope=read_tsv("data/accn.nonepitope.tsv")%>%
  select(Accession,ncbi_genbank_assembly_accession.x,cds_region)
```

    ## Rows: 46 Columns: 15
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (15): Accession, ncbi_genbank_assembly_accession.x, cds_region, Genus, S...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
temp2=right_join(taxonomy.lloyd.healthy,accn.nonepitope,by="Accession")#%>%
  #select(1:15)%>%
  #rename(Genus="Genus.x")%>%
  #rename(Species="Species.x")

#Export final taxonomy table of final candidates
write_tsv(temp2,"results/taxonomy.final.abundantcandidates.tsv")
```

``` r
pdf(file="graphics/genus_count_finalcandidates.pdf",height=4,width=6)
plot.temp2=ggplot(data=temp2%>%count(Genus),aes(x=n,y=reorder(Genus,n)))+
  geom_bar(stat="identity",color="black", fill="#8491B4FF")+
  labs(x="Genus count",y="GTDB_Genus")+
  theme_bw()
dev.off()
```

    ## png 
    ##   2

``` r
plot.temp2
```

![](1_HealthyHumanFlagellome_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

``` r
sessionInfo()
```

    ## R version 4.1.2 (2021-11-01)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 18.04.5 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas/libblas.so.3
    ## LAPACK: /usr/lib/x86_64-linux-gnu/libopenblasp-r0.2.20.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] ggsci_2.9        conflicted_1.1.0 taxonomizr_0.9.3 tidytext_0.3.4  
    ##  [5] scales_1.2.0     forcats_0.5.1    stringr_1.4.0    dplyr_1.0.9     
    ##  [9] purrr_0.3.4      readr_2.1.2      tidyr_1.2.0      tibble_3.1.7    
    ## [13] ggplot2_3.3.6    tidyverse_1.3.1 
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] httr_1.4.3        bit64_4.0.5       vroom_1.5.7       jsonlite_1.8.0   
    ##  [5] modelr_0.1.8      assertthat_0.2.1  highr_0.9         blob_1.2.3       
    ##  [9] cellranger_1.1.0  yaml_2.3.5        pillar_1.7.0      RSQLite_2.2.14   
    ## [13] backports_1.4.1   lattice_0.20-45   glue_1.6.2        digest_0.6.29    
    ## [17] rvest_1.0.2       colorspace_2.0-3  htmltools_0.5.2   Matrix_1.3-4     
    ## [21] pkgconfig_2.0.3   broom_0.8.0       haven_2.5.0       tzdb_0.3.0       
    ## [25] generics_0.1.2    farver_2.1.0      ellipsis_0.3.2    cachem_1.0.6     
    ## [29] withr_2.5.0       cli_3.3.0         magrittr_2.0.3    crayon_1.5.1     
    ## [33] readxl_1.4.0      memoise_2.0.1     evaluate_0.15     tokenizers_0.2.3 
    ## [37] janeaustenr_1.0.0 fs_1.5.2          fansi_1.0.3       SnowballC_0.7.0  
    ## [41] xml2_1.3.3        tools_4.1.2       data.table_1.14.2 hms_1.1.1        
    ## [45] lifecycle_1.0.1   munsell_0.5.0     reprex_2.0.1      compiler_4.1.2   
    ## [49] rlang_1.0.3       grid_4.1.2        rstudioapi_0.13   labeling_0.4.2   
    ## [53] rmarkdown_2.14    gtable_0.3.0      DBI_1.1.3         R6_2.5.1         
    ## [57] lubridate_1.8.0   knitr_1.39        fastmap_1.1.0     bit_4.0.4        
    ## [61] utf8_1.2.2        stringi_1.7.6     parallel_4.1.2    Rcpp_1.0.8.3     
    ## [65] vctrs_0.4.1       dbplyr_2.2.1      tidyselect_1.1.2  xfun_0.31
