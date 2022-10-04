2_RhFlaB-like search
================
Andrea Borbon

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
library(gplots)
```

    ## 
    ## Attaching package: 'gplots'

    ## The following object is masked from 'package:stats':
    ## 
    ##     lowess

``` r
library(conflicted)
```

## Set working directory

``` r
knitr::opts_knit$set(root.dir = "/ebio/abt3_projects/small_projects/aborbon/TLR5/SilentFlagellins_code/")
```

##### Declare conflict preferences

``` r
conflict_prefer("filter","dplyr")
```

    ## [conflicted] Will prefer dplyr::filter over any other package

## 1. Blast results of RhFlaB-like

##### 1. DIAMOND Blastp of truncated C-terminus of RhFlaB as query against the database of full-length flagellins

``` r
blastp_flab_to_fullLength=read.delim("data/blastp_truncated_FlaB_to_fullLength_flag",header=F)
names(blastp_flab_to_fullLength)=c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
head(blastp_flab_to_fullLength)
```

    ##                            qseqid         sseqid pident length mismatch gapopen
    ## 1 FlaB_Rhominis_truncated_CD1_CD0 WP_014081191.1  100.0    117        0       0
    ## 2 FlaB_Rhominis_truncated_CD1_CD0     RGS38368.1   98.3    117        2       0
    ## 3 FlaB_Rhominis_truncated_CD1_CD0     CUN88233.1   98.3    117        2       0
    ## 4 FlaB_Rhominis_truncated_CD1_CD0 WP_055040199.1   77.8    117       26       0
    ## 5 FlaB_Rhominis_truncated_CD1_CD0     CUN28412.1   77.8    117       26       0
    ## 6 FlaB_Rhominis_truncated_CD1_CD0     RHE92765.1   77.8    117       26       0
    ##   qstart qend sstart send  evalue bitscore
    ## 1      1  117    390  506 6.9e-61    230.3
    ## 2      1  117    390  506 3.4e-60    228.0
    ## 3      1  117    390  506 7.6e-60    226.9
    ## 4      1  117    378  494 9.7e-47    183.3
    ## 5      1  117    377  493 9.7e-47    183.3
    ## 6      1  117    377  493 9.7e-47    183.3

##### Descriptive stats of blast results

``` r
#Histogram of length
ggplot(data=blastp_flab_to_fullLength,aes(x=length))+
  geom_histogram(bins=50)+
  theme_bw()+xlab("Length")+ylab("Frequency")+
  geom_vline(aes(xintercept = median(length)),col='red',size=1)+ 
  geom_vline(aes(xintercept = 117),col='blue',size=1) #Add the line indicating the length of the query (truncated CD of FlaB)
```

![](2_FlaB_like_flagellins_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
#Histogram of mismatch
ggplot(data=blastp_flab_to_fullLength,aes(x=mismatch))+
  geom_histogram(bins=50)+
  theme_bw()+xlab("Mismatch")+ylab("Frequency")+
  geom_vline(aes(xintercept = median(mismatch)),col='red',size=1) 
```

![](2_FlaB_like_flagellins_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
#Histogram of pident
ggplot(data=blastp_flab_to_fullLength,aes(x=pident))+
  geom_histogram(bins=50)+
  theme_bw()+xlab("% Identity")+ylab("Frequency")+
   geom_vline(aes(xintercept = median(pident)),col='red',size=1) 
```

![](2_FlaB_like_flagellins_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->

##### 2. Filtering results of first blast

Filtered results of blastp using the median of percentage of identity
(pident), length and mismatches. Plot the new distribution of those
parameters in the set of filtered data.

``` r
#Filter hits according to above criteria
blastp_flab_to_fullLength_filt = blastp_flab_to_fullLength%>%
  filter(length>median(length))%>%
  filter(pident>median(pident))%>%
  filter(mismatch<median(mismatch)) 

#Histogram of length
ggplot(data=blastp_flab_to_fullLength_filt,aes(x=length))+
  geom_histogram(bins=50)+
  theme_bw()+xlab("Length")+ylab("Frequency")+
  geom_vline(aes(xintercept = 117),col='blue',size=1) #Add the line indicating the length of the query (truncated CD of FlaB)
```

![](2_FlaB_like_flagellins_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
#Histogram of mismatch
ggplot(data=blastp_flab_to_fullLength_filt,aes(x=mismatch))+
  geom_histogram(bins=50)+
  theme_bw()+xlab("Mismatch")+ylab("Frequency")
```

![](2_FlaB_like_flagellins_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

``` r
#Histogram of pident
ggplot(data=blastp_flab_to_fullLength_filt,aes(x=pident))+
  geom_histogram(bins=50)+
  theme_bw()+xlab("% Identity")+ylab("Frequency")
```

![](2_FlaB_like_flagellins_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->

``` r
write_tsv(blastp_flab_to_fullLength_filt,"results/blastp_flab_to_fullLength_filt.tsv")
write_lines(blastp_flab_to_fullLength_filt$sseqid,"results/flab_tofulllength_hits_filt.tsv")
```

##### 3. Results of DIAMOND Blastp of truncated FliC against the filtered FlaB hits

``` r
#Input data: Blast result of truncated ND of FliC against filtered hits from previous step: blastp_tflic_toflab
blastp_tflic_toflab=read.delim("data/blastp_truncatedFlic_to_FlaBhits.txt",header=F)
names(blastp_tflic_toflab)=c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")

#remove duplicated hits
derep.blastp_tflic_toflab=blastp_tflic_toflab%>%
  distinct(sseqid,.keep_all=T)
  
#Histogram of length
ggplot(data=derep.blastp_tflic_toflab,aes(x=length))+
  geom_histogram(bins=30)+
  theme_bw()+xlab("Length")+ylab("Frequency")+
   geom_vline(aes(xintercept = median(length)),col='red',size=1) #median=170
```

![](2_FlaB_like_flagellins_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
#Histogram of mismatch
ggplot(data=derep.blastp_tflic_toflab,aes(x=mismatch))+
  geom_histogram(bins=30)+
  theme_bw()+xlab("Mismatch")+ylab("Frequency")+
  geom_vline(aes(xintercept = median(mismatch)),col='red',size=1) #median=83
```

![](2_FlaB_like_flagellins_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

``` r
#Histogram of pident
ggplot(data=derep.blastp_tflic_toflab,aes(x=pident))+
  geom_histogram(bins=30)+
  theme_bw()+xlab("% Identity")+ylab("Frequency")+
   geom_vline(aes(xintercept = median(pident)),col='red',size=1) #median=46.7
```

![](2_FlaB_like_flagellins_files/figure-gfm/unnamed-chunk-6-3.png)<!-- -->

##### 4. Write final sequences after de-replicating duplicated sseqids from the blast result

``` r
write_lines(derep.blastp_tflic_toflab$sseqid,"results/final_hits_flab_flic.txt")
nrow(derep.blastp_tflic_toflab)
```

    ## [1] 917

##### 5. Manual filter

In Geneious, the 917 accessions were aligned, and sorted by the residue
in the position 467. All sequences with an R (n=178) or K (n=544)
residue in that position were discarded. The list was shortened to 195
sequences

##### 6. Import accessions and summarize taxonomy

``` r
#Import GTDB taxonomy table
taxonomy.gtdb.full.db=read_tsv("data/taxonomy.gtdb.full.db.txt",col_names=T)%>%
  separate(gtdb_taxonomy,into=c("Domain","Phylum","Class","Order","Family","Genus","Species"),sep=";")
```

    ## Rows: 26392 Columns: 117
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (73): Accession, ncbi_phylum, ncbi_class, ncbi_order, ncbi_family, ncbi...
    ## dbl  (37): ncbi_taxid, ambiguous_bases, checkm_completeness, checkm_contamin...
    ## lgl   (5): gtdb_representative, gtdb_type_species_of_genus, mimag_high_quali...
    ## date  (2): ncbi_date, ncbi_seq_rel_date
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
#Import the accession numbers of the sequences kept after the manual filter
finalhits_no.rk=as_tibble(read.delim("data/3_finalhits_ids_no467RK.txt",header=F))%>%
  rename(Accession=V1)

#Taxonomic summary
taxonomy.nork.finalhits=right_join(taxonomy.gtdb.full.db,finalhits_no.rk,by="Accession")%>%
  select(contains(c("Accession","Domain","Phylum","Class","Order","Family","Genus","Species","ncbi_taxid")))%>%
  rename(gtdb_assembly_accesion="accession")%>%
  mutate_all(~replace(., is.na(.), "Unknown"))

#Plot genus in the final set
ggplot(data=taxonomy.nork.finalhits%>%count(Genus),
       aes(x=n,y=reorder(Genus,n)))+geom_bar(stat="identity")
```

![](2_FlaB_like_flagellins_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
#Plot family in the final set
ggplot(data=(group_by(taxonomy.nork.finalhits,ncbi_family)%>%
               summarise(N=n())),
       aes(x=N,y=reorder(ncbi_family,N)))+geom_bar(stat="identity")
```

![](2_FlaB_like_flagellins_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->

``` r
#Plot phylum in the final set
ggplot(data=(group_by(taxonomy.nork.finalhits,ncbi_phylum)%>%
               summarise(N=n())),
       aes(x=N,y=reorder(ncbi_phylum,N)))+geom_bar(stat="identity")
```

![](2_FlaB_like_flagellins_files/figure-gfm/unnamed-chunk-9-3.png)<!-- -->

#Dereplicate at 99% identity the final list, to reduce redundancy of
flagellins from the same species. I ran cd-hit at 99% identity on the
final list of flagellins to remove those ones that are almost identical.
In parallel, I “esearched” the accessions in the IPG to obtain their
CDSs. So, I’ll intersect the resulting table with the headers of the
dereplicated fasta. After, I can map the assembly identifiers in the CDS
file. From 131 RhFlaB-like and abundant flagellins, 125 were succesfully
associated with a CDS.

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
    ##  [1] conflicted_1.1.0 gplots_3.1.3     forcats_0.5.1    stringr_1.4.0   
    ##  [5] dplyr_1.0.9      purrr_0.3.4      readr_2.1.2      tidyr_1.2.0     
    ##  [9] tibble_3.1.7     ggplot2_3.3.6    tidyverse_1.3.1 
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] lubridate_1.8.0    gtools_3.9.3       assertthat_0.2.1   digest_0.6.29     
    ##  [5] utf8_1.2.2         R6_2.5.1           cellranger_1.1.0   backports_1.4.1   
    ##  [9] reprex_2.0.1       evaluate_0.15      highr_0.9          httr_1.4.3        
    ## [13] pillar_1.7.0       rlang_1.0.3        readxl_1.4.0       rstudioapi_0.13   
    ## [17] rmarkdown_2.14     labeling_0.4.2     bit_4.0.4          munsell_0.5.0     
    ## [21] broom_0.8.0        compiler_4.1.2     modelr_0.1.8       xfun_0.31         
    ## [25] pkgconfig_2.0.3    htmltools_0.5.2    tidyselect_1.1.2   fansi_1.0.3       
    ## [29] crayon_1.5.1       tzdb_0.3.0         dbplyr_2.2.1       withr_2.5.0       
    ## [33] bitops_1.0-7       grid_4.1.2         jsonlite_1.8.0     gtable_0.3.0      
    ## [37] lifecycle_1.0.1    DBI_1.1.3          magrittr_2.0.3     scales_1.2.0      
    ## [41] KernSmooth_2.23-20 vroom_1.5.7        cli_3.3.0          stringi_1.7.6     
    ## [45] cachem_1.0.6       farver_2.1.0       fs_1.5.2           xml2_1.3.3        
    ## [49] ellipsis_0.3.2     generics_0.1.2     vctrs_0.4.1        tools_4.1.2       
    ## [53] bit64_4.0.5        glue_1.6.2         hms_1.1.1          parallel_4.1.2    
    ## [57] fastmap_1.1.0      yaml_2.3.5         colorspace_2.0-3   caTools_1.18.2    
    ## [61] rvest_1.0.2        memoise_2.0.1      knitr_1.39         haven_2.5.0
