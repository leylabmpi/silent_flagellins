3_FlaB-like flagellins in human gut metagenome
================
Andrea Borbon

This notebook describes how we obtained the intersect between flagellins
present in human gut metagenomes and FlaB-like candidates previously
selected

#### Load libraries

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
library(tidytext)
library(conflicted)
```

``` r
knitr::opts_knit$set(root.dir = "/ebio/abt3_projects/small_projects/aborbon/TLR5/SilentFlagellins_code/")
```

``` r
conflict_prefer("filter","dplyr")
```

    ## [conflicted] Will prefer dplyr::filter over any other package

#### 1. Import files

##### 1.1 Import the list of accession IDs of flab-like flagellins with non-R/K in residue 467

``` r
#This list was produced with Geneious Pro after manually checking the flagellin alignments and selecting those accessions that did not have R/K in the position 467
tofilter=read.delim("data/4_ToFilterFlaBFliCfinalhits.txt",header=F)
names(tofilter)="Accession"
tofilter
```

    ##                        Accession
    ## 1                      Accession
    ## 2    FlaB_tr|G2T5W2|G2T5W2_ROSHA
    ## 3                     RGZ90314.1
    ## 4                     CBK92329.1
    ## 5                 WP_055224501.1
    ## 6                 WP_118392529.1
    ## 7                     RGW39946.1
    ## 8                 WP_118372266.1
    ## 9                     RGT75551.1
    ## 10                WP_117998021.1
    ## 11                WP_117719277.1
    ## 12                    RGM51729.1
    ## 13                    RHG23314.1
    ## 14                    RGN25659.1
    ## 15                    CBK91820.1
    ## 16                    CDC36194.1
    ## 17                    EFF67400.1
    ## 18                    OKZ37889.1
    ## 19                    CCY76345.1
    ## 20                WP_067199939.1
    ## 21                WP_056307372.1
    ## 22                    SEB70901.1
    ## 23                WP_042537481.1
    ## 24                    KQQ64734.1
    ## 25                    AXL12300.1
    ## 26                    KQZ25061.1
    ## 27                    PFG75253.1
    ## 28                WP_074718906.1
    ## 29                    SFG05749.1
    ## 30                    SER88836.1
    ## 31                WP_029067570.1
    ## 32                WP_014081191.1
    ## 33                    CUN88233.1
    ## 34                    RGS38368.1
    ## 35                WP_055193362.1
    ## 36                    RHM07267.1
    ## 37                    RHG26511.1
    ## 38                    CDA55210.1
    ## 39                    OLA54098.1
    ## 40                    RHC18549.1
    ## 41                    CBL13527.1
    ## 42                    RGX91734.1
    ## 43                    RHN06330.1
    ## 44                    EEV02820.1
    ## 45                WP_118209975.1
    ## 46                WP_118591627.1
    ## 47                    CDF42232.1
    ## 48                WP_004058364.1
    ## 49                    OLA45146.1
    ## 50                WP_049918078.1
    ## 51                WP_044908852.1
    ## 52                WP_024866244.1
    ## 53                WP_026650949.1
    ## 54                WP_074756729.1
    ## 55                    ADL33235.1
    ## 56                WP_027206649.1
    ## 57                    SCY34380.1
    ## 58                WP_029232837.1
    ## 59                WP_022777566.1
    ## 60                WP_026488695.1
    ## 61                    CUN28412.1
    ## 62                    RHE92765.1
    ## 63                    EEG95358.1
    ## 64                WP_118582665.1
    ## 65                WP_055040199.1
    ## 66                    RHF86582.1
    ## 67                    RGQ46029.1
    ## 68                WP_055300807.1
    ## 69                WP_075680391.1
    ## 70                    RGF42798.1
    ## 71                    RGF56743.1
    ## 72                WP_117841879.1
    ## 73                WP_118057603.1
    ## 74                    RGG47594.1
    ## 75                    CDC13483.1
    ## 76                    OLA77709.1
    ## 77                    CDA25855.1
    ## 78                    CCZ79617.1
    ## 79                    OLA60300.1
    ## 80                WP_055262733.1
    ## 81                    CRL36647.1
    ## 82                WP_117830620.1
    ## 83                WP_031542518.1
    ## 84                WP_115480276.1
    ## 85                    SCP96375.1
    ## 86                WP_092563040.1
    ## 87                WP_090016620.1
    ## 88                    KSV59214.1
    ## 89                    CDE45837.1
    ## 90                    SHL09268.1
    ## 91                    SCJ80861.1
    ## 92                WP_118374225.1
    ## 93                WP_021984889.1
    ## 94                WP_118545770.1
    ## 95                    CDE67975.1
    ## 96                    RHV03021.1
    ## 97                    PWM04712.1
    ## 98                    CDC95742.1
    ## 99                WP_118560469.1
    ## 100                   PWL71449.1
    ## 101                   CDA85633.1
    ## 102                   CCZ89920.1
    ## 103                   CDD73438.1
    ## 104                   RGF56945.1
    ## 105               WP_092454342.1
    ## 106                   CCZ43036.1
    ## 107               WP_055944303.1
    ## 108                   OKZ81168.1
    ## 109               WP_118744492.1
    ## 110               WP_027105261.1
    ## 111                   RHP21408.1
    ## 112               WP_058257058.1
    ## 113                   ADZ82773.1
    ## 114               WP_105618991.1
    ## 115                   CRZ35733.1
    ## 116               WP_054742738.1
    ## 117               WP_073290153.1
    ## 118               WP_113672800.1
    ## 119               WP_114375048.1
    ## 120                   OGH96420.1
    ## 121               WP_077529038.1
    ## 122                   PAU80189.1
    ## 123 AHA06007.1_FliC_STyphimurium
    ## 124               WP_092872173.1
    ## 125                   RGS33867.1
    ## 126               WP_012740256.1
    ## 127                   RHV15618.1
    ## 128                   CCZ53492.1
    ## 129                   RHC12202.1
    ## 130               WP_118344389.1
    ## 131                   RHI64475.1
    ## 132                   CDA38946.1
    ## 133                   RHM11737.1
    ## 134                   RGW88779.1
    ## 135                   RHK43607.1
    ## 136                   SDN11166.1
    ## 137               WP_027431085.1
    ## 138                   CDD55527.1
    ## 139               WP_008116478.1
    ## 140               WP_027437934.1
    ## 141                   CDB65654.1
    ## 142                   CUP37306.1
    ## 143                   OLA13988.1
    ## 144                   CDB68961.1
    ## 145               WP_027428254.1
    ## 146               WP_018921994.1
    ## 147                   SMF32740.1
    ## 148               WP_042214857.1
    ## 149                   CDE55985.1
    ## 150                   RHT20893.1
    ## 151                   RHU62160.1
    ## 152                   RHP09176.1
    ## 153                   CDC99374.1
    ## 154                   OLA53353.1
    ## 155                   CDA99895.1
    ## 156                   SCH11803.1
    ## 157                   RHP98655.1
    ## 158                   RHU92598.1
    ## 159               WP_118551969.1
    ## 160                   RHV48110.1
    ## 161                   SCH62551.1
    ## 162               WP_027435651.1
    ## 163                   SOY27460.1
    ## 164               WP_022777320.1
    ## 165               WP_026517790.1
    ## 166               WP_022783963.1
    ## 167                   CDF07358.1
    ## 168               WP_074462912.1
    ## 169                   ADL32814.1
    ## 170               WP_026663022.1
    ## 171                   CDD48964.1
    ## 172               WP_089839980.1
    ## 173               WP_022772103.1
    ## 174               WP_026656233.1
    ## 175               WP_092336810.1
    ## 176               WP_092329363.1
    ## 177                   SEG35387.1
    ## 178               WP_092244991.1
    ## 179               WP_029231470.1
    ## 180               WP_022766863.1
    ## 181               WP_034445183.1
    ## 182               WP_024866128.1
    ## 183               WP_026513482.1
    ## 184               WP_026521459.1
    ## 185               WP_016293258.1
    ## 186               WP_031391632.1
    ## 187               WP_026489210.1
    ## 188               WP_026519205.1
    ## 189               WP_022769365.1
    ## 190                   SCX92535.1
    ## 191               WP_027207228.1
    ## 192               WP_027216912.1
    ## 193               WP_022758790.1
    ## 194               WP_027205020.1
    ## 195                   SEQ53841.1
    ## 196                   SES34383.1
    ## 197               WP_016300279.1
    ## 198                   SDL67829.1

##### 1.2. Import Franzosa mapping results

``` r
flagellin_diamond_franzosa=read.delim("data/mergedFlagFranzosa.txt",header=F,sep="")
names(flagellin_diamond_franzosa)=c("Study","Sample","qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
```

##### 1.3. Import Lloyd-Price mapping results

``` r
#Import diamond blastx results - reads mapped to flagellin accessions  
flagellin_diamond_lloyd=read.delim("data/mergedFlagellinLloyd.txt",header=F,sep="")%>%
  mutate(Study="LloydPrice_2019")
names(flagellin_diamond_lloyd)=c("Sample","qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","Study")
```

##### 1.4. Merge LloydPrice and Franzosa results

``` r
mg.mapping.results=rbind(flagellin_diamond_franzosa,flagellin_diamond_lloyd)%>%
  rename(Accession=sseqid)

mg.mapping.results=mg.mapping.results%>%
  group_by(Accession)%>%
  count()

mg.mapping.results %>%
  head()
```

    ## # A tibble: 6 × 2
    ## # Groups:   Accession [6]
    ##   Accession      n
    ##   <chr>      <int>
    ## 1 AAA17857.1     2
    ## 2 AAA17858.1     1
    ## 3 AAA17865.1     1
    ## 4 AAA23797.1    83
    ## 5 AAA23798.1   165
    ## 6 AAA23799.1    30

#### 2. Intersect of FlaB-like in human metagenomes

Here all flagellins present in human metagenomes will be intersected
with the set of the FlaB-like flagellins from 1.1

``` r
intersect.flab.mgs=inner_join(tofilter,mg.mapping.results)
```

    ## Joining, by = "Accession"

``` r
nrow(intersect.flab.mgs)
```

    ## [1] 145

``` r
intersect.flab.mgs %>%
  head()
```

    ##        Accession     n
    ## 1     RGZ90314.1 58714
    ## 2     CBK92329.1 21274
    ## 3 WP_055224501.1  4912
    ## 4 WP_118392529.1  3916
    ## 5     RGW39946.1 24402
    ## 6 WP_118372266.1  4240

#### 3. Export accessions of FlaB-like flagellins present in human metagenomes

This list will we used to pull out FlaB-like flagellins present in human
metagenomes from the comlpete flagellin database

``` r
write_lines(intersect.flab.mgs$Accession,"data/intersect.flab.metagones_accn.tsv")
```

#### 4. Remove FlgL sequences and de-replicate final sequences

``` r
#subset intersect.flab.lloyd from flagellinDB
# seqtk subseq /ebio/abt3_projects/small_projects/aborbon/TLR5/Dalong_Flagellin/flic_flab_search/Dalong_flagellinDB_100 /ebio/abt3_projects/small_projects/aborbon/TLR5/Dalong_Flagellin/flic_flab_search/diamond/1_flab_2_flic/intersect.flab.metagones_accn.tsv > intersect.flab.metagones_accn.fasta
# 
# #Manually removed all proteins annotated as "FlgL"
# grep "FlgL" intersect.flab.metagones_accn.fasta | cat > toRemove.tsv
# seqtk subseq intersect.flab.metagones_accn.fasta toRemove.tsv > finalhits_ibd_tlr5epitope_flab_flic.fasta #Revisar aquí como uso seqtk para esto
# #This reduced the dataset to 129 sequences
```

#### 5. Dereplication of final sequences

Dereplication using CD-HIT at 0.99 identity.

``` r
#Dereplication CD-HIT
# cd-hit -i finalhits_ibd_tlr5epitope_flab_flic.fasta -o derep0.99.flablike.humanmetagenomes.fasta -c 0.99 -T 10 -M 16000 -n 5
```

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
    ##  [1] conflicted_1.1.0 tidytext_0.3.4   forcats_0.5.1    stringr_1.4.0   
    ##  [5] dplyr_1.0.9      purrr_0.3.4      readr_2.1.2      tidyr_1.2.0     
    ##  [9] tibble_3.1.7     ggplot2_3.3.6    tidyverse_1.3.1 
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.8.3      lubridate_1.8.0   lattice_0.20-45   assertthat_0.2.1 
    ##  [5] digest_0.6.29     utf8_1.2.2        R6_2.5.1          cellranger_1.1.0 
    ##  [9] backports_1.4.1   reprex_2.0.1      evaluate_0.15     httr_1.4.3       
    ## [13] pillar_1.7.0      rlang_1.0.3       readxl_1.4.0      rstudioapi_0.13  
    ## [17] Matrix_1.3-4      rmarkdown_2.14    bit_4.0.4         munsell_0.5.0    
    ## [21] broom_0.8.0       compiler_4.1.2    janeaustenr_1.0.0 modelr_0.1.8     
    ## [25] xfun_0.31         pkgconfig_2.0.3   htmltools_0.5.2   tidyselect_1.1.2 
    ## [29] fansi_1.0.3       crayon_1.5.1      tzdb_0.3.0        dbplyr_2.2.1     
    ## [33] withr_2.5.0       SnowballC_0.7.0   grid_4.1.2        jsonlite_1.8.0   
    ## [37] gtable_0.3.0      lifecycle_1.0.1   DBI_1.1.3         magrittr_2.0.3   
    ## [41] scales_1.2.0      tokenizers_0.2.3  vroom_1.5.7       cachem_1.0.6     
    ## [45] cli_3.3.0         stringi_1.7.6     fs_1.5.2          xml2_1.3.3       
    ## [49] ellipsis_0.3.2    generics_0.1.2    vctrs_0.4.1       tools_4.1.2      
    ## [53] bit64_4.0.5       glue_1.6.2        hms_1.1.1         parallel_4.1.2   
    ## [57] fastmap_1.1.0     yaml_2.3.5        colorspace_2.0-3  rvest_1.0.2      
    ## [61] memoise_2.0.1     knitr_1.39        haven_2.5.0
