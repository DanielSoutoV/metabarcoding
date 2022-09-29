01 Merging mBrave data
================
Daniel
16/09/2022

I am now changing the markdown file to include my own script for BCI
metabarcoding data. This analysis relies hevily on tutorial for the
metacoder package from the Grunwald Lab.
<https://grunwaldlab.github.io/metacoder_documentation/example.html>

``` r
#the small green box to insert chunks in different formats. 
#```{r} for R.
rm(list=ls()) #reoves all objects in workspace, start with clean environment
#load packages
library('tidyr')
```

    ## Warning: package 'tidyr' was built under R version 4.1.2

``` r
library('dplyr')
```

    ## Warning: package 'dplyr' was built under R version 4.1.3

``` r
library('metacoder')
```

First thing is to merge and clean the samples from the format downloaded
from mbrave as right now its rather chaotic with column 1 sample, 2
BIN#, 3 #of seqs and 4 classificaiton.

``` r
all_merged <- read.csv('data/light_trap_19_21.csv', header=T)
head(all_merged)
```

    ##                sampleId      bin_uri sequences
    ## 1 DNA-LT-WHE2-MAR2021-B BOLD:ADU3670         2
    ## 2 DNA-LT-WHE2-MAR2021-B BOLD:AAA0085         2
    ## 3 DNA-LT-WHE2-MAR2021-B BOLD:ADF1141         2
    ## 4 DNA-LT-WHE2-MAR2021-B BOLD:AEF0071         2
    ## 5 DNA-LT-WHE2-MAR2021-B BOLD:ACK5291         2
    ## 6 DNA-LT-WHE2-MAR2021-B BOLD:ADA7150         2
    ##                                                                                                             classification
    ## 1                                                                                  Arthropoda;Insecta;Diptera;Sciaridae;;;
    ## 2                                  Arthropoda;Insecta;Lepidoptera;Choreutidae;Choreutinae;Rhobonda;Rhobonda gaurisanaDHJ04
    ## 3 Arthropoda;Insecta;Lepidoptera;Crambidae;Spilomelinae;spiloBioLep01, Hymenia;spiloBioLep01 BioLep219, Hymenia BioAlfa215
    ## 4                                     Arthropoda;Insecta;Lepidoptera;Depressariidae;Stenomatinae;Stenoma;Stenoma Janzen364
    ## 5                                                                             Arthropoda;Insecta;Coleoptera;Nitidulidae;;;
    ## 6                                                                                Arthropoda;Insecta;Coleoptera;Aderidae;;;

We need to group (and add) reads for the same BIN/classification in each
sample

``` r
all_merged %>%
  group_by(bin_uri, classification, sampleId) %>%
  summarise_all(sum) %>%
  data.frame() -> all_merged
write.csv(all_merged, 'data/dataset_merged.csv', row.names =FALSE)
```

This document now needs to be transposed (sample ID in columns, BIN_uri
and classificaiton in rows). I did in excel, surely a way to do in R.

``` r
merged_colnames <-read.csv('data/dataset_merged_samplecols.csv', header = T)
merged_colnames %>% #remove duplicated rows while adding reads to each 'sample' column
  group_by(bin_uri, classification) %>% #No more sampleID because they are the columns now
  summarise_all(sum) %>%
  data.frame() -> all_merged # so that newdf can further be used, if needed
write.csv(all_merged, 'data/merged_colnames_final.csv', row.names = F)
```

Now we have the dataset, cleaned and formatted for use with
Metacoder/MicroBiotaProcess Column 1 BIN_uri, Col2 classification (note
issues with unique names and match) cols 3-46 number of reads per
sample. The data could be further cleaned (e.g. include only
classifications beyond genus, or match to a single species name) but we
first will remove contaminants (e.g. species which ocurr with less than
X amount of reads.)

``` r
head(all_merged)
```

    ##        bin_uri
    ## 1 BOLD:AAA0009
    ## 2 BOLD:AAA0012
    ## 3 BOLD:AAA0016
    ## 4 BOLD:AAA0081
    ## 5 BOLD:AAA0082
    ## 6 BOLD:AAA0085
    ##                                                                                                                                    classification
    ## 1                               Arthropoda;Insecta;Lepidoptera;Gelechiidae;;gelJanzen1, gelJanzen01;gelJanzen01 Janzen1358, gelJanzen1 Janzen1358
    ## 2                                                                          Arthropoda;Insecta;Lepidoptera;Gelechiidae, Crambidae, Elachistidae;;;
    ## 3                                                        Arthropoda;Insecta;Lepidoptera;Depressariidae;Stenomatinae;Cerconota;Cerconota Janzen136
    ## 4                                                            Arthropoda;Insecta;Lepidoptera;Depressariidae;Stenomatinae;Stenoma;Stenoma Janzen148
    ## 5 Arthropoda;Insecta;Lepidoptera;Gelechiidae, Pterophoridae;;gelJanzen01, pterophoridJanzen01;gelJanzen01 Janzen485, pterophoridJanzen01 Janzen01
    ## 6                                                         Arthropoda;Insecta;Lepidoptera;Choreutidae;Choreutinae;Rhobonda;Rhobonda gaurisanaDHJ04
    ##   BIOLOT.39864 BIOLOT.39865 BIOLOT.39866 BIOLOT.39867 BIOLOT.39868 BIOLOT.39869
    ## 1            0            0            0            0            0            0
    ## 2            0            0            0            0            0            0
    ## 3            0            0            0            0            0            0
    ## 4            0            0            2            0            0            0
    ## 5            0            0            0            0            0            0
    ## 6            0            0            0            0            0            0
    ##   BIOLOT.39870 BIOLOT.39871 BIOLOT.39872 BIOLOT.39873 BIOLOT.39874 BIOLOT.39875
    ## 1            0            0            0            0            0            0
    ## 2            0            0            0            0            0            0
    ## 3            0            0            0            0            0            0
    ## 4            0            0            0            0            0            0
    ## 5            0            0            0            0            0            0
    ## 6            0            0            0            0            0            0
    ##   BIOLOT.39876 BIOLOT.39877 BIOLOT.39878 BIOLOT.39879 BIOLOT.39880 BIOLOT.39881
    ## 1            0            0            0            0            0            0
    ## 2            0            0            0            0            0            0
    ## 3            0            0            0            0            0            0
    ## 4            0            0            0            0            0            0
    ## 5            0            0            0            0            0            0
    ## 6            0            0            0            0            0            0
    ##   BIOLOT.39882 BIOLOT.39883 BLANKDNA_CTRL BLANKPCR_CTRL DNA.LT.ARM1.MAR2021.A
    ## 1            0            0             0             0                     0
    ## 2            0            0             0             0                     0
    ## 3            0            0             0             0                     0
    ## 4            0            0             0             0                     0
    ## 5            0            0             0             0                     0
    ## 6            0            0             0             0                     0
    ##   DNA.LT.ARM1.MAR2021.B DNA.LT.ARM2.MAR2021.A DNA.LT.ARM2.MAR2021.B
    ## 1                     0                    36                     0
    ## 2                     0                     0                     0
    ## 3                     0                     0                     0
    ## 4                     0                     0                     0
    ## 5                     0                     0                     0
    ## 6                     0                     0                     0
    ##   DNA.LT.ARM3.MAR2021.A DNA.LT.ARM3.MAR2021.B DNA.LT.ARM4.MAR2021.A
    ## 1                     0                     0                     0
    ## 2                     0                     0                     0
    ## 3                     0                     0                    62
    ## 4                     0                     3                     0
    ## 5                     0                     0                     0
    ## 6                     0                     0                     0
    ##   DNA.LT.ARM4.MAR2021.B DNA.LT.BAL1.MAR2021.A DNA.LT.BAL1.MAR2021.B
    ## 1                     0                     0                     0
    ## 2                     0                     0                     7
    ## 3                     0                     0                     0
    ## 4                     0                     0                     0
    ## 5                     0                     0                     0
    ## 6                     0                     0                     0
    ##   DNA.LT.DRA1.MAR2021.A DNA.LT.DRA1.MAR2021.B DNA.LT.WHE1.MAR2021.A
    ## 1                     0                     0                     0
    ## 2                     0                     0                     8
    ## 3                     0                     0                     0
    ## 4                     0                     0                     0
    ## 5                     0                     0                     0
    ## 6                     0                     0                     0
    ##   DNA.LT.WHE1.MAR2021.B DNA.LT.WHE2.MAR2021.A DNA.LT.WHE2.MAR2021.B
    ## 1                    40                     0                     0
    ## 2                     0                     0                     5
    ## 3                     0                     0                     0
    ## 4                     0                     0                     3
    ## 5                     0                     0                     2
    ## 6                     0                     0                     2
    ##   DNA.LT.ZET1.MAR2021.A DNA.LT.ZET1.MAR2021.B DNA.LT.ZET2.MAR2021.A
    ## 1                     0                     0                     0
    ## 2                     0                     0                     0
    ## 3                     0                     0                     0
    ## 4                     0                     0                     0
    ## 5                     9                     0                     0
    ## 6                     0                     0                     0
    ##   DNA.LT.ZET2.MAR2021.B GMP.07682_CTRL AMPtk_CTRL
    ## 1                     0              0          0
    ## 2                     0              0          0
    ## 3                     0              0          0
    ## 4                     0              0          0
    ## 5                     0              0          0
    ## 6                     0              0          0
