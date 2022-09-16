Rmarkdown
================
Daniel
16/09/2022

<<<<<<< HEAD
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
=======
## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax
for authoring HTML, PDF, and MS Word documents. For more details on
using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that
includes both content as well as the output of any embedded R code
chunks within the document. You can embed an R code chunk like this:

``` r
summary(cars)
```

    ##      speed           dist       
    ##  Min.   : 4.0   Min.   :  2.00  
    ##  1st Qu.:12.0   1st Qu.: 26.00  
    ##  Median :15.0   Median : 36.00  
    ##  Mean   :15.4   Mean   : 42.98  
    ##  3rd Qu.:19.0   3rd Qu.: 56.00  
    ##  Max.   :25.0   Max.   :120.00

## Including Plots

You can also embed plots, for example:

![](markdown_practice_files/figure-gfm/pressure-1.png)<!-- -->

Note that the `echo = FALSE` parameter was added to the code chunk to
prevent printing of the R code that generated the plot.
>>>>>>> a538aa559a469d9b789491e2ee180409f4616f62
