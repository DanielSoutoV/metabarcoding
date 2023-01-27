01 Merging all required data
================
Daniel
16/09/2022

There is a lot of manipulation that needs to be done with the data. In
general, what we need to analyze metabarcoding data is a data frame
organized as: Bold Bin (column A) and classification (separated by
taxonomic rank by “;” as: k\_\_ ; p\_\_ ; etc. in column B) in rows, and
samples in columns (c - onwards). Some of the main issues with the data
is that we have multiple duplicated BINs as they are present in multiple
sample columns. The scripts below use tidyr and dplyr to merge these
duplicates and sum the number of reads in the columns.

The scripts create multiple files which will be used in the subsequent
scripts. Sometimes they are redundant because I had to manually check
these files for more duplicates. We have issues with BOLD bins being the
same even though classified as different species, and several cases
where we have duplicated species with different BINs.

There is one big problem with nomenclature. Particularly when comparing
traditional taxonomic records vs metabarcoded records. The problem is
that duplicate BOLD BINs are named different in the BOLD database (the
name used by mBrave) and the ForestGEO (the names used for the
traditional data set). The names can be completely different, but
sometimes include spaces or periods, such as sp. 1YB vs sp1YB or sp. 1
YB.

Be aware that this is a very Frankenstein-ish pipeline which heavily
depends on other tutorials I’ve found. Specifically the metacoder
package from the Grunwald Lab
(<https://grunwaldlab.github.io/metacoder_documentation/example.html>)
and the Microbiotaprocess package from the Yu Lab
(<https://github.com/YuLab-SMU/MicrobiotaProcess>). Some of the scripts
are my own (the dirtier ones), and some of the manipulation has been
done in Excel (since I don’t know the best way to do so in R -
particularly the checking for duplicates and editing the species names).

############ Warning: These scripts have been ran already, and are meant only for information on the data cleaning/merging process. DO NOT RUN THEM as you risk overwriting the final datasets (it does not really matter since you will arrive to the same final data set but it will save you time).

Any questions or comments can be asked directly on the GitHub repo or to
my email: <daniel.souto.v@gmail.com>

First thing is to merge and clean the samples from the format downloaded
from mbrave as right now its rather chaotic with column 1 sample, 2
BIN#, 3 #of READS and 4 classification. Note that I am using the
formatted spreadsheet as directly downloaded from mBrave. When in doubt
of what I mean, look at the original file called (light_trap_19_21.csv),
and then the written new file (dataset_merged.csv).

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
    ##                                                                                                                              classification
    ## 1                                                                      k__Animalia;p__Arthropoda;c__Insecta;o__Diptera;f__Sciaridae;g__;s__
    ## 2                                 k__Animalia;p__Arthropoda;c__Insecta;o__Lepidoptera;f__Choreutidae;g__Rhobonda;s__Rhobonda gaurisanaDHJ04
    ## 3 k__Animalia;p__Arthropoda;c__Insecta;o__Lepidoptera;f__Crambidae;g__spiloBioLep01, Hymenia;s__spiloBioLep01 BioLep219, Hymenia BioAlfa215
    ## 4                                     k__Animalia;p__Arthropoda;c__Insecta;o__Lepidoptera;f__Depressariidae;g__Stenoma;s__Stenoma Janzen364
    ## 5                                                                 k__Animalia;p__Arthropoda;c__Insecta;o__Coleoptera;f__Nitidulidae;g__;s__
    ## 6                                                                    k__Animalia;p__Arthropoda;c__Insecta;o__Coleoptera;f__Aderidae;g__;s__

``` r
all_merged %>%
  group_by(bin_uri, classification, sampleId) %>%
  summarise_all(sum) %>%
  data.frame() -> all_merged
write.csv(all_merged, 'data/dataset_merged.csv', row.names =FALSE)
```

This document now needs to be transposed (sample ID in columns, BIN_uri
and classification in rows). I did this in excel, surely a way to do in
R with dplyr probably with the “transpose_long” function. Then once
again we group the bins and classifications and add the columns to avoid
repeating BINs.

``` r
merged_colnames <-read.csv('data/dataset_merged_samplecols.csv', header = T) #this is a modified version of dataset_merged.csv above but now samples are in columns.
merged_colnames %>% #remove duplicated rows while adding reads to each 'sample' column
  group_by(bin_uri, classification) %>% #No more sampleID because they are the columns now
  summarise_all(sum) %>%
  data.frame() -> all_merged # so that newdf can further be used, if needed
write.csv(all_merged, 'data/merged_colnames_final.csv', row.names = F) ##Be aware of column names - check that max row numbers is >0, delete those which are 0
```

Now we have the dataset, cleaned and formatted for use with
Metacoder/MicroBiotaProcess Column 1 BIN_uri (this column is later
renamed to bold_bin), Col2: classification (note issues with unique
names and match) cols 3-46 number of reads per sample.

NOTE: The file used for metabarcoding data beyond this point is
######finaldat_merged.csv######### - this file is based off the
merged_colnames_final.csv but I went one by one on the classification
column changing the names to have a single species name (look at both
files if in doubt). This was necessary because some of BOLD BINs have
several species names. As a rule of thumb, I picked the first name of
the list which corresponds (generally) to the name of the DS-BCIARTH,
the database within BOLD supplied by ForestGEO Arthrpod Initiative which
corresponds to Barro Colorado Island, where this data comes from.

The following chunks follow a similar process but now to merge and clean
the traditional data, and add the metabarcoidng data.

``` r
traditional <-read.csv('data/traditionaldatabinary.csv') #this .csv was created by taking all ForestGEO recods which have a BIN.
traditional %>%
  group_by(bold_bin, classification) %>%
  summarise_all(sum) %>%
  data.frame() -> traditional
write.csv(traditional, 'data/tradbinary_merged.csv', row.names = FALSE) #This new file now has all BINs only once, for each 'sample' in which it ocurrs.

#########
traditional2 <- read.csv('data/tradbinary_merged.csv')
traditional2 %>%
  group_by(bold_bin, classification) %>%
  summarise_all(sum) %>%
  data.frame() -> traditional2
write.csv(traditional2, 'data/tradbinary_merged.csv', row.names = FALSE)

#This second run was to make sure that no BOLD_BINs were repeated, I manually checked the first tradbinary_merged.csv file and changed names which needed to be double checked. I did this by hihglighting duplicate BOLD_BINs in excel and check any which was duplicated whit a different name.
########################################
```

The next chunk is run to include the metabarcoding data to the newly
created ‘tradbinary_merged.csv’ To do this, there is probably a
right_join function to run in R, but instead I did so manually in excel,
adding all metabarcoding BINs to the final row of the spreadsheet plus
all the columns from the mtabarcoding to the final column of the
tradbinary_merged.csv. (look at trad_and_meta_toMerge.csv)

``` r
tradmeta_toMerge <- read.csv('data/trad_and_meta_toMerge.csv') #this new file is the tradbinary_merged.csv with the metabarcoding data added to the end.
tradmeta_toMerge %>%
  group_by(bold_bin, classification) %>%
  summarise_all(sum) %>%
  data.frame() -> tradmeta_toMerge
write.csv(tradmeta_toMerge, 'data/trad_meta_merged1.csv', row.names = FALSE)

trad_meta01 <- read.csv('data/trad_meta_merged1.csv')
trad_meta01 %>%
  group_by(bold_bin, classification) %>%
  summarise_all(sum) %>%
  data.frame() -> trad_meta01
write.csv(trad_meta01, 'data/final_tradmeta_merged.csv', row.names = FALSE)

tradmetabr <- read.csv('data/final_tradmeta_merged.csv')
sample <- read.csv('data/metadata_trad_metabr.csv') #this includes info on location/collection method/classification method and any more info can be added
head(tradmetabr) #Data already organized by samples as columns but need to merge repeated BINs
```

    ##       bold_bin
    ## 1 BOLD:AAA0009
    ## 2 BOLD:AAA0012
    ## 3 BOLD:AAA0016
    ## 4 BOLD:AAA0081
    ## 5 BOLD:AAA0082
    ## 6 BOLD:AAA0085
    ##                                                                                                    classification
    ## 1              k__Animalia;p__Arthropoda;c__Insecta;o__Lepidoptera;f__Gelechiidae;s__;g__gelJanzen1;s__Janzen1358
    ## 2                                  k__Animalia;p__Arthropoda;c__Insecta;o__Lepidoptera;f__Gelechiidae;s__;g__;s__
    ## 3 k__Animalia;p__Arthropoda;c__Insecta;o__Lepidoptera;f__Depressariidae;s__Stenomatinae;g__Cerconota;s__Janzen136
    ## 4   k__Animalia;p__Arthropoda;c__Insecta;o__Lepidoptera;f__Depressariidae;s__Stenomatinae;g__Stenoma;s__Janzen148
    ## 5              k__Animalia;p__Arthropoda;c__Insecta;o__Lepidoptera;f__Gelechiidae;s__;g__gelJanzen01;s__Janzen485
    ## 6 k__Animalia;p__Arthropoda;c__Insecta;o__Lepidoptera;f__Choreutidae;s__Choreutinae;g__Rhobonda;s__gaurisanaDHJ04
    ##   ALTOSCAMP_na_LT ARGOS_dry_BA ARGOS_wet_BA ARM1_dry_B ARM1_wet_B B05_wet_BB
    ## 1               0            0            0          0          0          0
    ## 2               0            0            0          0          0          0
    ## 3               0            0            0          0          0          0
    ## 4               0            0            0          0          0          0
    ## 5               0            0            0          0          0          0
    ## 6               0            0            0          0          0          0
    ##   ARM1_dry_BB ARM1_wet_BB ARM1_dry_BT ARM1_wet_BT ARM1_wet_LT ARM1_dry_LT
    ## 1           0           0           0           0           0           0
    ## 2           0           0           0           0           0           0
    ## 3           0           0           0           0           0           0
    ## 4           0           0           0           0           0           0
    ## 5           0           0           0           0           0           0
    ## 6           0           0           0           0           0           0
    ##   ARM1_dry_W ARM2_dry_B ARM2_wet_B TEST_dry_BB ARM2_dry_BB B06_wet_BB
    ## 1          0          0          0           0           0          0
    ## 2          0          0          0           0           0          0
    ## 3          0          0          0           0           0          0
    ## 4          0          0          0           0           0          0
    ## 5          0          0          0           0           0          0
    ## 6          0          0          0           0           0          0
    ##   ARM2_wet_BB ARM2_wet_BT ARM2_dry_BT ARM2_wet_LT ARM2_dry_LT ARM2_dry_W
    ## 1           0           0           0           0           0          0
    ## 2           0           0           0           0           0          0
    ## 3           0           0           0           0           0          0
    ## 4           0           0           0           0           0          0
    ## 5           0           0           0           0           0          0
    ## 6           0           0           0           0           0          0
    ##   ARM3_dry_B ARM3_wet_B ARM3_wet_BB ARM3_dry_BB ARM3_dry_BT ARM3_wet_BT
    ## 1          0          0           0           0           0           0
    ## 2          0          0           0           0           0           0
    ## 3          0          0           0           0           0           0
    ## 4          0          0           0           0           0           0
    ## 5          0          0           0           0           0           0
    ## 6          0          0           0           0           0           0
    ##   ARM3_wet_LT ARM3_dry_LT L02_wet_LT L02TEST_dry_LT L02_dry_LT M07_dry_MP
    ## 1           0           0          0              0          0          0
    ## 2           0           0          0              0          0          0
    ## 3           0           0          0              0          0          0
    ## 4           0           0          0              0          0          0
    ## 5           0           0          0              0          0          0
    ## 6           0           0          0              0          0          0
    ##   ARM3_dry_W ARM4_wet_B ARM4_dry_B B08_wet_BB ARM4_dry_BB ARM4_wet_BB
    ## 1          0          0          0          0           0           0
    ## 2          0          0          0          0           0           0
    ## 3          0          0          0          0           0           0
    ## 4          0          0          0          0           0           0
    ## 5          0          0          0          0           0           0
    ## 6          0          0          0          0           0           0
    ##   ARM4_dry_BT ARM4_wet_BT ARM4_wet_LT ARM4_dry_LT M08_dry_MP ARM4_dry_W
    ## 1           0           0           0           0          0          0
    ## 2           0           0           0           0          0          0
    ## 3           0           0           0           0          0          0
    ## 4           0           0           0           0          0          0
    ## 5           0           0           0           0          0          0
    ## 6           0           0           0           0          0          0
    ##   ATO_LESM_na_B B_na_BM B_na_TBA BAL1_dry_B BAL1_wet_B BAL1_dry_BB BAL1_wet_BB
    ## 1             0       0        0          0          0           0           0
    ## 2             0       0        0          0          0           0           0
    ## 3             0       0        0          0          0           0           0
    ## 4             0       0        0          0          0           0           0
    ## 5             0       0        0          0          0           0           0
    ## 6             0       0        0          0          0           0           0
    ##   BAL1_wet_BT BAL1_dry_BT BAL1_dry_LT BAL1_wet_LT M01_dry_MP BAL1_dry_W
    ## 1           0           0           0           0          0          0
    ## 2           0           0           0           0          0          0
    ## 3           0           0           0           0          0          0
    ## 4           0           0           0           0          0          0
    ## 5           0           0           0           0          0          0
    ## 6           0           0           0           0          0          0
    ##   ARBOREALBUJAN_na_AM ASPAN_na_B BCI_dry_BE BCI_wet_BM BCI_dry_BM BCI_na_HC
    ## 1                   0          0          0          0          0         0
    ## 2                   0          0          0          0          0         0
    ## 3                   0          0          0          0          0         0
    ## 4                   0          0          0          0          0         0
    ## 5                   0          0          0          0          0         0
    ## 6                   0          0          0          0          0         0
    ##   BCI_wet_HC BCI_dry_HC ALATESHIK_na_LT MANUALJESSE_na_MA ARBOREALENDARA_na_MA
    ## 1          0          0               0                 0                    0
    ## 2          0          0               0                 0                    0
    ## 3          0          0               0                 0                    0
    ## 4          0          0               0                 0                    0
    ## 5          0          0               0                 0                    0
    ## 6          0          0               0                 0                    0
    ##   MANUALDONOSOLIT_na_MA MANUALDONOSOUND_na_MA BCI_wet_MM BCI_dry_MM
    ## 1                     0                     0          0          0
    ## 2                     0                     0          0          0
    ## 3                     0                     0          0          0
    ## 4                     0                     0          0          0
    ## 5                     0                     0          0          0
    ## 6                     0                     0          0          0
    ##   FARAOC_dry_RP FARAOC_na_RP FROGDIETDONOSO_na_SC TR8_dry_T TR10_dry_T
    ## 1             0            0                    0         0          0
    ## 2             0            0                    0         0          0
    ## 3             0            0                    0         0          0
    ## 4             0            0                    0         0          0
    ## 5             0            0                    0         0          0
    ## 6             0            0                    0         0          0
    ##   TR11_dry_T TR2_dry_T TR7_dry_T TR9_dry_T TR1_dry_T TR6_dry_T BL_dry_BM
    ## 1          0         0         0         0         0         0         0
    ## 2          0         0         0         0         0         0         0
    ## 3          0         0         0         0         0         0         0
    ## 4          0         0         0         0         0         0         0
    ## 5          0         0         0         0         0         0         0
    ## 6          0         0         0         0         0         0         0
    ##   BL_dry_TH C37_dry_BA C37_wet_BA C42_wet_BA C42_dry_BA C48_dry_BA C48_wet_BA
    ## 1         0          0          0          0          0          0          0
    ## 2         0          0          0          0          0          0          0
    ## 3         0          0          0          0          0          0          0
    ## 4         0          0          0          0          0          0          0
    ## 5         0          0          0          0          0          0          0
    ## 6         0          0          0          0          0          0          0
    ##   CED_dry_T C57_wet_BA C57_dry_BA C64_wet_BA C64_dry_BA CAC_wet_HC
    ## 1         0          0          0          0          0          0
    ## 2         0          0          0          0          0          0
    ## 3         0          0          0          0          0          0
    ## 4         0          0          0          0          0          0
    ## 5         0          0          0          0          0          0
    ## 6         0          0          0          0          0          0
    ##   ALDERCOTTE_na_LA CHAC_LESM_na_B CHAN_LESM_na_B COM_LESM_na_B CUED_LESM_na_B
    ## 1                0              0              0             0              0
    ## 2                0              0              0             0              0
    ## 3                0              0              0             0              0
    ## 4                0              0              0             0              0
    ## 5                0              0              0             0              0
    ## 6                0              0              0             0              0
    ##   D_na_BM D_na_TBA DRA1_dry_B DRA1_wet_B DRA1_dry_BB DRA1_wet_BB DRA1_wet_BT
    ## 1       0        0          0          0           0           0           0
    ## 2       0        0          0          0           0           0           0
    ## 3       0        0          0          0           0           0           0
    ## 4       0        0          0          0           0           0           0
    ## 5       0        0          0          0           0           0           0
    ## 6       0        0          0          0           0           0           0
    ##   DRA1_dry_BT DRA1_dry_LT DRA1_wet_LT L03_wet_LT L03_dry_LT DRA1_dry_W
    ## 1           0           0           0          0          0          0
    ## 2           0           0           0          0          0          0
    ## 3           0           0           0          0          0          0
    ## 4           0           0           0          0          0          0
    ## 5           0           0           0          0          0          0
    ## 6           0           0           0          0          0          0
    ##   EMI_LESM_na_B FB_dry_BM FB_dry_TBA FM1_dry_BM FM2_dry_BM GIG_wet_TN
    ## 1             0         0          0          0          0          0
    ## 2             0         0          0          0          0          0
    ## 3             0         0          0          0          0          0
    ## 4             0         0          0          0          0          0
    ## 5             0         0          0          0          0          0
    ## 6             0         0          0          0          0          0
    ##   HER_LESM_na_B IXT_LESM_na_B JDH_dry_BM JVT_dry_BM JVT_dry_TBA JZ_dry_TBA
    ## 1             0             0          0          0           0          0
    ## 2             0             0          0          0           0          0
    ## 3             0             0          0          0           0          0
    ## 4             0             0          0          0           0          0
    ## 5             0             0          0          0           0          0
    ## 6             0             0          0          0           0          0
    ##   L_dry_BM NRG_LESM_na_B OPEN_na_HC PAN18FG062_na_HC PAN18FG009_na_HC
    ## 1        0             0          0                0                0
    ## 2        0             0          0                0                0
    ## 3        0             0          0                0                0
    ## 4        0             0          0                0                0
    ## 5        0             0          0                0                0
    ## 6        0             0          0                0                0
    ##   PAN18FG078_na_HC PAN18FG060_na_HC PAN18FG024_na_HC PAN18FG004_na_HC
    ## 1                0                0                0                0
    ## 2                0                0                0                0
    ## 3                0                0                0                0
    ## 4                0                0                0                0
    ## 5                0                0                0                0
    ## 6                0                0                0                0
    ##   PAN18FG080_na_HC PNM_wet_HC PRD_dry_HC PRD_wet_HC PRD_na_HC PRD_dry_LT
    ## 1                0          0          0          0         0          0
    ## 2                0          0          0          0         0          0
    ## 3                0          0          0          0         0          0
    ## 4                0          0          0          0         0          0
    ## 5                0          0          0          0         0          0
    ## 6                0          0          0          0         0          0
    ##   PUC_LESM_na_B RIOFRIO2_na_HC RIOFRIO1_na_MM SOL_dry_BM T38_wet_BA T38_dry_BA
    ## 1             0              0              0          0          0          0
    ## 2             0              0              0          0          0          0
    ## 3             0              0              0          0          0          0
    ## 4             0              0              0          0          0          0
    ## 5             0              0              0          0          0          0
    ## 6             0              0              0          0          0          0
    ##   T39_wet_BA T39_dry_BA TER_dry_T T46_dry_BA T46_wet_BA T47_wet_BA T47_dry_BA
    ## 1          0          0         0          0          0          0          0
    ## 2          0          0         0          0          0          0          0
    ## 3          0          0         0          0          0          0          0
    ## 4          0          0         0          0          0          0          0
    ## 5          0          0         0          0          0          0          0
    ## 6          0          0         0          0          0          0          0
    ##   T61_wet_BA T61_dry_BA TAP_LESM_na_B TAPIJ_LESM_na_B TB_dry_TBA TB1_dry_BM
    ## 1          0          0             0               0          0          0
    ## 2          0          0             0               0          0          0
    ## 3          0          0             0               0          0          0
    ## 4          0          0             0               0          0          0
    ## 5          0          0             0               0          0          0
    ## 6          0          0             0               0          0          0
    ##   TB2_dry_BM TB2_dry_TBA TEX_LESM_na_B TGP_dry_TBA TEK_dry_T TK2_wet_BA
    ## 1          0           0             0           0         0          0
    ## 2          0           0             0           0         0          0
    ## 3          0           0             0           0         0          0
    ## 4          0           0             0           0         0          0
    ## 5          0           0             0           0         0          0
    ## 6          0           0             0           0         0          0
    ##   TK3_wet_BA TK4_wet_BA TK5_wet_BA TUP_dry_HC TUX_LESM_na_B WHE1_wet_B
    ## 1          0          0          0          0             0          0
    ## 2          0          0          0          0             0          0
    ## 3          0          0          0          0             0          0
    ## 4          0          0          0          0             0          0
    ## 5          0          0          0          0             0          0
    ## 6          0          0          0          0             0          0
    ##   WHE1_dry_B WHE1_dry_BB WHE1_wet_BB B02_wet_BB WHE1_dry_BT WHE1_wet_BT
    ## 1          0           0           0          0           0           0
    ## 2          0           0           0          0           0           0
    ## 3          0           0           0          0           0           0
    ## 4          0           0           0          0           0           0
    ## 5          0           0           0          0           0           0
    ## 6          0           0           0          0           0           0
    ##   WHE1_wet_LT WHE1_dry_LT M03_dry_MP WHE1_dry_W WHE2_dry_B WHE2_wet_B
    ## 1           0           0          0          0          0          0
    ## 2           0           0          0          0          0          0
    ## 3           0           0          0          0          0          0
    ## 4           0           0          0          0          0          0
    ## 5           0           0          0          0          0          0
    ## 6           0           0          0          0          0          0
    ##   WHE2_dry_BB WHE2_dry_BT WHE2_wet_BT L05_wet_LT L05TEST_wet_LT WHE2_dry_LT
    ## 1           0           0           0          0              0           0
    ## 2           0           0           0          0              0           0
    ## 3           0           0           0          0              0           0
    ## 4           0           0           0          0              0           0
    ## 5           0           0           0          0              0           0
    ## 6           0           0           0          0              0           0
    ##   WHE2_wet_LT L05TEST_dry_LT WHE2_dry_W WMW_dry_BM WMW_dry_TBA ZET1_wet_B
    ## 1           0              0          0          0           0          0
    ## 2           0              0          0          0           0          0
    ## 3           0              0          0          0           0          0
    ## 4           0              0          0          0           0          0
    ## 5           0              0          0          0           0          0
    ## 6           0              0          0          0           0          0
    ##   ZET1_dry_B ZET1_dry_BB B09_wet_BB ZET1_wet_BB ZET1_wet_BT ZET1_dry_BT
    ## 1          0           0          0           0           0           0
    ## 2          0           0          0           0           0           0
    ## 3          0           0          0           0           0           0
    ## 4          0           0          0           0           0           0
    ## 5          0           0          0           0           0           0
    ## 6          0           0          0           0           0           0
    ##   L04_wet_LT ZET1_wet_LT ZET1_dry_LT L04TEST_dry_LT ZET1_dry_W ZET2_wet_B
    ## 1          0           0           0              0          0          0
    ## 2          0           0           0              0          0          0
    ## 3          0           0           0              0          0          0
    ## 4          0           0           0              0          0          0
    ## 5          0           0           0              0          0          0
    ## 6          0           0           0              0          0          0
    ##   ZET2_dry_B ZET2_dry_BB B10_wet_BB ZET2_wet_BT ZET2_dry_BT ZET2_wet_LT
    ## 1          0           0          0           0           0           0
    ## 2          0           0          0           0           0           0
    ## 3          0           0          0           0           0           0
    ## 4          0           0          0           0           0           0
    ## 5          0           0          0           0           0           0
    ## 6          0           0          0           0           0           0
    ##   ZET2_dry_LT L01_wet_LT L01TEST_dry_LT L01_dry_LT M10_dry_MP ZET2_dry_W
    ## 1           0          0              0          0          0          0
    ## 2           0          0              0          0          0          0
    ## 3           0          0              0          0          0          0
    ## 4           0          0              0          0          0          0
    ## 5           0          0              0          0          0          0
    ## 6           0          0              0          0          0          0
    ##   WINKLERDONOSO_na_W total39864 total39865 total39866 total39867 total39868
    ## 1                  0          0          0          0          0          0
    ## 2                  0          0          0          0          0          0
    ## 3                  0          0          0          0          0          0
    ## 4                  0          0          0          1          0          0
    ## 5                  0          0          0          0          0          0
    ## 6                  0          0          0          0          0          0
    ##   total39869 total39870 total39871 total39872 total39873 total39874 total39875
    ## 1          0          0          0          0          0          0          0
    ## 2          0          0          0          0          0          0          0
    ## 3          0          0          0          0          0          0          0
    ## 4          0          0          0          0          0          0          0
    ## 5          0          0          0          0          0          0          0
    ## 6          0          0          0          0          0          0          0
    ##   total39876 total39877 total39878 total39879 total39880 total39881 total39882
    ## 1          0          0          0          0          0          0          0
    ## 2          0          0          0          0          0          0          0
    ## 3          0          0          0          0          0          0          0
    ## 4          0          0          0          0          0          0          0
    ## 5          0          0          0          0          0          0          0
    ## 6          0          0          0          0          0          0          0
    ##   total39883 ARM1A ARM1B ARM2A ARM2B ARM3A ARM3B ARM4A ARM4B BAL1A BAL1B DRA1A
    ## 1          0     0     0     1     0     0     0     0     0     0     0     0
    ## 2          0     0     0     0     0     0     0     0     0     0     1     0
    ## 3          0     0     0     0     0     0     0     1     0     0     0     0
    ## 4          0     0     0     0     0     0     1     0     0     0     0     0
    ## 5          0     0     0     0     0     0     0     0     0     0     0     0
    ## 6          0     0     0     0     0     0     0     0     0     0     0     0
    ##   DRA1B WHE1A WHE1B WHE2A WHE2B ZET1A ZET1B ZET2A ZET2B
    ## 1     0     0     1     0     0     0     0     0     0
    ## 2     0     1     0     0     1     0     0     0     0
    ## 3     0     0     0     0     0     0     0     0     0
    ## 4     0     0     0     0     1     0     0     0     0
    ## 5     0     0     0     0     1     1     0     0     0
    ## 6     0     0     0     0     1     0     0     0     0

``` r
tradmetabr %>%
  group_by(bold_bin, classification) %>%
  summarise_all(sum) %>%
  mutate_if(is.numeric, ~1 * (. > 0)) %>% #change values >1 to 1
  data.frame() -> tradmetabr
```

    ## `mutate_if()` ignored the following grouping variables:
    ## * Column `bold_bin`

``` r
head(tradmetabr)
```

    ##       bold_bin
    ## 1 BOLD:AAA0009
    ## 2 BOLD:AAA0012
    ## 3 BOLD:AAA0016
    ## 4 BOLD:AAA0081
    ## 5 BOLD:AAA0082
    ## 6 BOLD:AAA0085
    ##                                                                                                    classification
    ## 1              k__Animalia;p__Arthropoda;c__Insecta;o__Lepidoptera;f__Gelechiidae;s__;g__gelJanzen1;s__Janzen1358
    ## 2                                  k__Animalia;p__Arthropoda;c__Insecta;o__Lepidoptera;f__Gelechiidae;s__;g__;s__
    ## 3 k__Animalia;p__Arthropoda;c__Insecta;o__Lepidoptera;f__Depressariidae;s__Stenomatinae;g__Cerconota;s__Janzen136
    ## 4   k__Animalia;p__Arthropoda;c__Insecta;o__Lepidoptera;f__Depressariidae;s__Stenomatinae;g__Stenoma;s__Janzen148
    ## 5              k__Animalia;p__Arthropoda;c__Insecta;o__Lepidoptera;f__Gelechiidae;s__;g__gelJanzen01;s__Janzen485
    ## 6 k__Animalia;p__Arthropoda;c__Insecta;o__Lepidoptera;f__Choreutidae;s__Choreutinae;g__Rhobonda;s__gaurisanaDHJ04
    ##   ALTOSCAMP_na_LT ARGOS_dry_BA ARGOS_wet_BA ARM1_dry_B ARM1_wet_B B05_wet_BB
    ## 1               0            0            0          0          0          0
    ## 2               0            0            0          0          0          0
    ## 3               0            0            0          0          0          0
    ## 4               0            0            0          0          0          0
    ## 5               0            0            0          0          0          0
    ## 6               0            0            0          0          0          0
    ##   ARM1_dry_BB ARM1_wet_BB ARM1_dry_BT ARM1_wet_BT ARM1_wet_LT ARM1_dry_LT
    ## 1           0           0           0           0           0           0
    ## 2           0           0           0           0           0           0
    ## 3           0           0           0           0           0           0
    ## 4           0           0           0           0           0           0
    ## 5           0           0           0           0           0           0
    ## 6           0           0           0           0           0           0
    ##   ARM1_dry_W ARM2_dry_B ARM2_wet_B TEST_dry_BB ARM2_dry_BB B06_wet_BB
    ## 1          0          0          0           0           0          0
    ## 2          0          0          0           0           0          0
    ## 3          0          0          0           0           0          0
    ## 4          0          0          0           0           0          0
    ## 5          0          0          0           0           0          0
    ## 6          0          0          0           0           0          0
    ##   ARM2_wet_BB ARM2_wet_BT ARM2_dry_BT ARM2_wet_LT ARM2_dry_LT ARM2_dry_W
    ## 1           0           0           0           0           0          0
    ## 2           0           0           0           0           0          0
    ## 3           0           0           0           0           0          0
    ## 4           0           0           0           0           0          0
    ## 5           0           0           0           0           0          0
    ## 6           0           0           0           0           0          0
    ##   ARM3_dry_B ARM3_wet_B ARM3_wet_BB ARM3_dry_BB ARM3_dry_BT ARM3_wet_BT
    ## 1          0          0           0           0           0           0
    ## 2          0          0           0           0           0           0
    ## 3          0          0           0           0           0           0
    ## 4          0          0           0           0           0           0
    ## 5          0          0           0           0           0           0
    ## 6          0          0           0           0           0           0
    ##   ARM3_wet_LT ARM3_dry_LT L02_wet_LT L02TEST_dry_LT L02_dry_LT M07_dry_MP
    ## 1           0           0          0              0          0          0
    ## 2           0           0          0              0          0          0
    ## 3           0           0          0              0          0          0
    ## 4           0           0          0              0          0          0
    ## 5           0           0          0              0          0          0
    ## 6           0           0          0              0          0          0
    ##   ARM3_dry_W ARM4_wet_B ARM4_dry_B B08_wet_BB ARM4_dry_BB ARM4_wet_BB
    ## 1          0          0          0          0           0           0
    ## 2          0          0          0          0           0           0
    ## 3          0          0          0          0           0           0
    ## 4          0          0          0          0           0           0
    ## 5          0          0          0          0           0           0
    ## 6          0          0          0          0           0           0
    ##   ARM4_dry_BT ARM4_wet_BT ARM4_wet_LT ARM4_dry_LT M08_dry_MP ARM4_dry_W
    ## 1           0           0           0           0          0          0
    ## 2           0           0           0           0          0          0
    ## 3           0           0           0           0          0          0
    ## 4           0           0           0           0          0          0
    ## 5           0           0           0           0          0          0
    ## 6           0           0           0           0          0          0
    ##   ATO_LESM_na_B B_na_BM B_na_TBA BAL1_dry_B BAL1_wet_B BAL1_dry_BB BAL1_wet_BB
    ## 1             0       0        0          0          0           0           0
    ## 2             0       0        0          0          0           0           0
    ## 3             0       0        0          0          0           0           0
    ## 4             0       0        0          0          0           0           0
    ## 5             0       0        0          0          0           0           0
    ## 6             0       0        0          0          0           0           0
    ##   BAL1_wet_BT BAL1_dry_BT BAL1_dry_LT BAL1_wet_LT M01_dry_MP BAL1_dry_W
    ## 1           0           0           0           0          0          0
    ## 2           0           0           0           0          0          0
    ## 3           0           0           0           0          0          0
    ## 4           0           0           0           0          0          0
    ## 5           0           0           0           0          0          0
    ## 6           0           0           0           0          0          0
    ##   ARBOREALBUJAN_na_AM ASPAN_na_B BCI_dry_BE BCI_wet_BM BCI_dry_BM BCI_na_HC
    ## 1                   0          0          0          0          0         0
    ## 2                   0          0          0          0          0         0
    ## 3                   0          0          0          0          0         0
    ## 4                   0          0          0          0          0         0
    ## 5                   0          0          0          0          0         0
    ## 6                   0          0          0          0          0         0
    ##   BCI_wet_HC BCI_dry_HC ALATESHIK_na_LT MANUALJESSE_na_MA ARBOREALENDARA_na_MA
    ## 1          0          0               0                 0                    0
    ## 2          0          0               0                 0                    0
    ## 3          0          0               0                 0                    0
    ## 4          0          0               0                 0                    0
    ## 5          0          0               0                 0                    0
    ## 6          0          0               0                 0                    0
    ##   MANUALDONOSOLIT_na_MA MANUALDONOSOUND_na_MA BCI_wet_MM BCI_dry_MM
    ## 1                     0                     0          0          0
    ## 2                     0                     0          0          0
    ## 3                     0                     0          0          0
    ## 4                     0                     0          0          0
    ## 5                     0                     0          0          0
    ## 6                     0                     0          0          0
    ##   FARAOC_dry_RP FARAOC_na_RP FROGDIETDONOSO_na_SC TR8_dry_T TR10_dry_T
    ## 1             0            0                    0         0          0
    ## 2             0            0                    0         0          0
    ## 3             0            0                    0         0          0
    ## 4             0            0                    0         0          0
    ## 5             0            0                    0         0          0
    ## 6             0            0                    0         0          0
    ##   TR11_dry_T TR2_dry_T TR7_dry_T TR9_dry_T TR1_dry_T TR6_dry_T BL_dry_BM
    ## 1          0         0         0         0         0         0         0
    ## 2          0         0         0         0         0         0         0
    ## 3          0         0         0         0         0         0         0
    ## 4          0         0         0         0         0         0         0
    ## 5          0         0         0         0         0         0         0
    ## 6          0         0         0         0         0         0         0
    ##   BL_dry_TH C37_dry_BA C37_wet_BA C42_wet_BA C42_dry_BA C48_dry_BA C48_wet_BA
    ## 1         0          0          0          0          0          0          0
    ## 2         0          0          0          0          0          0          0
    ## 3         0          0          0          0          0          0          0
    ## 4         0          0          0          0          0          0          0
    ## 5         0          0          0          0          0          0          0
    ## 6         0          0          0          0          0          0          0
    ##   CED_dry_T C57_wet_BA C57_dry_BA C64_wet_BA C64_dry_BA CAC_wet_HC
    ## 1         0          0          0          0          0          0
    ## 2         0          0          0          0          0          0
    ## 3         0          0          0          0          0          0
    ## 4         0          0          0          0          0          0
    ## 5         0          0          0          0          0          0
    ## 6         0          0          0          0          0          0
    ##   ALDERCOTTE_na_LA CHAC_LESM_na_B CHAN_LESM_na_B COM_LESM_na_B CUED_LESM_na_B
    ## 1                0              0              0             0              0
    ## 2                0              0              0             0              0
    ## 3                0              0              0             0              0
    ## 4                0              0              0             0              0
    ## 5                0              0              0             0              0
    ## 6                0              0              0             0              0
    ##   D_na_BM D_na_TBA DRA1_dry_B DRA1_wet_B DRA1_dry_BB DRA1_wet_BB DRA1_wet_BT
    ## 1       0        0          0          0           0           0           0
    ## 2       0        0          0          0           0           0           0
    ## 3       0        0          0          0           0           0           0
    ## 4       0        0          0          0           0           0           0
    ## 5       0        0          0          0           0           0           0
    ## 6       0        0          0          0           0           0           0
    ##   DRA1_dry_BT DRA1_dry_LT DRA1_wet_LT L03_wet_LT L03_dry_LT DRA1_dry_W
    ## 1           0           0           0          0          0          0
    ## 2           0           0           0          0          0          0
    ## 3           0           0           0          0          0          0
    ## 4           0           0           0          0          0          0
    ## 5           0           0           0          0          0          0
    ## 6           0           0           0          0          0          0
    ##   EMI_LESM_na_B FB_dry_BM FB_dry_TBA FM1_dry_BM FM2_dry_BM GIG_wet_TN
    ## 1             0         0          0          0          0          0
    ## 2             0         0          0          0          0          0
    ## 3             0         0          0          0          0          0
    ## 4             0         0          0          0          0          0
    ## 5             0         0          0          0          0          0
    ## 6             0         0          0          0          0          0
    ##   HER_LESM_na_B IXT_LESM_na_B JDH_dry_BM JVT_dry_BM JVT_dry_TBA JZ_dry_TBA
    ## 1             0             0          0          0           0          0
    ## 2             0             0          0          0           0          0
    ## 3             0             0          0          0           0          0
    ## 4             0             0          0          0           0          0
    ## 5             0             0          0          0           0          0
    ## 6             0             0          0          0           0          0
    ##   L_dry_BM NRG_LESM_na_B OPEN_na_HC PAN18FG062_na_HC PAN18FG009_na_HC
    ## 1        0             0          0                0                0
    ## 2        0             0          0                0                0
    ## 3        0             0          0                0                0
    ## 4        0             0          0                0                0
    ## 5        0             0          0                0                0
    ## 6        0             0          0                0                0
    ##   PAN18FG078_na_HC PAN18FG060_na_HC PAN18FG024_na_HC PAN18FG004_na_HC
    ## 1                0                0                0                0
    ## 2                0                0                0                0
    ## 3                0                0                0                0
    ## 4                0                0                0                0
    ## 5                0                0                0                0
    ## 6                0                0                0                0
    ##   PAN18FG080_na_HC PNM_wet_HC PRD_dry_HC PRD_wet_HC PRD_na_HC PRD_dry_LT
    ## 1                0          0          0          0         0          0
    ## 2                0          0          0          0         0          0
    ## 3                0          0          0          0         0          0
    ## 4                0          0          0          0         0          0
    ## 5                0          0          0          0         0          0
    ## 6                0          0          0          0         0          0
    ##   PUC_LESM_na_B RIOFRIO2_na_HC RIOFRIO1_na_MM SOL_dry_BM T38_wet_BA T38_dry_BA
    ## 1             0              0              0          0          0          0
    ## 2             0              0              0          0          0          0
    ## 3             0              0              0          0          0          0
    ## 4             0              0              0          0          0          0
    ## 5             0              0              0          0          0          0
    ## 6             0              0              0          0          0          0
    ##   T39_wet_BA T39_dry_BA TER_dry_T T46_dry_BA T46_wet_BA T47_wet_BA T47_dry_BA
    ## 1          0          0         0          0          0          0          0
    ## 2          0          0         0          0          0          0          0
    ## 3          0          0         0          0          0          0          0
    ## 4          0          0         0          0          0          0          0
    ## 5          0          0         0          0          0          0          0
    ## 6          0          0         0          0          0          0          0
    ##   T61_wet_BA T61_dry_BA TAP_LESM_na_B TAPIJ_LESM_na_B TB_dry_TBA TB1_dry_BM
    ## 1          0          0             0               0          0          0
    ## 2          0          0             0               0          0          0
    ## 3          0          0             0               0          0          0
    ## 4          0          0             0               0          0          0
    ## 5          0          0             0               0          0          0
    ## 6          0          0             0               0          0          0
    ##   TB2_dry_BM TB2_dry_TBA TEX_LESM_na_B TGP_dry_TBA TEK_dry_T TK2_wet_BA
    ## 1          0           0             0           0         0          0
    ## 2          0           0             0           0         0          0
    ## 3          0           0             0           0         0          0
    ## 4          0           0             0           0         0          0
    ## 5          0           0             0           0         0          0
    ## 6          0           0             0           0         0          0
    ##   TK3_wet_BA TK4_wet_BA TK5_wet_BA TUP_dry_HC TUX_LESM_na_B WHE1_wet_B
    ## 1          0          0          0          0             0          0
    ## 2          0          0          0          0             0          0
    ## 3          0          0          0          0             0          0
    ## 4          0          0          0          0             0          0
    ## 5          0          0          0          0             0          0
    ## 6          0          0          0          0             0          0
    ##   WHE1_dry_B WHE1_dry_BB WHE1_wet_BB B02_wet_BB WHE1_dry_BT WHE1_wet_BT
    ## 1          0           0           0          0           0           0
    ## 2          0           0           0          0           0           0
    ## 3          0           0           0          0           0           0
    ## 4          0           0           0          0           0           0
    ## 5          0           0           0          0           0           0
    ## 6          0           0           0          0           0           0
    ##   WHE1_wet_LT WHE1_dry_LT M03_dry_MP WHE1_dry_W WHE2_dry_B WHE2_wet_B
    ## 1           0           0          0          0          0          0
    ## 2           0           0          0          0          0          0
    ## 3           0           0          0          0          0          0
    ## 4           0           0          0          0          0          0
    ## 5           0           0          0          0          0          0
    ## 6           0           0          0          0          0          0
    ##   WHE2_dry_BB WHE2_dry_BT WHE2_wet_BT L05_wet_LT L05TEST_wet_LT WHE2_dry_LT
    ## 1           0           0           0          0              0           0
    ## 2           0           0           0          0              0           0
    ## 3           0           0           0          0              0           0
    ## 4           0           0           0          0              0           0
    ## 5           0           0           0          0              0           0
    ## 6           0           0           0          0              0           0
    ##   WHE2_wet_LT L05TEST_dry_LT WHE2_dry_W WMW_dry_BM WMW_dry_TBA ZET1_wet_B
    ## 1           0              0          0          0           0          0
    ## 2           0              0          0          0           0          0
    ## 3           0              0          0          0           0          0
    ## 4           0              0          0          0           0          0
    ## 5           0              0          0          0           0          0
    ## 6           0              0          0          0           0          0
    ##   ZET1_dry_B ZET1_dry_BB B09_wet_BB ZET1_wet_BB ZET1_wet_BT ZET1_dry_BT
    ## 1          0           0          0           0           0           0
    ## 2          0           0          0           0           0           0
    ## 3          0           0          0           0           0           0
    ## 4          0           0          0           0           0           0
    ## 5          0           0          0           0           0           0
    ## 6          0           0          0           0           0           0
    ##   L04_wet_LT ZET1_wet_LT ZET1_dry_LT L04TEST_dry_LT ZET1_dry_W ZET2_wet_B
    ## 1          0           0           0              0          0          0
    ## 2          0           0           0              0          0          0
    ## 3          0           0           0              0          0          0
    ## 4          0           0           0              0          0          0
    ## 5          0           0           0              0          0          0
    ## 6          0           0           0              0          0          0
    ##   ZET2_dry_B ZET2_dry_BB B10_wet_BB ZET2_wet_BT ZET2_dry_BT ZET2_wet_LT
    ## 1          0           0          0           0           0           0
    ## 2          0           0          0           0           0           0
    ## 3          0           0          0           0           0           0
    ## 4          0           0          0           0           0           0
    ## 5          0           0          0           0           0           0
    ## 6          0           0          0           0           0           0
    ##   ZET2_dry_LT L01_wet_LT L01TEST_dry_LT L01_dry_LT M10_dry_MP ZET2_dry_W
    ## 1           0          0              0          0          0          0
    ## 2           0          0              0          0          0          0
    ## 3           0          0              0          0          0          0
    ## 4           0          0              0          0          0          0
    ## 5           0          0              0          0          0          0
    ## 6           0          0              0          0          0          0
    ##   WINKLERDONOSO_na_W total39864 total39865 total39866 total39867 total39868
    ## 1                  0          0          0          0          0          0
    ## 2                  0          0          0          0          0          0
    ## 3                  0          0          0          0          0          0
    ## 4                  0          0          0          1          0          0
    ## 5                  0          0          0          0          0          0
    ## 6                  0          0          0          0          0          0
    ##   total39869 total39870 total39871 total39872 total39873 total39874 total39875
    ## 1          0          0          0          0          0          0          0
    ## 2          0          0          0          0          0          0          0
    ## 3          0          0          0          0          0          0          0
    ## 4          0          0          0          0          0          0          0
    ## 5          0          0          0          0          0          0          0
    ## 6          0          0          0          0          0          0          0
    ##   total39876 total39877 total39878 total39879 total39880 total39881 total39882
    ## 1          0          0          0          0          0          0          0
    ## 2          0          0          0          0          0          0          0
    ## 3          0          0          0          0          0          0          0
    ## 4          0          0          0          0          0          0          0
    ## 5          0          0          0          0          0          0          0
    ## 6          0          0          0          0          0          0          0
    ##   total39883 ARM1A ARM1B ARM2A ARM2B ARM3A ARM3B ARM4A ARM4B BAL1A BAL1B DRA1A
    ## 1          0     0     0     1     0     0     0     0     0     0     0     0
    ## 2          0     0     0     0     0     0     0     0     0     0     1     0
    ## 3          0     0     0     0     0     0     0     1     0     0     0     0
    ## 4          0     0     0     0     0     0     1     0     0     0     0     0
    ## 5          0     0     0     0     0     0     0     0     0     0     0     0
    ## 6          0     0     0     0     0     0     0     0     0     0     0     0
    ##   DRA1B WHE1A WHE1B WHE2A WHE2B ZET1A ZET1B ZET2A ZET2B
    ## 1     0     0     1     0     0     0     0     0     0
    ## 2     0     1     0     0     1     0     0     0     0
    ## 3     0     0     0     0     0     0     0     0     0
    ## 4     0     0     0     0     1     0     0     0     0
    ## 5     0     0     0     0     1     1     0     0     0
    ## 6     0     0     0     0     1     0     0     0     0

``` r
write.csv(tradmetabr, 'data/tradmetabr_merged.csv', row.names =FALSE) #save final dataset
```

tradmetabr_merged.csv is now the final file we will use in script 05.
