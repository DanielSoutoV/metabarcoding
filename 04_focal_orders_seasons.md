04_focal_order_seasons
================
Daniel
06/10/2022

``` r
obj <- readRDS('data/taxmap_object.rds') #loads the taxmap object created in script 02_metacoder_heat_trees

obj %>%  metacoder::filter_taxa(taxon_names %in% c("Lepidoptera"),#here is to fliter the figure by groups
              subtaxa = TRUE) -> leps #we will create separate files for each order as this simplifies downstream analysis in microbiotaprocess - until I find a way to filter taxa on the mpse object directly.

obj %>%  metacoder::filter_taxa(taxon_names %in% c("Coleoptera"),#here is to fliter the figure by groups
              subtaxa = TRUE) -> coleo
obj %>%  metacoder::filter_taxa(taxon_names %in% c("Diptera"),#here is to fliter the figure by groups
              subtaxa = TRUE) -> dips
obj %>%  metacoder::filter_taxa(taxon_names %in% c("Hymenoptera"),#here is to fliter the figure by groups
              subtaxa = TRUE) -> bees #i know its hymenoptera
obj %>%  metacoder::filter_taxa(taxon_names %in% c("Hemiptera"),#here is to fliter the figure by groups
              subtaxa = TRUE) -> hemi
obj %>%  metacoder::filter_taxa(taxon_names %in% c("Blattodea"),#here is to fliter the figure by groups
              subtaxa = TRUE) -> blats 

sample <- read.csv('data/location_ctrl.csv')
```

The script above gave me the taxmap object separated by order, however,
it removes all taxonomic ranks up until “Order”, so the script below is
firstly, changing the taxmap object (e.g. leps) into a phyloseq object
(as_phyloseq) and then using the microViz package to add the taxonomic
ranks. These are necessary for some of the microbiotaprocess analyses
otherwise the package complains that there is not enough taxonomic
ranks.

There must be certainly a cleaner way how to do this without so many
steps, but I then had to change the classification rank names again
(“k\_\_” instead of Kingdom) as this is more friendly with
MicroBiotaProcess. Within the metabarcoding package community, there
seems to have been several issues with the ranks, some have root, some
don’t, some have subfamily, some don’t, so take note of the nomenclature
you are using. I am using
Kingdom;Phylum;Order;Class;Family;genus;species (1 - 7).

If someone knows how to write a loop to make this faster, I’ll buy you a
beer.

Ok we now have datasets separated for each order (ForestGEO focal plus
diptera - if we want more, we can adapt the script above to include more
focal groups, or even to refine them - e.g. family level) But below we
will create now some taxonomic trees showing differences between seasons
(as in the previous script).

``` r
set.seed(1024)
deresleps <- diff_analysis(obj = ps_leps, classgroup = "SEASON",
                       mlfun = "lda",
                       filtermod = "pvalue",
                       firstcomfun = "kruskal_test",
                       firstalpha = 0.05,
                       strictmod = TRUE,
                       secondcomfun = "wilcox_test",
                       subclmin = 3,
                       subclwilc = TRUE,
                       secondalpha = 0.01,
                       lda=3,
                       action = "add")
deresleps
```

    ## The original data: 1597 features and 40 samples
    ## The sample data: 1 variables and 40 samples
    ## The taxda contained 1170 by 7 rank
    ## after first test (kruskal_test) number of feature (pvalue<=0.05):352
    ## after second test (wilcox_test and generalizedFC) number of significantly discriminative feature:219
    ## after lda, Number of discriminative features: 159 (certain taxonomy classification:143; uncertain taxonomy classication: 16)

``` r
diffclade_leps <- ggdiffclade(
  obj=deresleps, 
  alpha=0.5, 
  linewd=0.01,
  skpointsize=0.05, 
  layout="radial",
  cladetext = 0.4,
  taxlevel=7, #taxonomy level from 1 to 7 kingdom:phylum:class:order:family:genus:species
  removeUnkown=TRUE,
  reduce=TRUE # This argument is to remove the branch of unknown taxonomy.
) +
  scale_fill_manual(
    values=c("goldenrod", "steelblue")
  ) +
  guides(color = guide_legend(
    keywidth = 0.1, 
    keyheight = 0.6,
    order = 5,
    ncol=3)
  ) +
  theme(
    panel.background=element_rect(fill=NA),
    legend.position="right", 
    plot.margin=margin(0,0,0,0),
    legend.spacing.y=unit(0.02, "cm"), 
    legend.title=element_text(size=7),
    legend.text=element_text(size=6), 
    legend.box.spacing=unit(0.02,"cm")
  )
```

    ## The `removeUnkown` has been deprecated, Please use `removeUnknown` instead!

    ## The color has been set automatically, you can reset it manually by adding scale_fill_manual(values=yourcolors)

    ## Scale for 'fill' is already present. Adding another scale for 'fill', which
    ## will replace the existing scale.

``` r
pdf("./04_focal_orders_seasons_files/leps_cladogram.pdf")
diffclade_leps
dev.off()
```

    ## png 
    ##   2

``` r
classtaxa_leps <- get_taxadf(obj=ps_leps, taxlevel=5) #taxonomy level from 1 to 7 kingdom:phylum:class:order:family:genus:species
# The 10 most abundant taxonomy will be visualized by default (parameter `topn=10`). 
pclass_leps <- ggbartax(obj=classtaxa_leps, facetNames="SEASON", topn=10) +
  xlab(NULL) +
  ylab("relative abundance (%)") +
  scale_fill_manual(values=c(colorRampPalette(RColorBrewer::brewer.pal(11,"Set3"))(11))) +
  guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5))
```

    ## The color has been set automatically, you can reset it 
    ##             manually by adding scale_fill_manual(values=yourcolors)
    ## Scale for 'fill' is already present. Adding another scale for 'fill', which
    ## will replace the existing scale.

``` r
pdf("./04_focal_orders_seasons_files/leps_relative_difference_season.pdf")
pclass_leps
dev.off()
```

    ## png 
    ##   2

``` r
ps_leps %>% as.MPSE() %>% mp_rrarefy() %>% mp_diff_analysis(.abundance=RareAbundance, .group=SEASON, action='get') %>% dplyr::filter(grepl("^g_", f)) %>% ggdiffbox(colorlist=c("goldenrod", "steelblue"), notch = FALSE) -> ggdiffbox_leps
```

    ## The otutree is empty in the MPSE object!

    ## The color has been set automatically, you can reset it manually by adding scale_color_manual(values=yourcolors)

    ## Scale for 'colour' is already present. Adding another scale for 'colour',
    ## which will replace the existing scale.

``` r
pdf("./04_focal_orders_seasons_files/leps_diffbox.pdf") 
ggdiffbox_leps
```

    ## notch went outside hinges. Try setting notch=FALSE.

    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.

``` r
dev.off()
```

    ## png 
    ##   2

``` r
es_leps <- ggeffectsize(obj=deresleps, 
                     lineheight=1,
                     linewidth=0.01, ytextsize = 0.1) + 
  scale_color_manual(values=c("goldenrod", "steelblue")) 
```

    ## The color has been set automatically, you can reset it manually by adding scale_color_manual(values=yourcolors)

    ## Scale for 'colour' is already present. Adding another scale for 'colour',
    ## which will replace the existing scale.

``` r
pdf("./04_focal_orders_seasons_files/leps_effectsize.pdf") 
es_leps + theme(axis.text.y = element_text(size=8, angle=30))
dev.off()
```

    ## png 
    ##   2

``` r
#I will buy a beer to whoever manages to fix the y axis lables to be smaller and better spaced.
```

``` r
set.seed(1024)
derescoleo <- diff_analysis(obj = ps_coleo, classgroup = "SEASON",
                       mlfun = "lda",
                       filtermod = "pvalue",
                       firstcomfun = "kruskal_test",
                       firstalpha = 0.05,
                       strictmod = TRUE,
                       secondcomfun = "wilcox_test",
                       subclmin = 3,
                       subclwilc = TRUE,
                       secondalpha = 0.01,
                       lda=3,
                       action = "add")
derescoleo
```

    ## The original data: 275 features and 40 samples
    ## The sample data: 1 variables and 40 samples
    ## The taxda contained 272 by 7 rank
    ## after first test (kruskal_test) number of feature (pvalue<=0.05):87
    ## after second test (wilcox_test and generalizedFC) number of significantly discriminative feature:62
    ## after lda, Number of discriminative features: 62 (certain taxonomy classification:43; uncertain taxonomy classication: 19)

``` r
diffclade_coleo <- ggdiffclade(
  obj=derescoleo, 
  alpha=0.5, 
  linewd=0.01,
  skpointsize=0.5, 
  layout="radial",
  cladetext = 0.4,
  taxlevel=7, #taxonomy level from 1 to 7 kingdome:phylum:class:order:family:genus:species
  removeUnkown=TRUE,
  reduce=TRUE # This argument is to remove the branch of unknown taxonomy.
) +
  scale_fill_manual(
    values=c("goldenrod", "steelblue")
  ) +
  guides(color = guide_legend(
    keywidth = 0.1, 
    keyheight = 0.6,
    order = 5,
    ncol=3)
  ) +
  theme(
    panel.background=element_rect(fill=NA),
    legend.position="right", 
    plot.margin=margin(0,0,0,0),
    legend.spacing.y=unit(0.02, "cm"), 
    legend.title=element_text(size=7),
    legend.text=element_text(size=6), 
    legend.box.spacing=unit(0.02,"cm")
  )
```

    ## The `removeUnkown` has been deprecated, Please use `removeUnknown` instead!

    ## The color has been set automatically, you can reset it manually by adding scale_fill_manual(values=yourcolors)

    ## Scale for 'fill' is already present. Adding another scale for 'fill', which
    ## will replace the existing scale.

``` r
pdf("./04_focal_orders_seasons_files/coleo_cladogram.pdf")
diffclade_coleo
dev.off
```

    ## function (which = dev.cur()) 
    ## {
    ##     if (which == 1) 
    ##         stop("cannot shut down device 1 (the null device)")
    ##     .External(C_devoff, as.integer(which))
    ##     dev.cur()
    ## }
    ## <bytecode: 0x000000001435caf8>
    ## <environment: namespace:grDevices>

``` r
classtaxa_coleo <- get_taxadf(obj=ps_coleo, taxlevel=5)
# The 10 most abundant taxonomy will be visualized by default (parameter `topn=10`). 
pclass_coleo <- ggbartax(obj=classtaxa_coleo, facetNames="SEASON", topn=10) +
  xlab(NULL) +
  ylab("relative abundance (%)") +
  scale_fill_manual(values=c(colorRampPalette(RColorBrewer::brewer.pal(11,"Set3"))(11))) +
  guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5))
```

    ## The color has been set automatically, you can reset it 
    ##             manually by adding scale_fill_manual(values=yourcolors)
    ## Scale for 'fill' is already present. Adding another scale for 'fill', which
    ## will replace the existing scale.

``` r
pdf("./04_focal_orders_seasons_files/coleo_diffclass_top5.pdf")
pclass_coleo
dev.off()
```

    ## png 
    ##   2

``` r
ps_coleo %>% as.MPSE() %>% mp_rrarefy() %>% mp_diff_analysis(.abundance=RareAbundance, .group=SEASON, action='get') %>% dplyr::filter(grepl("^f_", f)) %>% ggdiffbox(colorlist=c("steelblue", "goldenrod"), notch = FALSE) -> coleo_ggdiffbox
```

    ## The otutree is empty in the MPSE object!

    ## The color has been set automatically, you can reset it manually by adding scale_color_manual(values=yourcolors)

    ## Scale for 'colour' is already present. Adding another scale for 'colour',
    ## which will replace the existing scale.

``` r
pdf("./04_focal_orders_seasons_files/coleo_ggdiffbox")
coleo_ggdiffbox
```

    ## notch went outside hinges. Try setting notch=FALSE.

    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.

``` r
dev.off()
```

    ## png 
    ##   2

``` r
es_coleo <- ggeffectsize(obj=derescoleo, 
                     lineheight=0.1,
                     linewidth=0.3, ytextsize = 0.3) + 
  scale_color_manual(values=c("goldenrod", "steelblue")) 
```

    ## The color has been set automatically, you can reset it manually by adding scale_color_manual(values=yourcolors)

    ## Scale for 'colour' is already present. Adding another scale for 'colour',
    ## which will replace the existing scale.

``` r
pdf("./04_focal_orders_seasons_files/coleo_effectsize.pdf")
es_coleo
dev.off()
```

    ## png 
    ##   2

``` r
set.seed(1024)
deresdips <- diff_analysis(obj = ps_dips, classgroup = "SEASON",
                       mlfun = "lda",
                       filtermod = "pvalue",
                       firstcomfun = "kruskal_test",
                       firstalpha = 0.05,
                       strictmod = TRUE,
                       secondcomfun = "wilcox_test",
                       subclmin = 3,
                       subclwilc = TRUE,
                       secondalpha = 0.01,
                       lda=3,
                       action = "add")
deresdips
```

    ## The original data: 259 features and 40 samples
    ## The sample data: 1 variables and 40 samples
    ## The taxda contained 511 by 7 rank
    ## after first test (kruskal_test) number of feature (pvalue<=0.05):73
    ## after second test (wilcox_test and generalizedFC) number of significantly discriminative feature:55
    ## after lda, Number of discriminative features: 51 (certain taxonomy classification:28; uncertain taxonomy classication: 23)

``` r
diffclade_dips <- ggdiffclade(
  obj=deresdips, 
  alpha=0.5, 
  linewd=0.01,
  skpointsize=0.5, 
  layout="radial",
  cladetext = 0.4,
  taxlevel=5, #taxonomy level from 1 to 7 kingdome:phylum:class:order:family:subfamily:genus:species
  removeUnkown=TRUE,
  reduce=TRUE # This argument is to remove the branch of unknown taxonomy.
) +
  scale_fill_manual(
    values=c("goldenrod", "steelblue")
  ) +
  guides(color = guide_legend(
    keywidth = 0.1, 
    keyheight = 0.6,
    order = 5,
    ncol=3)
  ) +
  theme(
    panel.background=element_rect(fill=NA),
    legend.position="right", 
    plot.margin=margin(0,0,0,0),
    legend.spacing.y=unit(0.02, "cm"), 
    legend.title=element_text(size=7),
    legend.text=element_text(size=6), 
    legend.box.spacing=unit(0.02,"cm")
  )
```

    ## The `removeUnkown` has been deprecated, Please use `removeUnknown` instead!

    ## The color has been set automatically, you can reset it manually by adding scale_fill_manual(values=yourcolors)

    ## Scale for 'fill' is already present. Adding another scale for 'fill', which
    ## will replace the existing scale.

``` r
pdf("./04_focal_orders_seasons_files/dips_clades.pdf")
diffclade_dips
dev.off()
```

    ## png 
    ##   2

``` r
classtaxa_dips <- get_taxadf(obj=ps_dips, taxlevel=5)
# The 10 most abundant taxonomy will be visualized by default (parameter `topn=10`). 
pclass_dips <- ggbartax(obj=classtaxa_dips, facetNames="SEASON", topn=10) +
  xlab(NULL) +
  ylab("relative abundance (%)") +
  scale_fill_manual(values=c(colorRampPalette(RColorBrewer::brewer.pal(11,"Set3"))(11))) +
  guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5))
```

    ## The color has been set automatically, you can reset it 
    ##             manually by adding scale_fill_manual(values=yourcolors)
    ## Scale for 'fill' is already present. Adding another scale for 'fill', which
    ## will replace the existing scale.

``` r
pdf("./04_focal_orders_seasons_files/dips_top5.pdf")
pclass_dips
dev.off()
```

    ## png 
    ##   2

``` r
ps_dips %>% as.MPSE() %>% mp_rrarefy() %>% mp_diff_analysis(.abundance=RareAbundance, .group=SEASON, action='get') %>% dplyr::filter(grepl("^f_", f)) %>% ggdiffbox(colorlist=c("goldenrod", "steelblue"), notch = FALSE) -> diptera_ggdiffbox
```

    ## The otutree is empty in the MPSE object!

    ## The color has been set automatically, you can reset it manually by adding scale_color_manual(values=yourcolors)

    ## Scale for 'colour' is already present. Adding another scale for 'colour',
    ## which will replace the existing scale.

``` r
pdf("./04_focal_orders_seasons_files/diptera_diffbox.pdf")
diptera_ggdiffbox
```

    ## notch went outside hinges. Try setting notch=FALSE.

    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.

``` r
dev.off()
```

    ## png 
    ##   2

``` r
es_dips <- ggeffectsize(obj=deresdips, 
                     lineheight=0.1,
                     linewidth=0.3, ytextsize = 0.3) + 
  scale_color_manual(values=c("goldenrod", "steelblue")) 
```

    ## The color has been set automatically, you can reset it manually by adding scale_color_manual(values=yourcolors)

    ## Scale for 'colour' is already present. Adding another scale for 'colour',
    ## which will replace the existing scale.

``` r
pdf("./04_focal_orders_seasons_files/diptera_effectsize.pdf")
es_dips
dev.off()
```

    ## png 
    ##   2

``` r
set.seed(1024)
dereshemi <- diff_analysis(obj = ps_hemi, classgroup = "SEASON",
                       mlfun = "lda",
                       filtermod = "pvalue",
                       firstcomfun = "kruskal_test",
                       firstalpha = 0.05,
                       strictmod = TRUE,
                       secondcomfun = "wilcox_test",
                       subclmin = 3,
                       subclwilc = TRUE,
                       secondalpha = 0.01,
                       lda=3,
                       action = "add")
dereshemi
```

    ## The original data: 158 features and 40 samples
    ## The sample data: 1 variables and 40 samples
    ## The taxda contained 114 by 7 rank
    ## after first test (kruskal_test) number of feature (pvalue<=0.05):39
    ## after second test (wilcox_test and generalizedFC) number of significantly discriminative feature:29
    ## after lda, Number of discriminative features: 29 (certain taxonomy classification:17; uncertain taxonomy classication: 12)

``` r
diffclade_hemi <- ggdiffclade(
  obj=dereshemi, 
  alpha=0.5, 
  linewd=0.01,
  skpointsize=0.5, 
  layout="radial",
  cladetext = 0.4,
  taxlevel=1, #taxonomy level from 1 to 7 kingdome:phylum:class:order:family:subfamily:genus:species
  removeUnkown=TRUE,
  reduce=TRUE # This argument is to remove the branch of unknown taxonomy.
) +
  scale_fill_manual(
    values=c("steelblue", "goldenrod" ) #Note that I had to change the order here, there are no hemiptera differences in DRY season so they did not come out in plot
  ) +
  guides(color = guide_legend(
    keywidth = 0.1, 
    keyheight = 0.6,
    order = 5,
    ncol=3)
  ) +
  theme(
    panel.background=element_rect(fill=NA),
    legend.position="right", 
    plot.margin=margin(0,0,0,0),
    legend.spacing.y=unit(0.02, "cm"), 
    legend.title=element_text(size=7),
    legend.text=element_text(size=6), 
    legend.box.spacing=unit(0.02,"cm")
  )
```

    ## The `removeUnkown` has been deprecated, Please use `removeUnknown` instead!

    ## The color has been set automatically, you can reset it manually by adding scale_fill_manual(values=yourcolors)

    ## Scale for 'fill' is already present. Adding another scale for 'fill', which
    ## will replace the existing scale.

``` r
pdf("./04_focal_orders_seasons_files/hemiptera_cladogram.pdf")
diffclade_hemi
dev.off()
```

    ## png 
    ##   2

``` r
classtaxa_hemi <- get_taxadf(obj=ps_hemi, taxlevel=5)
# The 10 most abundant taxonomy will be visualized by default (parameter `topn=10`). 
pclass_hemi <- ggbartax(obj=classtaxa_hemi, facetNames="SEASON", topn=10) +
  xlab(NULL) +
  ylab("relative abundance (%)") +
  scale_fill_manual(values=c(colorRampPalette(RColorBrewer::brewer.pal(11,"Set3"))(11))) +
  guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5))
```

    ## The color has been set automatically, you can reset it 
    ##             manually by adding scale_fill_manual(values=yourcolors)
    ## Scale for 'fill' is already present. Adding another scale for 'fill', which
    ## will replace the existing scale.

``` r
pdf("./04_focal_orders_seasons_files/hemiptera_top5.pdf")
pclass_hemi
dev.off()
```

    ## png 
    ##   2

``` r
ps_hemi %>% as.MPSE() %>% mp_rrarefy() %>% mp_diff_analysis(.abundance=RareAbundance, .group=SEASON, action='get') %>% dplyr::filter(grepl("^f_", f)) %>% ggdiffbox(colorlist=c("goldenrod", "steelblue"), notch = FALSE) -> hemiptera_ggdiffbox
```

    ## The otutree is empty in the MPSE object!

    ## The color has been set automatically, you can reset it manually by adding scale_color_manual(values=yourcolors)

    ## Scale for 'colour' is already present. Adding another scale for 'colour',
    ## which will replace the existing scale.

``` r
pdf("./04_focal_orders_seasons_files/hemiptera_diffbox.pdf")
hemiptera_ggdiffbox
```

    ## notch went outside hinges. Try setting notch=FALSE.

    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.

``` r
dev.off()
```

    ## png 
    ##   2

``` r
es_hemi <- ggeffectsize(obj=dereshemi, 
                     lineheight=0.1,
                     linewidth=0.3, ytextsize = 0.3) + 
  scale_color_manual(values=c( "steelblue", "goldenrod")) #same change due to wrong plotting 
```

    ## The color has been set automatically, you can reset it manually by adding scale_color_manual(values=yourcolors)

    ## Scale for 'colour' is already present. Adding another scale for 'colour',
    ## which will replace the existing scale.

``` r
pdf("./04_focal_orders_seasons_files/hemi_effectsize.pdf")
es_hemi
dev.off()
```

    ## png 
    ##   2

``` r
classtaxa_bees <- get_taxadf(obj=ps_bees, taxlevel=5)
# The 10 most abundant taxonomy will be visualized by default (parameter `topn=10`). 
pclass_bees <- ggbartax(obj=classtaxa_bees, facetNames="SEASON", topn=10) +
  xlab(NULL) +
  ylab("relative abundance (%)") +
  scale_fill_manual(values=c(colorRampPalette(RColorBrewer::brewer.pal(11,"Set3"))(11))) +
  guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5))
```

    ## The color has been set automatically, you can reset it 
    ##             manually by adding scale_fill_manual(values=yourcolors)

    ## Scale for 'fill' is already present. Adding another scale for 'fill', which
    ## will replace the existing scale.

``` r
pdf("./04_focal_orders_seasons_files/hymenoiptera_top10.pdf")
pclass_bees
```

    ## Warning: Removed 10 rows containing missing values (position_stack).

``` r
dev.off()
```

    ## png 
    ##   2

``` r
classtaxa_blats <- get_taxadf(obj=ps_blats, taxlevel=5)
# The 10 most abundant taxonomy will be visualized by default (parameter `topn=10`). 
pclass_blats <- ggbartax(obj=classtaxa_blats, facetNames="SEASON", topn=5) +
  xlab(NULL) +
  ylab("relative abundance (%)") +
  scale_fill_manual(values=c(colorRampPalette(RColorBrewer::brewer.pal(11,"Set3"))(11))) +
  guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5))
```

    ## The color has been set automatically, you can reset it 
    ##             manually by adding scale_fill_manual(values=yourcolors)
    ## Scale for 'fill' is already present. Adding another scale for 'fill', which
    ## will replace the existing scale.

``` r
pdf("./04_focal_orders_seasons_files/blats_top5.pdf")
pclass_blats
```

    ## Warning: Removed 5 rows containing missing values (position_stack).

``` r
dev.off()
```

    ## png 
    ##   2

############################ 

There are problems with the Hymenoptera and Blattodea datasets. I still
must figure it out. I suspect its because there is missing data for one
site for each group DRA1B has no Hymenoptera and ZET1A no termites

one more beer to whoever can figure it out!!!!!!!!!!!!!!

################################## 

r Hymenoptera set.seed(1024) deresbees \<- diff_analysis(obj = ps_bees,
classgroup = “SEASON”, mlfun = “lda”, filtermod = “pvalue”, firstcomfun
= “kruskal_test”, firstalpha = 0.05, strictmod = TRUE, secondcomfun =
“wilcox_test”, subclmin = 3, subclwilc = TRUE, secondalpha = 0.01,
lda=3, action = “add”) deresbees

diffclade_bees \<- ggdiffclade( obj=deresbees, alpha=0.5, linewd=0.15,
skpointsize=0.2, layout=“radial”, cladetext = 0.7, taxlevel=5, #taxonomy
level from 1 to 7
kingdome:phylum:class:order:family:subfamily:genus:species
removeUnkown=TRUE, reduce=TRUE # This argument is to remove the branch
of unknown taxonomy. ) + scale_fill_manual( values=c(“goldenrod”,
“steelblue”) ) + guides(color = guide_legend( keywidth = 0.1, keyheight
= 0.6, order = 5, ncol=3) ) + theme(
panel.background=element_rect(fill=NA), legend.position=“right”,
plot.margin=margin(0,0,0,0), legend.spacing.y=unit(0.02, “cm”),
legend.title=element_text(size=7), legend.text=element_text(size=6),
legend.box.spacing=unit(0.02,“cm”) ) diffclade_bees

ps_bees %>% as.MPSE() %>% mp_rrarefy() %>%
mp_diff_analysis(.abundance=RareAbundance, .group=SEASON, action=‘get’)
%>% dplyr::filter(grepl(“^f\_”, f)) %>%
ggdiffbox(colorlist=c(“steelblue”, “goldenrod”), notch = FALSE)

es_bees \<- ggeffectsize(obj=deresbees, lineheight=0.1, linewidth=0.3,
ytextsize = 0.3) + scale_color_manual(values=c(“goldenrod”,
“steelblue”))

es_bees

r Blattodea set.seed(1024) deresblats \<- diff_analysis(obj = ps_blats,
classgroup = “SEASON”, mlfun = “lda”, filtermod = “pvalue”, firstcomfun
= “kruskal_test”, firstalpha = 0.05, strictmod = TRUE, secondcomfun =
“wilcox_test”, subclmin = 3, subclwilc = TRUE, secondalpha = 0.01,
lda=3, action = “add”) deresblats

diffclade_blats \<- ggdiffclade( obj=deresblats, alpha=0.5, linewd=0.15,
skpointsize=0.2, layout=“radial”, cladetext = 0.7, taxlevel=6, #taxonomy
level from 1 to 8
kingdome:phylum:class:order:family:subfamily:genus:species
removeUnkown=TRUE, reduce=TRUE # This argument is to remove the branch
of unknown taxonomy. ) + scale_fill_manual( values=c(“goldenrod”,
“steelblue”) ) + guides(color = guide_legend( keywidth = 0.1, keyheight
= 0.6, order = 5, ncol=3) ) + theme(
panel.background=element_rect(fill=NA), legend.position=“right”,
plot.margin=margin(0,0,0,0), legend.spacing.y=unit(0.02, “cm”),
legend.title=element_text(size=7), legend.text=element_text(size=6),
legend.box.spacing=unit(0.02,“cm”) ) diffclade_blats

ps_blats %>% as.MPSE() %>% mp_rrarefy() %>%
mp_diff_analysis(.abundance=RareAbundance, .group=SEASON, action=‘get’)
%>% dplyr::filter(grepl(“^st\_”, f)) %>%
ggdiffbox(colorlist=c(“steelblue”, “goldenrod”), notch = FALSE)

es_blats \<- ggeffectsize(obj=deresblats, lineheight=0.1, linewidth=0.3,
ytextsize = 0.3) + scale_color_manual(values=c(“goldenrod”,
“steelblue”))

es_blats
