# Visualization with R
Here I've included a number of R scripts for working with R data.

- [Venn diagrams](#venn-diagrams)
  - [1. Read in the data](#1-read-in-the-data)
  - [2. Summarize by environment (site-year combination)](#2-summarize-by-environment--site-year-combination-)
  - [3. Transpose, relabel column names, and convert to presence-absence](#3-transpose--relabel-column-names--and-convert-to-presence-absence)
  - [4. Venn diagram for the canola root microbiome.](#4-venn-diagram-for-the-canola-root-microbiome)
  
- [Radial phylogenetic trees with relative abundances](#radial-phylogenetic-trees-with-relative-abundances)
   - [1. Read in the data](#1-read-in-the-data-1)
   - [2. Deal with R adding Xs to the abundance df during processing.](#2-deal-with-r-adding-xs-to-the-abundance-df-during-processing)
   - [3. Create a subsetting vector to use to prune the tree to taxa of interest](#3-create-a-subsetting-vector-to-use-to-prune-the-tree-to-taxa-of-interest)
   - [4. Prepare a df to use for tree-plotting](#4-prepare-a-df-to-use-for-tree-plotting)
   - [5. Now to create the phylogenetic tree to plot and add sample information](#5-now-to-create-the-phylogenetic-tree-to-plot-and-add-sample-information)

---

### Venn diagrams

``````
library(vegan)          # Used to convert abundance data to presence-absence
library(tidyverse)      # Load a bunch of tidy packages (here we're using dplyr)
library(eulerr)         # This is the best package I've found for generating Venn diagrams
library(scales)         # Get transparency in colours

rm(list = ls())         # I like to start with a clean workspace
``````

#### 1. Read in the data

`
root.asv <- read.csv("~/Dropbox/r code repository/Data/root.asv.sub.csv", header = T, row.names = 1)
`

The data are available [here.](https://www.dropbox.com/s/zm2l1p7j0yoe5x7/root.asv.sub.csv?dl=0)

These data are a subset of 16S sequence counts for bacteria from canola roots. The row names are sample information. The first column is a factor that indicates the site and year of the sample. Each site-year combination represents an environment we are interested in (n = 4 environments).

#### 2. Summarize by environment (site-year combination)

Changes from 360 x 1001 to 4 x 1001

````
venn.df <- root.asv %>% 
  group_by(SiteYear) %>% 
  summarise_each(list(mean))
dim(venn.df) # 4 x 7562
````

In this step, we've used a tidy notation "piping" (`%>%`). Use `%>%` to emphasise a sequence of actions, rather than the object that the actions are being performed on. Here, we've told R to take the root.asv df, group by the SiteYear factor, and calculate the mean abundance for each taxon. This yields a df with 4 rows (site-year) and 1000 columns (1000 ASVs).

Now we can strip the SiteYear column, transpose, rename the remaining columns, and convert to presence-absence to use for the Venn diagram

#### 3. Transpose, relabel column names, and convert to presence-absence
```
venn.df2 <- as.data.frame(t(venn.df[-1]))
names(venn.df2) <- c("L.2016","L.2017","M.2017","S.2017")
venn.df3 <- decostand(venn.df2, "pa")
```

#### 4. Venn diagram for the canola root microbiome.
`````
set.seed(10)
venn.fit <- euler(venn.df3)
eulerr_options(pointsize = 15)
plot(venn.fit, fills = scales::alpha(c("cadetblue2","darkorchid1","darkseagreen2","khaki2"),0.5), quantities = T, lty = 0, legend = list(labels = c("L.2016","L.2017","M.2017","S.2017")))
`````
<img src="https://user-images.githubusercontent.com/44586553/68962285-d39e9e80-0799-11ea-834b-ba9aa8669793.jpg" width="400" height="300">

---

### Radial phylogenetic trees with relative abundances

`````````````
library(ape)          # drop.tip function
library(picante)      # read.tree
library(ggplot2)      # Required for tree presentation and manipulation
library(ggimage)      # Required for tree presentation and manipulation
library(ggtree)       # Required for tree presentation and manipulation
library(picante)      # Required for tree presentation and manipulation
library(phytools)     # Required for tree presentation and manipulation
library(adephylo)     # Required for tree presentation and manipulation
library(phylobase)    # Required for tree presentation and manipulation
library(phylosignal)  # Needed for phylo4d plotting
library(vegan)        # decostand

rm(list = ls())
`````````````

#### 1. Read in the data

Tree constructed using fragment insertion in [QIIME2:](https://library.qiime2.org/plugins/q2-fragment-insertion/16/)

`
root.tree <- read.tree("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Rooted_tree/InsertionTree/insertion-tree.nwk")
`

A list of ASVs we'd like to use to subset the tree:

`
asv.int.imp1 <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Bal.Tab.intersect.imp1.csv", header = T, row.names = 1)
`

16S bacterial abundance table
`
root.asv <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/canola.asv.root.csv", header = T, row.names = 1)
`

Taxonomy for the 16S bacterial abundances:

`
root.tax <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/canola.tax.root.csv", header = T, row.names = 1)
`

Rename the columns in the subset df to something useful:

`
names(asv.int.imp1) <- c("Taxa","Group")
`
#### 2. Deal with R adding Xs to the abundance df during processing.

R doesn't like when column names start with a number. So any ASV identifiers that began with a digit now have an 'X' prepended. I found it was just easier to add X's to the taxonomy and tree files so we can link taxonomy and abundance.

Tree:

```
x <- root.tree$tip.label
x[grepl("^[[:digit:]]", x)] <- paste("X", x[grepl("^[[:digit:]]", x)], sep = "")
root.tree$tip.label <- x
```

Taxonomy:

```
x <- rownames(root.tax)
x[grepl("^[[:digit:]]", x)] <- paste("X", x[grepl("^[[:digit:]]", x)], sep = "")
rownames(root.tax) <- x
```

#### 3. Create a subsetting vector to use to prune the tree to taxa of interest

Be sure to make the data strings rather than factors:

`
asvs.28 <- as.character(asv.int.imp1$Taxa)
`

Use the vector created above to isolate the branches of interest:

```
pruned.tree <- drop.tip(root.tree, root.tree$tip.label[-match(asvs.28, root.tree$tip.label)])
tax.28 <- root.tax[pruned.tree$tip.label,]          # Use the pruned tree to subset the taxonomy file
identical(pruned.tree$tip.label, rownames(tax.28))  # Always double-check to make sure the process did what you think it did
```

#### 4. Prepare a df to use for tree-plotting

Convert relative abundances to %:

`
root.abu <- decostand(root.asv, "total")*100
`

Subset the relative abundances to the 28 we're interested in here:

`
asvs.28.abu <- root.abu[,asvs.28]
`

Calculate the taxa means, convert to df:

`
asvs.28.abu.means <- as.data.frame(colMeans(asvs.28.abu))
`

Rename the abundance column to something useful:

`
names(asvs.28.abu.means) <- "Abundance"
`

Merge the taxonomy to the 28 relative abundances:

`
merge1 <- merge(tax.28, asvs.28.abu.means, by = 0)
`

Assign the ASV IDs as the row names:

`
rownames(merge1) <- merge1$Row.names
`

Assign the ASV IDs as the row names of a grouping df we will merge below:

`
rownames(asv.int.imp1) <- asv.int.imp1$Taxa
`

Merge the taxonomy-abundance df to a grouping df:

`
merge2 <- merge(merge1, asv.int.imp1, by = 0)
`

Assign the ASV IDs as row names once again:

`
rownames(merge2) <- merge2$Row.names
`

There are two columns restuling from the merge that are identical - get rid of those:

`
merge2$Row.names <- NULL; merge2$Row.names <- NULL
`

Here is a cheeky way of sorting the merged df to match the pruned tree tip labels:

`
merge3 <- merge2[pruned.tree$tip.label,]
`

There's our double check to make sure they match:

`
identical(rownames(merge3), pruned.tree$tip.label)
`

Create columns for the abundances of the numerator and denominator groups:

``
merge3$abu.den <- 0
merge3$abu.num <- 0
``

Populate these columns:

``
merge3$abu.den[merge3$Group == "DEN"] <- merge3$Abundance[merge3$Group == "DEN"]
merge3$abu.num[merge3$Group == "NUM"] <- merge3$Abundance[merge3$Group == "NUM"]
``

Only keep the columns needed for plotting, remove relict factor levels. Then create a column that will contain the colour vectors for the numerator and denominator:

````
merge4 <- droplevels(merge3[,c(6,10:12)])
merge4$col <- NA
merge4$col[merge4$Group == "DEN"] <- "salmon"
merge4$col[merge4$Group == "NUM"] <- "skyblue"
````

Always double check the labels match:

`
identical(rownames(merge4), pruned.tree$tip.label)
`

Create unique labels for non-unique genera:

`
merge4$genus <- make.unique(as.character(merge4$genus))
`

#### 5. Now to create the phylogenetic tree to plot and add sample information

Creating a new df to work with, backfill the abundance column, and subset to only what is needed for plotting:

```
merge5 <- merge4
merge5$abundance <- merge5$abu.den + merge5$abu.num
merge5 <- merge5[,-c(3:5)]
```

Plot the tree. See https://yulab-smu.github.io/treedata-book/chapter1.html for help. I've provided two options here. 1. you can plot with extra space at the root of the tree if you'd like to add text, etc. there. 2. plot the tree normally.

Option 1: plot the tree with extra space at the root:

`
p <- ggtree(pruned.tree, layout='circular') %<+% merge5 + xlim(-0.5, NA)  # Plot the tree with extra space at the root
`

Option 2: normal plotting:

`
p <- ggtree(pruned.tree, layout='circular') %<+% merge5 + xlim(0, NA)
`

Add group information for point colour, abundance information for point size, and genus information for tip labels:

```
p$data$Group <- c(as.character(merge5$Group), rep(NA,27))
p$data$abundance <- c(merge5$abundance, rep(NA,27))
p$data$genus <- c(merge5$genus, rep(NA,27))
```

Here you can plot just the tree:

`
p
`

Obsessively double-check the plotting object and then add the genus information:

``
identical(p$data$label[1:28], rownames(merge5))
p$data$label[1:28] <- merge5$genus
``

Finally, we can plot the final tree:

``````
p2 <- p + geom_tiplab(aes(angle = angle, cex = 0.05, label = paste0('italic(', genus,')')), offset = 0.12, parse = T) +
  geom_tippoint(aes(x = p$data$x+0.07, color = as.factor(Group), size = abundance)) +
  scale_size_continuous(range = c(0.5, 5)) + 
  geom_text(aes(x = 1, y = 2, label = ""),show.legend=FALSE) +
  labs(size = "Relative abundance (%)", colour = "Group")
p2
``````

![Rplot](https://user-images.githubusercontent.com/44586553/69073421-068f9f00-09f3-11ea-9755-df3f09e1631f.jpg)

