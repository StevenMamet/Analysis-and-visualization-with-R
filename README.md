# Visualization with R
Here I've included a number of R scripts for working with R data.

- [Venn diagrams](#venn-diagrams)
  - [1. Read in the data](#1-read-in-the-data)
  - [2. Summarize by environment (site-year combination)](#2-summarize-by-environment--site-year-combination-)
  - [3. Transpose, relabel column names, and convert to presence-absence](#3-transpose--relabel-column-names--and-convert-to-presence-absence)
  - [4. Venn diagram for the canola root microbiome.](#4-venn-diagram-for-the-canola-root-microbiome)

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
