# FilterZone

R Package For Detecting the Anomaly and Erroneous Zones and Alignment Filtering in Large Phylogenomic Datasets

![](/pics/header-plot.svg)

This R package is for detecting the anomaly zone and facilitating alignment and gene tree filtering of large phylogenomic datasets. 

In extreme cases of ILS, it is possible that the most common gene tree topology will not match the true species tree, a phenomenon that has been termed “the anomaly zone”. For species trees in the anomaly zone, concatenation methods can provide strong support for the most common anomalous topology (i.e. anomalous gene trees: “AGTs”) while species tree methods can recover the correct species tree as ILS is into account. However, erroneous gene trees ("EGTs") has been shown to lead to erroneous species tree topologies when gene tree estimation error is high. EGTs resulting from non-biological properties of alignments (e.g. missing data, informative sites, alignment length) produces discordant EGTs from the true species tree, and filtering based on the informativeness of the alignments can lead to more robust species tree estimation (Hutter and Duellman, in review). 

Main features of the package:
  1) Testing for anomaly zones
  2) Calculates comprehensive alignment statistics for a single or folder of alignments 
  3) Creates filtered datasets from alignment statistics (alignment length, informativeness, etc). Generates filtered datasets for:
      - Alignments (concatenated or a folder)
      - Gene trees (using the alignment statistics)
  4) Assess your filtration results through concordance factors (Minh et al. 2020) directly from R (iqtree2 program required)
  5) Plot your filtration concordance factor and anomaly zone results

This package is still in beta testing phase, and more features and expanded functionality will be added in the future. If you find any issues with something not working, or you would like features to be added, go to issues in the top menu bar and submit them. 


# Citation

Publication is in review. Hutter and Duellman, in review. 

For now, you can cite the R package by linking to this GitHub if you use it. 


# Contents

1) Installation
2) Setting up R environment
3) Detecting the anomaly zone
4) Alignment and genetree dataset filtration 
5) Testing effectiveness of filtering on anomaly zone
6) Plot filtering anomaly zone results


# 1) Installation

For the full analysis pipeline, the following programs are needed:
  1) ASTRAL-III is available on GitHub here: https://github.com/smirarab/ASTRAL
  2) IQTREE 2, specifically version 2 or higher: http://www.iqtree.org
  
Instructions for installation and testing ASTRAL-III and IQTREE 2 are included in the respective documentation.

FilterZone depends on two other R packages:
  1) AstralPlane (>=0.1): https://github.com/chutter/AstralPlane
  2) data.table (>=1.12)

The package has two additional R package dependencies, which are treated as imports (i.e. you need them installed, but library(ape) and library(stringr) not needed: 
  ape (>= 5.0)
  stringr (>= 1.4)
  
To install FilterZone, you can use the R package devtools. Here are step-by-step instructions for installation:

1) Install devtools by typing in your R console: install.packages("devtools", dependencies = T)

2) Install AstralPlane by typing in your R console: devtools::install_github("chutter/AstralPlane")

3) Devtools will ask you to install the package dependecies (ape and stringr), select "Yes". If devtools asks you to update packages, you may choose to do so. I would recommend not to install packages from source if devtools asks you. Ape is problemic from source and I could not get it to install on my machine. If devtools is still giving you trouble, you can install the dependencies with "install.packages(c("ape", "stringr"))". Then rerun Step 2 and skip package updates. 

4) Install FilterZone by typing in your R console: devtools::install_github("chutter/FilterZone"), and repeat Step 3. 

5) Devtools should finish and say the packages loaded properly. Load the packages with library(AstralPlane), library(FilterZone), and library(data.table). 

And installation should be complete. 



# 2) Setting up R environment

I have included an R script in the main repository with some examples. It is also described here in detail. 

1) first install and load the R package. Its a good idea to install new (or check) every time as this package is being updated frequently. Functions may also be modified and stop working, so check back here for updated tutorial instructions. 

```r
devtools::install_github("chutter/AstralPlane")
library(AstralPlane)

devtools::install_github("chutter/FilterZone")
library(FilterZone)
library(data.table)

```

2) You will want a character variable that includes your full path to the astral and iqtree jar files. NOTE: if you move the astral jar file, you will need to move the lib/ directory along with it, as astral depends on it. 


```r
astral.path = "/usr/local/bin/Astral-5-14/astral.5.14.2.jar"
iqtree.path = "/usr/local/bin/IQTREE/bin/iqtree2"
```

3)Setup your working directory and create if necessary

```r
work.dir = "/Test_FilterZone"
dir.create(work.dir)
setwd(work.dir)
```


# 3) Detecting the anomaly zone

1) To detect the anomaly zone, you need a species tree estimated with coalescent branch lengths, where ASTRAL-III will provide this for you. As input into ASTRAL-III, you will need gene trees estimated separately for each alignment marker in your dataset. The R package AstralPlane provides some R functions that will streamline your data analysis pipeline:
  a. alignment and gene concatenation
  b. Within gene tree filtering (collapsing nodes with low support, taxa removal)
  c. Prepare gene trees for input into ASTRAL-III
  d. A wrapper to run ASTRAL-III and import results into R
  e. Plotting and results viewing

Instructions can be found here: https://github.com/chutter/AstralPlane

Once you have a species tree through AstralPlane or on your own, you can import this species tree into R. 

2) Create a set of character variables with the path to your tree file. Also indicate your outgroups for rooting the tree. Finally, the save.name is the desired output save name. 

```r
tree.file.path = "/Trees/test-tree.tre"
outgroups = c("Species_A", "Species_B")
save.name = "test-dataset"
```

3) Next, read the tree file into R, where the read.tree function from ape works to read in trees from ASTRAL-III. Alternatively, the file path to the tree file can be input directly into the "tree" parameter in the anomalyZone function and the function will read the tree. 


```r
test.tree = ape::read.tree(tree.file.path)
anom.data = anomalyZone(tree = test.tree,
                        outgroups = outgroup.taxa)
```

Alternatively, to read the file in directly from a file path:

```r
anom.data = anomalyZone(tree = tree.file.path,
                        outgroups = outgroup.taxa)
```

Parameter explanations: 

```r
tree = tree file from ASTRAL-III read into R as a phylo object or a file path to this tree. 
outgroups = your outgroup taxa for rooting the tree
```

4) When you have collected the data for the anomaly zone across the tree, you can view the data.frame that contains the nodes and branches where the anomaly zone was detected. Additionally, you can plot the results on the phylogenetic tree: 

```r
#Estimated run time: 1 second
plot.anomalyZone(tree = uce.tree,
                 data = anom.data,
                 outgroups = outgroup.taxa,
                 save.file = NULL,
                 tip.label.size = 0.5,
                 node.label.size = 1,
                 edge.width = 3)
```

Parameter explanations: 

```r
tree = tree file from ASTRAL-III read into R as a phylo object or a file path to a tree file
data = the output data.frame from the anomalyZone function
outgroups = your outgroup taxa for rooting the tree
save.file = NULL or blank to not save a file; otherwise file name to save PDF
tip.label.size = size of the tip labels, passed to the cex function of ape::plot
node.label.size = size of the node labels, passed to the cex function of ape::nodelabels
```

![](/pics/az-example-plot.svg)


# 4) Alignment and genetree dataset filtration 

The anomaly zone occurs when there are extreme cases of ILS and the most common gene tree topology does not match the true species tree. Species tree methods are designed to take into account ILS, however, they were not designed to take gene tree estimation error into account. A recent study this package was designed for dubbed the "erroneous zone", where common gene tree estimation error can estimate an incorrect species tree while concatenation provides the correct topology (Hutter & Duellman, in review). The erroneous zone can be detected and avoided through extensive filtration of the alignments and resulting gene trees prior to phylogeny estimation using concatenation and summary species tree method. If the species tree is within the erroneous zone, after filtration of EGTs the anomaly zone will not be detected; however, under ILS AGTs are expected to occur randomly and filtration would have no impact on the detection of the anomaly zone. The results of this study are critically important to systematists, because it could provide clarity on why species tree methods provide different results than concatenation methods. 

1) To begin, you will first need a folder of alignments in phylip format and a folder of gene trees from IQTREE (other programs will probably work; if not, let me know and I can add them in). Create your working directory first (or use an existing directory). tree.files and align.files link to the gene tree files and alignments that estimated them. The names must match between the genes and alignments (except for the file extension). 


```r
work.dir = "WorkingDirectory"
tree.files = "WorkingDirectory/gene-trees"
align.files = "WorkingDirectory/alignments"
```

2) Next, you will want to select your filters to use, adjusted to the features of your dataset. Here is an example: 

```r
filter.length = c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
                  1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000,
                  2100, 2200, 2300, 2400, 2500) #number of base pairs
filter.sample = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1) #proportion samples
filter.prop.pis = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1) #proportion pis
filter.count.pis = c(10, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500) #count of pis
```

3) To obtain a table of statistics for each alignment, run the summarizeAlignments function. The inputs are the alignment directory path and the file export name. This function is general and can be used for other purposes. 


```r
#Estimated run time: 30 minutes
align.summary = summarizeAlignments(alignment.path = align.dir,
                                    file.export = "alignment_stats",
                                    alignment.type = "phylip")
```

Parameter explanations: 

```
alignment.path: path to a folder of multiple sequence alignments in phylip format
file.export: if a name is provided, the table is saved to file
overwrite: if TRUE overwrites file if it exists; FALSE the dataset is skipped
dataset.name: a unique name for your dataset. i.e. exons, introns, UCEs
alignment.type: select the format of the alignment. Phylip is avaialble for now, will be expanded in the future.
```

The summary table created by the function has the following columns: 

```
alignment name = the file name of the alignment
number_samples = number of samples or taxa in the alignment (number of rows) 
proportion_samples = the proportion of samples in the alignment (number_samples / max samples)
alignment_length = total number of base pairs long for the alignment
count_pis = number of parsimony informative sites in the alignment
proportion_pis = proportion of sites that are informative (count_pis / alignment_length)
count_missing_bp = total number of bases missing from the alignment matrix
proportion_missing_bp = proportion of bases missing from the alignment matrix (count_missing_bp / total bp)
```


4) Now that alignment statistics have been calculated, the filterSummary function can be used to obtain a quick summary of the datasets that will be generated using your selected filters (from above) and the alignment statistics. 

```r
#Estimated run time: 10 seconds
filt.summary = filterSummary(alignment.data = align.summary,
                             alignment.folder = align.dir,
                             dataset.name  = "exons",
                             file.out = "filter_summary",
                             length.filters = filter.length,
                             sample.filters = filter.sample,
                             prop.pis.filters = filter.prop.pis,
                             count.pis.filters = filter.count.pis)
```


Parameter explanations: 

```r
alignment.data: Alignment summary stats calculcated from summarizeAlignments
alignment.folder: The alignment folder from which the stats were calculated from in alignment.data
dataset.name: The name of your dataset, where all filtered datasets will be placed in this folder
file.out: if you wish to save to file, provide a file name for the summary
overwrite: whether to overwrite an existing dataset
length.filters: Your selected length filters as a vector of values for the alignment length
sample.filters: Your selected sampling fraction filters as a vector of values between 0-1
prop.pis.filters: Your selected parsimony informatives sites filter as a vector of values between 0-1
count.pis.filters: Your selected parsimony informatives sites filter as a vector of base pair counts
```

5) After alignment and filtered datasets summary statistics have been calculated, these two can be combined together to create concatenated alignments and parittion files for each filtered dataset. These concatenated alignments can be used in concatenation-based phylogenetics software (e.g. IQTREE, RAXML) to test out the impact of filtering on concatenation analyses. 

```r
#Estimated run time: 10 minutes
filterAlignments(filter.summary = filt.summary,
                 alignment.data = align.summary,
                 alignment.folder = align.dir,
                 format = "concatenated",
                 min.alignments = 5,
                 overwrite = TRUE)
```


Parameter explanations: 

```r
filter.summary: summary data file from filterSummary
alignment.data: summary data file from alignmentSummary
alignment.folder: folder of alignments to be filtered
format: save format for alignments
min.alignments: minimum number of alignments for a filtered set of data
overwrite: if TRUE overwrites file if it exists; FALSE the dataset is skipped
```

6) Prior to filtered summary species tree analysis in ASTRAL-III, the next function creates gene tree datasets for each filtered dataset, prepared for input in ASTRAL-III. Each filtered gene tree dataset is concatenated together (or placed in a folder with format = "folder") and saved to a folder called "filtered-genetrees-concatenated" for the concatenated gene trees or "filtered-genetrees-folders" for directories of gene trees for each filtered dataset. 

```r
#Estimated run time: 10 minutes
filterGeneTrees(filter.summary = filt.summary,
                alignment.data = align.summary,
                genetree.folder = tree.dir,
                format = "concatenated",
                overwrite = TRUE,
                taxa.remove = NULL,
                min.trees = 5,
                min.n.samples = 4,
                min.sample.prop = NULL,
                make.polytomy = TRUE,
                polytomy.limit = 10,
                remove.node.labels = FALSE)
```

Parameter explanations: 

```
filter.summary: summary data file from filterSummary
alignment.data: summary data file from alignmentSummary
genetree.folder: your target folder of gene trees that correspond to the alignments being filtered
format: save format for genetrees
overwrite: if TRUE overwrites file if it exists; FALSE the dataset is skipped
taxa.remove: species that you would like removed from each gene tree
min.trees: mimimum number of trees to keep filtered set
min.n.samples: the minimum number of samples to keep a gene tree
min.sample.prop: the minimum proportion of samples to keep a gene tree
make.polytomy: collapses polytomies in the gene trees
polytomy.limit: the value at which to collapse a node into a polytomy
remove.node.labels: strips trees of node labels if downstream analyses give you trouble (not recommended)
```

7) Finally, the concatenated gene tree dataset can be provided to the astralRunner function, which runs ASTRAL-III for each concatenated set of gene trees from each filtered dataset in the "filtered-genetrees-concatenated" folder. Each summary species tree is saved in the output.dir. 

```r
#Estimated run time: hours, depending on how many filters selected
AstralPlane::astralRunner(concat.genetree.folder = "filtered-genetrees-concatenated",
                          output.dir = "filtered-astral",
                          overwrite = TRUE,
                          astral.path = astral.path,
                          astral.t = 2,
                          quiet = FALSE,
                          multi.thread = TRUE,
                          memory = "8g")               
```

Parameter explanations: 

```
concat.genetree.folder: a folder of genetree files that are concatenated.
output.dir: the output directory name for the astral file
overwrite: overwrite = TRUE to overwrite existing files
astral.path: the absolute file path to the ASTRAL-III jar file
astral.t: the ASTRAL-III "t" parameter for different annotations, t = 2 is all annotation
quiet: TRUE hides the screen output from astral
multi.thread: TRUE to use Astral-MP multithreading 
memory: memory value to be passed to java. Should be in "Xg" format, X = an integer
```


# 5) Testing effectiveness of filtering on anomaly zone

1) Now that alignment and filtration statistics have been calculated and filtered ASTRAL-III trees have been estimated, this collection of data can be analyzed together. First the necessary directory paths can be put into character vectors:

```r
work.dir = "WorkingDirectory"
astral.dir = "filtered-astral"
setwd(work.dir)
```

The working directory "work.dir" is the main project folder. The "astral.dir" is the directory of filtered astral datasets saved in the "output.dir" in Step 6 above ("filtered-astral"). 


2) The goal of the next step will be to target and obtain results for specific nodes or clades. To pull results from specific clades, first create an R list object. Include all taxa in the tree from each clade in each list object position, and name the list to correspond to the desired clade names. An example is show below: 

```r
#outgroups
outgroup.taxa = c("Outgroup_genus_1", "Outgroup_genus_2")

#Set up a list of your clades of interest
taxa.set = list()
taxa.set[[1]] = c("Genus_species_1", "Genus_species_2", "Genus_species_3")
taxa.set[[2]] = c("Genus_species_1", "Genus_species_2")
taxa.set[[3]] = c("Genus_species_1", "Genus_species_2", "Genus_species_3", "Genus_species_4")
names(taxa.set) = c("node1", "node2", "node3")
```

3) Load in the alignment and filtered summary data calculated in the previous section using "read.csv". 


```r
align.summary = read.csv("alignment_summary.csv")
filt.summary = read.csv("filter_summary.csv")
```

4) Once the input data is ready, run the filterAnomalies function to collect anomaly and erroneous zone data from all the filtered datasets. This data is calculated across all filtration replicates across all branches and nodes in the tree. 

```r
#Estimated run time: 1 minute
anomaly.data = filterAnomalies(astral.directory = astral.dir,
                               outgroups = outgroup.taxa,
                               filter.data = filt.summary)
```

Parameter explanations: 

```
astral.directory: directory of filtered astral results
outgroups: outgroups to root your tree
filter.data: your master filtered dataset summary stats

```

5) Next, to obtain concordance factors data from the filtered datasets, run the filterConcordance function. The resulting table will contain the site and gene concordance factors calculated for each node across all the filtration replicates. These data.frames from filterAnomalies and filterConcordance can be used for other analyses. 

```r
#Estimated run time: 1 minute
concord.data = filterConcordance(input.dir = "concordance-factors",
                                 clade.list = taxa.set,
                                 outgroups  = outgroup.taxa)
```

Parameter explanations: 

```
input.dir: directory of concordance factor data generated from the filtered datasets
clade.list: a named list of clades of interest to test for concordance factors
outgroups: outgroups to root the tree
```


# 6) Plot filtering anomaly zone results


6) The results from the previous two functions can be summarized from the tables to find the best filtered tree, or plotted out using the plot.filterZone function. This function will plot the gCF or sCF (on the y axis) for each filtration replicate (on the x axis). In addition, the points will be colored by anomaly zone calculation presence/absence (az.colors parameter). setting dataset.name = "all" will plot all datasets together on the same plot (e.g. exons, introns, UCEs) so that the impact of filtration on concordance factors can be compared across different data types and sets of analyses. 


```r
#Estimated run time: 1 second
plot.filterZone(anomaly.zone.data = anomaly.data,
                concordance.factors.data = concord.data,
                save.plots = TRUE,
                output.dir = "Filter-Plots",
                dataset.name = "all",
                plot.gcf = TRUE,
                plot.scf = TRUE,
                az.colors = c("#7BC143", "#DE3293"),
                m.shape = c(22, 21),
                min.trees = 10)
```


Parameter explanations: 

```
anomaly.zone.data: table output from the filterAnomalies function
concordance.factors.data: table output from the filterConcordance function
save.plots: if you wish to save to file select TRUE
output.dir: the name of the output directory to save the plots if TRUE above
filter.name: the filter to plot. Options include: alignment_length, count_pis, proportion_pis, proportion_sample
dataset.name: from your main sets of data (i.e. exons, introns). All = plots of all them.
plot.gcf: should the gene concordance factor be plotted?
plot.scf: should the site concordance factor be plotted?
az.colors: colors to indicate the anomaly zone. Default: Green: presence; Purple: absence
m.shape: monophyly shape on the graph; circle = monophyletic; square paraphyletic
min.trees: minimum number of trees to keep a filtration replicate. Default: 10
```

7) This next function can summarize the filtration replicates support on a single node at a time. This function will plot the gCF or sCF (on the y axis) for each filtration replicate (on the x axis) for the given taxon group delimited at the start of the script. In addition, the points will be colored by anomaly zone calculation presence/absence (az.colors parameter). The shape (circle or square) represents whether the focal clade was monophyletic in that analysis (circle) or not (square).


```r
#Estimated run time: 1 second
plot.filterNode(anomaly.zone.data = anomaly.data,
                concordance.factors.data = concord.data,
                output.dir = "Filter-Plots",
                focal.node = "node2",
                filter.name = "alignment_length",
                dataset.name = "all",
                plot.gcf = TRUE,
                plot.scf = TRUE,
                az.colors = c("#7BC143", "#DE3293"),
                m.shape = c(22, 21),
                min.trees = 10 )
```


Parameter explanations: 

```
anomaly.zone.data: table output from the filterAnomalies function
concordance.factors.data: table output from the filterConcordance function
save.plots: if you wish to save to file select TRUE
output.dir: the name of the output directory to save the plots if TRUE above
focal.node: the node annotated from the filterConcordance function
filter.name: the filter to plot. Options include: alignment_length, count_pis, proportion_pis, proportion_sample
dataset.name: from your main sets of data (i.e. exons, introns). All = plots of all them.
plot.gcf: should the gene concordance factor be plotted?
plot.scf: should the site concordance factor be plotted?
az.colors: colors to indicate the anomaly zone. Default: Green: presence; Purple: absence
m.shape: monophyly shape on the graph; circle = monophyletic; square paraphyletic
min.trees: minimum number of trees to keep a filtration replicate. Default: 10
```




##### Coming soon: guide for many datasets at once
