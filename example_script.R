
#Load in packages
devtools::install_github("chutter/AstralPlane")
library(AstralPlane)
library(data.table)

#To do:
### Make name simpler
### Use hyphen
#### Use astral results to calcualte anomaly zone across everything
#### Collect node data from anomaly zone

######################################################################################
##### Testing the anomaly zone

tree.file = "/Users/chutter/Dropbox/Research/3_Finished-Submitted/Chan_etal_Rhacophoridae/Trees_Alignments/astral_unf_uce.tre"
outgroup.taxa = "Scaphiophryne_marmorata_CRH920"

tree.file = "/Users/chutter/Dropbox/Research/2_WIP/Hylidae/Trees/Tree_Grid/Unfiltered/Astral/Astral_uce.tre"
outgroup.taxa = c("Phyllomedusa_tomopterna_WED_55506", "Nyctimystes_infrafrenatus_SLT_771")

uce.tree = read.tree(tree.file)
anom.data = anomalyZone(tree = uce.tree,
                        outgroups = outgroup.taxa)

plot.anomalyZone(tree = uce.tree,
                 data = anom.data,
                 outgroups = outgroup.taxa,
                 save.file = NULL,
                 tip.label.size = 0.5,
                 node.label.size = 1,
                 edge.width = 3)


#### ***** ADD DUMMY NUMBER TO OUTGROUPS SO ASTRAL ROOTS PROPERLY


######################################################################################
##### Filtering demonstration

#Set up your directories
align.dir = "/Volumes/Armored/Hylidae/Alignments/all-markers_trimmed"
tree.dir = "/Users/chutter/Dropbox/Research/2_WIP/Hylidae/Trees/Gene_Trees/all-markers_trimmed"
work.dir = "/Volumes/Armored/Hylidae/Dataset-filtering"
astral.path = "/usr/local/bin/Astral-5-14/astral.5.14.2.jar"
dir.create(work.dir)
setwd(work.dir)

#Set up your filters
filter.length = c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
                  1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000,
                  2100, 2200, 2300, 2400, 2500) #number of base pairs
filter.sample = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1) #proportion samples
filter.prop.pis = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1) #proportion pis
filter.count.pis = c(10, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500) #count of pis

#Alignment summary function [20 minutes]
align.summary = summarizeAlignments(alignment.path = align.dir,
                                    file.export = "Hylidae_Alignment_stats",
                                    alignment.type = "phylip")

#Apply filters and create summary table of filters [< 1 min]
filt.summary = filterSummary(alignment.data = align.summary,
                             alignment.folder = align.dir,
                             dataset.name  = "unpartitioned",
                             file.out = "filter_summary",
                             length.filters = filter.length,
                             sample.filters = filter.sample,
                             prop.pis.filters = filter.prop.pis,
                             count.pis.filters = filter.count.pis)

#Make filtered alignments datasets [10 minutes]
filterAlignments(filter.summary = filt.summary,
                 alignment.data = align.summary,
                 alignment.folder = align.dir,
                 format = "concatenated",
                 min.alignments = 5,
                 overwrite = TRUE)

#Make filtered gene trees datasets [5 minutes]
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

#Runs astral across all filtered gene tree sets
AstralPlane::astralRunner(concat.genetree.folder = "filtered-genetrees-concatenated",
                          output.dir = "filtered-Astral",
                          overwrite = TRUE,
                          astral.path = astral.path,
                          astral.t = 2,
                          quiet = FALSE,
                          multi.thread = TRUE,
                          memory = "8g")

######################################################################################
##### Analyze data all together

#Set up your directories
align.dir = "/Volumes/Armored/Hylidae/Alignments/all-markers_trimmed"
tree.dir = "/Users/chutter/Dropbox/Research/2_WIP/Hylidae/Trees/Gene_Trees/all-markers_trimmed"
work.dir = "/Volumes/Armored/Hylidae/Dataset-filtering"
astral.dir = "filtered-Astral"
setwd(work.dir)

align.summary = read.csv("Hylidae_Alignment_stats.csv")
filt.summary = read.csv("filter_summary.csv")

#Set up a list of your clades of interest
outgroup.taxa = c("Phyllomedusa_tomopterna_WED_55506", "Nyctimystes_infrafrenatus_SLT_771")
taxa.set = list()
taxa.set[[1]] = c("Colomascirtus_staufferorum_LAC_2153", "Hyloscirtus_phyllognathus_WED_58378", "Hypsiboas_boans_WED_57877")
taxa.set[[2]] = c("Sphaenorhynchus_lacteus_WED_54090", "Scinax_garbei_LAC_1787", "Scinax_ruber_WED_56140" )
taxa.set[[3]] = c("Trachycephalus_jordani_WED_53658", "Phyrnohyas_venulosa_WED_55450", "Osteocephalus_taurinus_WED_55452")
taxa.set[[4]] = c("Scarthyla_goinorum_WED_58246" , "Lysapsus_laevis_CAS_257655", "Pseudis_paradoxus_CAS_245053")
taxa.set[[5]] = c("Dendropsophus_leucophyllatus_WED_59288", "Dendropsophus_parviceps_MZUTI_1357", "Dendropsophus_koechlini_WED_57879")
taxa.set[[6]] = c("Plectrohyla_quecchi_MVZ_251534" , "Ptychohyla_salvadorensis_EBG_518","Smilisca_phaeota_LAC_2299",
                  "Hyla_sarda_WED_54544", "Hyla_walkeri_MVZ_263408", "Dryophtes_cinerea_WED_56355")
taxa.set[[7]] = c("Acris_blanchardi_DSM_2012", "Pseudacris_triseriata_BLO_006","Hyliola_cadaverina_WED_54461")
names(taxa.set) = c("node1", "node2", "node3", "node4", "node5", "node6", "node7")


#### Function start
combined.data = filterAnomalies(astral.directory = astral.dir,
                                outgroups = outgroup.taxa,
                                filter.data = filt.summary)


######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################
##### Filtering demonstration with multi-dataset loop

#Load in packages
devtools::install_github("chutter/AstralPlane")
library(AstralPlane)
library(data.table)

#Set up your directories
align.dir = "/Volumes/Armored/Hylidae/Alignments"
tree.dir = "/Users/chutter/Dropbox/Research/2_WIP/Hylidae/Trees/Gene_Trees"
work.dir = "/Volumes/Armored/Hylidae/Dataset-filtering"
astral.path = "/usr/local/bin/Astral-5-14/astral.5.14.2.jar"
dir.create(work.dir)
setwd(work.dir)

#Set up your filters
filter.length = c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
                  1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000,
                  2100, 2200, 2300, 2400, 2500) #number of base pairs
filter.sample = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1) #proportion samples
filter.prop.pis = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1) #proportion pis
filter.count.pis = c(10, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500) #count of pis

#Gather datasets
datasets = list.dirs(path = align.dir, full.names = FALSE)
datasets = datasets[grep("_trimmed", datasets)]
datasets = datasets[grep("legacy-markers", datasets, invert = T)]

master.align.summary = data.table()
master.filt.summary = data.table()
for (i in 1:length(datasets)){

  #makes the alignment directory name from the dataset
  dataset.align = paste0(align.dir, "/", datasets[i])
  dataset.trees = paste0(tree.dir, "/", datasets[i])

  #Alignment summary function [20 minutes]
  #Do not want to save a file. Do this at the end.
  #Use datasetalign here
  align.summary = summarizeAlignments(alignment.path = dataset.align,
                                      dataset.name = datasets[i],
                                      file.export = NULL,
                                      alignment.type = "phylip")

  #Apply filters and create summary table of filters [< 1 min]
  #Do not want to save a file. Do this at the end.
  filt.summary = filterSummary(alignment.data = align.summary,
                               alignment.folder = align.dir,
                               dataset.name  = datasets[i],
                               file.out = NULL,
                               length.filters = filter.length,
                               sample.filters = filter.sample,
                               prop.pis.filters = filter.prop.pis,
                               count.pis.filters = filter.count.pis)

  #Make filtered alignments datasets [10 minutes]
  #Use datasetalign here; overwrite = FALSE
  filterAlignments(filter.summary = filt.summary,
                   alignment.data = align.summary,
                   alignment.folder = dataset.align,
                   format = "concatenated",
                   min.alignments = 5,
                   overwrite = FALSE)

  #Make filtered gene trees datasets [5 minutes]
  #Use datasetalign here; overwrite = FALSE

  #Gets list of gene trees
  gene.trees = list.files(dataset.trees)


  filterGeneTrees(filter.summary = filt.summary,
                  alignment.data = align.summary,
                  genetree.folder = dataset.trees,
                  format = "concatenated",
                  overwrite = FALSE,
                  taxa.remove = NULL,
                  min.trees = 5,
                  min.n.samples = 4,
                  min.sample.prop = NULL,
                  make.polytomy = TRUE,
                  polytomy.limit = 10,
                  remove.node.labels = FALSE)

  master.align.summary = rbind(master.align.summary, align.summary)
  master.filt.summary = rbind(master.filt.summary, filt.summary)
  print(paste0(datasets[i], " complete!"))
}#end i loop

##### Save the files here since it too awhile to generate
write.csv(master.align.summary, file = "master_alignment_summary.csv", row.names = F)
write.csv(master.filt.summary, file = "master_filter_summary.csv", row.names = F)


#Runs astral across all filtered gene tree sets
AstralPlane::astralRunner(concat.genetree.folder = "filtered-genetrees-concatenated",
                          output.dir = "filtered-astral",
                          overwrite = TRUE,
                          astral.path = astral.path,
                          astral.t = 2,
                          quiet = FALSE,
                          multi.thread = TRUE,
                          memory = "8g")


#Concordance factor analysis across all filtered datasets
concordanceRunner(alignment.dir = paste0(work.dir, "/filtered-alignments-concatenated"),
                  species.tree.dir = paste0(work.dir, "/", astral.dir),
                  genetree.dir = paste0(work.dir, "/filtered-genetrees-concatenated"),
                  output.dir = "concordance-factors",
                  overwrite = TRUE,
                  quiet = TRUE,
                  threads  = 6)


######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################
##### Analyze data all together

#Load in packages
devtools::install_github("chutter/AstralPlane")
library(AstralPlane)
library(data.table)

#Set up your directories
align.dir = "/Volumes/Armored/Hylidae/Alignments"
tree.dir = "/Users/chutter/Dropbox/Research/2_WIP/Hylidae/Trees/Gene_Trees"
work.dir = "/Volumes/Armored/Hylidae/Dataset-filtering"
astral.dir = "filtered-astral"
setwd(work.dir)

#outgroups
outgroup.taxa = c("Phyllomedusa_tomopterna_WED_55506", "Nyctimystes_infrafrenatus_SLT_771")

#Set up a list of your clades of interest
taxa.set = list()
taxa.set[[1]] = c("Colomascirtus_staufferorum_LAC_2153", "Hyloscirtus_phyllognathus_WED_58378", "Hypsiboas_boans_WED_57877")
taxa.set[[2]] = c("Sphaenorhynchus_lacteus_WED_54090", "Scinax_garbei_LAC_1787", "Scinax_ruber_WED_56140" )
taxa.set[[3]] = c("Trachycephalus_jordani_WED_53658", "Phyrnohyas_venulosa_WED_55450", "Osteocephalus_taurinus_WED_55452")
taxa.set[[4]] = c("Scarthyla_goinorum_WED_58246" , "Lysapsus_laevis_CAS_257655", "Pseudis_paradoxus_CAS_245053")
taxa.set[[5]] = c("Dendropsophus_leucophyllatus_WED_59288", "Dendropsophus_parviceps_MZUTI_1357", "Dendropsophus_koechlini_WED_57879")
taxa.set[[6]] = c("Plectrohyla_quecchi_MVZ_251534" , "Ptychohyla_salvadorensis_EBG_518","Smilisca_phaeota_LAC_2299",
                  "Hyla_sarda_WED_54544", "Hyla_walkeri_MVZ_263408", "Dryophtes_cinerea_WED_56355")
taxa.set[[7]] = c("Pseudacris_triseriata_BLO_006","Hyliola_cadaverina_WED_54461")
names(taxa.set) = c("node1", "node2", "node3", "node4", "node5", "node6", "node7")

align.summary = read.csv("master_alignment_summary.csv")
filt.summary = read.csv("master_filter_summary.csv")

#### Collect anamaly zone data from all filtered replicates
anomaly.data = filterAnomalies(astral.directory = astral.dir,
                               outgroups = outgroup.taxa,
                               filter.data = filt.summary)

#Obtains the concordance factors data for all filtered replicates
concord.data = filterConcordance(input.dir = "concordance-factors",
                                 clade.list = taxa.set,
                                 outgroups  = outgroup.taxa,
                                 all.data = TRUE)

#Plot alignment length, node 2
plot.filterCFAZ(anomaly.zone.data = anomaly.data,
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

#Plot alignment length, node 6
plot.filterCFAZ(anomaly.zone.data = anomaly.data,
                concordance.factors.data = concord.data,
                output.dir = "Filter-Plots",
                focal.node = "node6",
                filter.name = "alignment_length",
                dataset.name = "all",
                plot.gcf = TRUE,
                plot.scf = TRUE,
                az.colors = c("#7BC143", "#DE3293"),
                m.shape = c(22, 21),
                min.trees = 10 )

#Plot count PIS node 2
plot.filterCFAZ(anomaly.zone.data = anomaly.data,
                concordance.factors.data = concord.data,
                output.dir = "Filter-Plots",
                focal.node = "node2",
                filter.name = "count_pis",
                dataset.name = "all",
                plot.gcf = TRUE,
                plot.scf = TRUE,
                az.colors = c("#7BC143", "#DE3293"),
                m.shape = c(22, 21),
                min.trees = 10 )

#Plot count PIS  node 6
plot.filterCFAZ(anomaly.zone.data = anomaly.data,
                concordance.factors.data = concord.data,
                output.dir = "Filter-Plots",
                focal.node = "node6",
                filter.name = "count_pis",
                dataset.name = "all",
                plot.gcf = TRUE,
                plot.scf = TRUE,
                az.colors = c("#7BC143", "#DE3293"),
                m.shape = c(22, 21),
                min.trees = 10 )


###############################################################################
###############################################################################
###############  PLOT THE NODES OF INTEREST!           ########################
###############################################################################
###############################################################################

#### Ideas: Randomization test
#Using concordance factors, we develop an alternative test of whether the species tree
#is in the erroneous zone, predicting that gCF will be less than sCF because EGTs will be
#unable estimate the true species tree because of insufficient informative sites.
#Conversely, to support the anomaly zone, we would expect gCF and sCF to be about equal
#because ILS would occur on the underlying sites and be reflected on gene trees given accurate
#estimation.



