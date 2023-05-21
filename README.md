# Aspergillus fumigatus het mapping
This repository contains codes and data to analyze the heterokaryon incompatibility genes in Aspergillus fumigatus, as part of the publication of Auxier et al. 2022

The mapping contains the R script and gentoypic data used to map heterokaryon incompatiblity. The script het_mapping.R is used to genetically map the trait, and the fine_mapping.sh is a bash script used to extract the respective regions from the genomes of the two parents as well as the reference genome Af293.

The related_genomes_downloader will download the representative genomes of A. fumigatus and three related species used for Figure 5

The six extractor files will pull the regions of the coding sequence used for the phylogenetic trees, and then align and produce the ML tree using iqtree.

The fumigatus_extractor.sh script will extract the regions from all the fumigatus de novo assemblies, used for Figure 5

For the analysis of population level variation, Supplemental Figure, we used two datasets. One was from Rhodes et al. 2022, and the VCF was obtained from that manuscript. For the other, a set of German isolates, we used the scripts in the barber.data folder to download the data from Barber et al. 2020 and Barber er al. 2021
