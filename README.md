# Aspergillus fumigatus het mapping
This repository contains codes and data to analyze the heterokaryon incompatibility genes in Aspergillus fumigatus, as part of the publication of Auxier et al. 2022

The file fine_mapping.sh is used to extract the homologus regions from both parents, for use both in gene prediction as well as for synteny analysis

The related_genomes_downloader will download the representative genomes of A. fumigatus and three related species used for Figure 5

The six extractor files will pull the regions of the coding sequence used for the phylogenetic trees, and then align and produce the ML tree using iqtree.

The fumigatus_extractor.sh script will extract the regions from all the fumigatus de novo assemblies, used for Figure 5

For the analysis of population level variation, Supplemental Figure, we used two datasets. One was from Rhodes et al. 2022, and the VCF was obtained from that manuscript. For the other, a set of German isolates, we used the scripts in the barber.data folder to download the data from Barber et al. 2020 and Barber er al. 2021
