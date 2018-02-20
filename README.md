# Climatic Gradient
Reference Packages for PPlacer used to analyze metagenomic data from climatic gradient:

Salton Sea (33.328 N, 115.843 W)

Sonoran desert (33.652 N, 116.372 W)

Pinyon-juniper scrubland (33.605 N, 116.455 W)

Coastal grassland (33.737 N, 117.695 W)

Pine-oak forest (33.808 N, 116.772 W)

Subalpine forest (33.824 N, 116.755 W)

# Updated from reference database created [here](https://github.com/alex-b-chase/LRGCE/)

Major Updates include using the "Representative Genomes" from the [PATRIC](https://www.patricbrc.org/) database. Genomes were curated for corrected nomentclature and used to construct a 12,271 concatenated alignment from 21 highly-conserved, single-copy marker genes. 

These reference packages can be used to "place" short metagenomic reads onto each marker gene tree. The software uses phylogenetic inference to place each read at the "best" node in the reference phylogeny. This approach has several advantages, including more conservative taxonomic information due to genomic databases primarily consisting of human pathogenic strains.
