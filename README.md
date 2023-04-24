# Pole2Pole_eDNA

Code to assemble all datasets of eDNA using MOTUs (Teleo), explore diversity patterns and produce analyses and figures. 

### Overview

This repo presents a set of functions to automatize the treatment of raw eDNA data from MiSeq sequencing and cleaned HiSeq sequencing (RapidRun).   
It is only available for outputs of ecotag and swarm.

The codes are used to perform ecological analyses on eDNA marine fish data.   

### Data

#### NGS data (swarm)

For each project / campaign / group output of NGS, you must store the data in [data/swarm/](data/swarm/).  
In each project folder there are two types of files:  
1. one \*teleo_table_motu.csv containing the occurences of each MOTU - output of *SWARM*
2. one \*teleo_ecotag_ncbi_motu.csv file containing the taxonomical assignation of each MOTU - output of *ecotag*

These files are available for samples belonging to the project, and sometimes for blank samples (empty wells on the same sequencing plate) and other samples (sequenced in the same batch as samples from the project).

For each project, the sequencing metadata file is stored in the data/swarm/\*/metadata folder, This file contains the association of each sample name (SPYGEN identification) to each sequencing run. This file is necessary to apply some quality filters.


#### Reference database

The [data/reference_database](data/reference_database/) folder contains the referencedatabase files used for the taxonomic assignment of our MOTUs.

### Metadata
 
The [metadata](metadata/) folder contains the field metadata for each sample : sample code, date, station, site, province, country, project, sample type, volume, depth of sampling, habitat and coordinates.


### Scripts

A set of functions is available in [scripts](scripts) to automatize the treatment. 
They rely on the data described in the last part. And also on some R packages, notably `taxize`.  
It is used to automatically filter only the fish species. To use `taxize`, you need an API key through a NCBI account. 
I provide here my API key. You must add it as an environmental variable in your local machine to use the functions, here is how to do it easily (to save some googling time):

```
library(rentrez)
set_entrez_key("YOUR_API_KEY") # My API key
Sys.getenv("ENTREZ_KEY") 
```

#### 01_read_cleaning

Contains the scripts to read and clean the data.
You obtain a list containing one dataframe per project and one big dataframe containing all projects cleaned with the thresholds you chose. (or you can keep the list). One of these scripts also calculate the number of reads after each step of cleaning.

#### 02_Spatial_mapping

Contains the scripts used for mapping the samples and measuring the distance to coast or to Coral triangle.

#### 03_explanatory_variables

Contains the scripts used to format explanatory variable, compute pairwise correlations and select uncorrelated variables for further analyses.


#### 04_Analysis

Contains all the scripts used for statistical analysis on the data : species and family accumulation curves, family presence proportions, MOTU distribution, diversity partitionning, dbRDA, rarity.

#### 05_Plots

Contains scripts used to produce the figures for the publication.

### outputs

Contains all the outputs of the previously mentioned analyses, sorted by folders, and the final figures used in the publication.


If you encouter any error, please contact me & fill a detailed issue on this repo. 
