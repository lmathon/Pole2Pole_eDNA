# Pole2Pole_eDNA

Code to assemble all datasets of eDNA using MOTUs (Teleo) and explore diversity patterns

### Overview

This repo presents a set of functions to automatize the treatment of raw eDNA data from MiSeq sequencing and cleaned HiSeq sequencing (RapidRun).   
It is only available for outputs of ecotag and swarm.   

### Data

#### NGS data

For each project / campaign / group output of NGS, you must store the data in [data/swarm/](data/All/).  
There must be two files:  
1. one Project\*.csv containing the occurences of each MOTU - output of *SWARM*
2. one Project\*.csv file containing the taxonomical assignation of each MOTU - output of *ecotag*

#### Metadata

Two types of metadata are necessary.  
1. the field metadata, from the big joint excel file located in [metadata](metadata), make sure it is updated with the most recent observations. All NGS data must have their field data in this file. 
2. the sequencing metadata, associating each sample name (SPYGEN identification) to each sequencing run. This file is necessary to apply some quality filters. One file is generated for each project. 

### Details

A set of functions is available in [scripts](scripts) to automatize the treatment. 
They rely on the data described in the last part. And also on some R packages, notably `taxize`.  
It is used to automatically filter only the fish species. To use `taxize`, you need an API key through a NCBI account. 
I provide here my API key. You must add it as an environmental variable in your local machine to use the functions, here is how to do it easily (to save some googling time):

```
library(rentrez)
set_entrez_key("e1b887b07de1764a6e68883fce0f9f69d108") # My API key
Sys.getenv("ENTREZ_KEY") 
```
