

##### Global EDNA

############ R version 4.0.2 (2020-06-22)


###HDI
Human development Index is a synthetic measure capturing elements of life expectancy, education and wealth. We used HDI values for the year 2019 from the Human Development Indicators and Indices


###Conflicts Scrore (1946_2018)
Conflicts : For conflict value we used the yearly dataset covering 1946-2018 from the Uppsala Conflict Data Program (UCDP) (UCDP/PRIO Armed Conflict Dataset version 19.1). UCDP is the world's main provider of data on organized violence and the oldest ongoing data collection project for civil war, with a history of almost 40 years. (https://www.pcr.uu.se/research/ucdp/). We calculate a conflict score for each country which correspond to the sum of year of conflict in which the country is engaged. For each cell we assigned the conflict score of its sovereign state.


### Population 
year 2020 - 30 arc second resolution
https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals-rev11/data-download


### EEZ from marine region
v11

### No Violence/ Control of corruption 
From wgi dataset mean 2012-2018

### gravity travel time
From maire et al 2016
neartt = travel time to nearest population (minuts). If more than 500 km we assumed linear distance x speed boat (19 km/h)
gravtot = total gravity in a buffer of 500 km. (grav = pop/tt^2)

### Marine Ecosystem dependency
From Selig et al.

### Distance à la cote
From gadm v3.6
nearest neighbor  - st_nn from nngeo R package (sf_0.9-5)
distance in m - st_distance from sf  R package (nngeo_0.3.8)


#### Bathymetry
From gebco 500m


### natural resources rents
From David