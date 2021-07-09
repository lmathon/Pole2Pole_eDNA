library(tidyverse)
library(sp)
library(rgeos)
library(rgdal)
library(maptools)
library(sf)
library(raster)
library(exactextractr)


metadata_sampling <- read.csv("metadata/Metadata_eDNA_Pole2Pole_v3.csv", header = T, sep = ";", stringsAsFactors = F, na.strings=c("","NA"))

## select samples with GPS data
metadata_dist  <- metadata_sampling %>%
  filter(!is.na(longitude_start))

metadata_na  <- metadata_sampling %>%
  filter(is.na(longitude_start))

bathy <- raster("data/Bathy.tif")

crswgs84=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
pts = metadata_dist[,c("longitude_start", "latitude_start")]
pts$latitude_start <- as.numeric(pts$latitude_start)
pts$longitude_start <- as.numeric(pts$longitude_start)
pts_sp <- SpatialPoints(pts,proj4string = crswgs84)
pts_sp <- data.frame(pts_sp)
names(pts_sp) <- c("long", "lat")

# extract exact bathy
prof <- cbind(extract(bathy, pts_sp, df=T), pts_sp)

prof$station <- metadata_dist$station

#---------------------------------------------------------------------------------------------------------------------------------
# extract bathy with buffer

# ---------------------------------------------------------------------------------------------------- # 
#                     FUNCTIONS 

# df is a sf object 
# radius: in meters, the buffer radius
# bathy is the bathy file
# Outputs df with the column with the depth
function_calc_mean_depth_buffer <- function(df, bathy, radius = 1000, crs=32662, value_calculation=c("min", "max", "mean", "weighted.mean")){
  # Convert to projected data 
  df_proj <- st_transform(df, 32662)
  
  # Calculate radius
  df_proj_buffer <- st_buffer(df_proj, radius) # 1000 metre de radius
  
  # Re-transform data 
  df_proj_buffer_wgs <- st_transform(df_proj_buffer, 4326)
  
  # Test no empty geom
  test <- df_proj %>%
    mutate(is_empty = ifelse(st_is_empty(geometry) == TRUE, "yes", "no"))
  
  if(isTRUE(st_is_empty(test))){message("There is empty geometry")}
  
  # Calculate pixels
  r1 <- bathy
  p1 <- df_proj_buffer_wgs 
  
  # Add to file 
  p1$depth_value <- exact_extract(r1, p1, value_calculation)
  #p1$all_values <- raster::extract(x = r1, y = p1, FUN=mean)
  #p1$depth_mean <- sapply(p1$all_values, mean, na.rm=T)
  
  # Return the value for the column 
  return(p1$depth_value)
  
}

# ---------------------------------------------------------------------------------------------------- # 
#                     CODE 

# Set > 0 to NA - easier for calculations later 
bathy_sea <- bathy
bathy_sea[bathy_sea > -2] <- NA

# Set as sf obj
metadata_sampling_sf <- st_as_sf(metadata_dist, 
                                 coords = c("longitude_start", "latitude_start"), 
                                 crs = 4326, agr = "constant")

# Create a df with unique samples for ease of calculations 
metadata_sampling_sf_unique <- metadata_sampling_sf %>%
  distinct(code_spygen, station, geometry)

# For debug and function application
df <- metadata_sampling_sf_unique
bathy=bathy_sea
radius = 1000
crs=32662

# Add the depth column
metadata_sampling_sf_unique_filled <- metadata_sampling_sf_unique %>%
  # Min depth for different radius 
  mutate(max_depth_1_km = function_calc_mean_depth_buffer(., radius = 1000, bathy = bathy_sea, value_calculation="max"), 
         max_depth_3_km = function_calc_mean_depth_buffer(., radius = 3000, bathy = bathy_sea,value_calculation="max"),
         max_depth_5_km = function_calc_mean_depth_buffer(., radius = 5000, bathy = bathy_sea,value_calculation="max"),
         max_depth_10_km = function_calc_mean_depth_buffer(., radius = 10000, bathy = bathy_sea,value_calculation="max"))

prof_buffer <- metadata_sampling_sf_unique_filled[, c("station", "max_depth_1_km", "max_depth_3_km")]

#------------------------------------------------------------------------------------------------------------------------------
# assemble real + buffer

prof_all <- left_join(prof, prof_buffer)

# chose appropriate depth

for (i in 1:nrow(prof_all)) {
  ifelse(prof_all[i, "bathy"] <0, prof_all[i, "depth_fin"] <- prof_all[i, "bathy"], prof_all[i, "depth_fin"] <- prof_all[i, "max_depth_1_km"])
}

for (i in 1:nrow(prof_all)) {
  if(is.na(prof_all[i, "depth_fin"])){
    prof_all[i, "depth_fin"] <- prof_all[i, "max_depth_3_km"]
  }
}


# assemble with metadata

metadata <- left_join(metadata_sampling, prof_all[,c("station", "depth_fin")])
metadata <- metadata %>% distinct(code_spygen, .keep_all=T)

# if a real measured depth exist, choose it instead of the one from bathy
for (i in 1:nrow(metadata)) {
  if(!is.na(metadata[i, "depth_seafloor"])){
    metadata[i, "depth_fin"] <- metadata[i, "depth_seafloor"]
  }
}

metadata$depth_sampling <- paste("-", metadata$depth_sampling, sep="")
metadata$depth_sampling <- as.numeric(metadata$depth_sampling)
metadata$depth_sampling[is.na(metadata$depth_sampling)] <- 0
metadata$depth_fin[is.na(metadata$depth_fin)] <- 0

# if depth sampling is lower than depth_fin, change depth_fin for depth_sampling-5
for (i in 1:nrow(metadata)) {
  if(metadata[i, "depth_sampling"] < metadata[i, "depth_fin"]){
    metadata[i, "depth_fin"] <- metadata[i, "depth_sampling"] - 5
  }
}

metadata$depth_sampling[metadata$depth_sampling == 0] <- NA
metadata$depth_fin[metadata$depth_fin == 0] <- NA


write.csv(metadata, "metadata/Metadata_eDNA_Pole2Pole_v4.csv", row.names = F)

