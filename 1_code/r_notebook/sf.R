library(rgdal)  # needed to work wih spatial data
library(sf) # work with shape files
library(raster) # work with rasters
library(tmap) # to map
library(dplyr)

#######################################
# Import layer from geodatabase
#######################################

# load gdb
fgdb <- "../downloads/MPB_AerialSurvey_2011toCurrent/Data/MPB_AerialSurvey_2011toCurrent.gdb"

# get a list of layers un the geodatabase
fc_list <- rgdal::ogrListLayers(fgdb) 

# import layer
mpb_11<-st_read("../downloads/MPB_AerialSurvey_2011toCurrent/Data/MPB_AerialSurvey_2011toCurrent.gdb", layer = "ab_0ufohn11x")
mpb_21<-st_read("../downloads/MPB_AerialSurvey_2011toCurrent/Data/MPB_AerialSurvey_2011toCurrent.gdb", layer = "ab_0ufohn21p")

#######################################
# Create a spatial object with point count locations
#######################################

# import bird data
load("0_data/manual/bird/bd_wide_2023-02-10.rData")

# select xy coordinates and station name
ss_xy<-bd_wide%>%
  ungroup()%>%
  dplyr::select(SS, x, y)%>%
  distinct()

# convert to spatial object. Use the crs argument to set the coordinate system
ss_xy<-st_as_sf(ss_xy, coords=c("x","y"), crs=4326)


#######################################
# Match point count object with polygons from the geodatabaswe
#######################################
# transform the crs of points to match that used in the geodatabse
ss_xy<-ss_xy %>% st_transform(crs=st_crs(mpb_21))

# use st_join to perform a spatial join between spatial objects
ss_xy_harvest<-st_join(ss_xy, mpb_21, join=st_intersects)


# Use dplyr to filter and explore data
ss_xy_harvest_filter<-ss_xy_harvest%>%
  dplyr::filter(!is.na(DMG_CODE))


#######################################
# map 
#######################################

# import a shape file of Alberta provincial boundary
alberta<-st_read("../Data/Canada/alberta_bc/alberta_bc.shp")%>%
  filter(PRENAME=="Alberta")

mainmap<-
  tm_shape(alberta)+
  tm_borders(col="black")+
  tm_shape(mpb_21)+
  tm_borders(col="blue")+
  # tm_shape(ss_xy)+
  # tm_symbols(shape = 4, alpha = .4, size = .02, col = "red")+
 # to color the satelite basemap
  tm_scale_bar(position=c("left", "BOTTOM"), text.color = "black", color.light="lightgrey")+
  tm_graticules(lines=FALSE)+
  tm_legend(outside=TRUE)


