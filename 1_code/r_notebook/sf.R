library(sf)
library(raster)
library(rgdal)
library(tmap)

# Import layer from geodatabase
# load gdb
/Volumes/Projects/downloads/MPB_AerialSurvey_2011toCurrent/Data/MPB_AerialSurvey_2011toCurrent.gdb

fgdb <- "../downloads/MPB_AerialSurvey_2011toCurrent/Data/MPB_AerialSurvey_2011toCurrent.gdb"
fc_list <- rgdal::ogrListLayers(fgdb) 

# import layer
mpb_11<-st_read("../downloads/MPB_AerialSurvey_2011toCurrent/Data/MPB_AerialSurvey_2011toCurrent.gdb", layer = "ab_0ufohn11x")
mpb_21<-st_read("../downloads/MPB_AerialSurvey_2011toCurrent/Data/MPB_AerialSurvey_2011toCurrent.gdb", layer = "ab_0ufohn21p")


# import bird data
load("0_data/manual/bird/bd_wide2023-02-10.rData")


# select xy coordinates and station name
ss_xy<-bd_wide%>%
  dplyr::select(SS, x, y)%>%
  distinct()

# convert to spatial object
ss_xy<-st_as_sf(ss_xy, coords=c("x","y"), crs=4326)

ss_xy<-ss_xy %>% st_transform(crs=st_crs(mpb_21))

ss_xy_harvest<-st_join(ss_xy, mpb_21, join=st_intersects)


## Use dplyr to filter
ss_xy_harvest_filter<-ss_xy_harvest%>%
  filter(!is.na(DMG_CODE))



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


