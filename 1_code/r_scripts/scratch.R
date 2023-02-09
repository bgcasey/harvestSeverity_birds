
# ABMI NBR

abmi_nbr<-st_read("0_data/external/ABMI_HarvAreaSpecRegn2019_v1/ABMI_HarvestArea_SpectralRegeneration_2019_v1.shp")


bd_xy<-bd[1:3]
bd_xy<-distinct(bd_xy)

bd_xy_sf<-st_as_sf(bd_xy, coords=c("x","y"), crs=4326)

bd_xy_sf_tr<-st_transform(bd_xy_sf, st_crs(abmi_nbr))

bd_xy_cas<-st_join(bd_xy_sf_tr, abmi_nbr)

bd_filter<-bd_xy_cas%>%
  filter(YEAR>hrvYr_m & hrvYr_m>0)%>%
  na.omit()

nrow(as.data.frame(unique(bd_cov_f$SS)))


bd_abmi<-left_join(bd_filter, bird_response_variables)%>%
  mutate(dNBR_perc=nbrDstb_m/preNBR_m)





