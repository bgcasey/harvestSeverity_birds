load("0_data/manual/for_models/bd_cov_processed_scaled_filtered_2023-03-20.rData")
load("/Volumes/Projects/harvestSeverity_birds/0_data/manual/bird/bd_wide_2023-03-20.rData")

#number of point count events
d<-semi_join(bd_wide, bd_cov_processed_scaled_filtered)
d1<-d%>%dplyr::select(SS, DATE, TIME)%>%distinct()
nrow(d1)
#9700 point counts
nrow(as.data.frame(unique(bd_cov_processed_scaled_filtered$SS)))
# 1739 - SS
nrow(as.data.frame(unique(dplyr::select(ungroup(bd_cov_processed_scaled_filtered), SS, survey_year))))
# 3219 - SS and years
nrow(as.data.frame(unique(dplyr::select(ungroup(bd_cov_processed_scaled_filtered), HFI_ID))))
# 1181 harvest blocks

# 641 harvest blocks

# Harvest area
d<-bd_cov_processed_scaled_filtered
sc_1 <- attr(d$area, 'scaled:scale')
ce_1 <- attr(d$area, 'scaled:center')

d<-d%>%
  mutate(area_us=(area*sc_1+ce_1)*0.0001) #scale and convert to hectares

summary(d$area_us)
mean(d$area_us)
#115.4039
sd(d$area_us)
#110.2337


