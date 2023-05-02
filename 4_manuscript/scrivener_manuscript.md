# Introduction

Canada's boreal region contains large intact forest spanning over 270 million hectares, providing critical breeding habitat and stopover sites for North American bird species that migrate to the tropics during winter months [@Blancher2005; @NaturalResourcesCanada2017]. Almost half of all bird species in North America rely on this region for food, shelter, and nesting sites [@wells2014boreal]. However, forestry and energy exploration have altered the composition, structure, and age distribution of boreal forests, altering bird communities and populations [@Brandt2013; @nortonSongbirdResponsePartialcut1997; @Schmiegelow1997; @Venier2014]. Therefore, understanding the impact of forestry activities on bird communities is crucial for the conservation of migratory birds in boreal forests.



Forestry companies are increasingly adopting sustainable practices, such as retention harvesting, to mitigate the effects of forestry on biodiversity and emulate natural disturbance regimes [@FedrowitzKoricheva2014; @galettoVariableRetentionHarvesting2019]. Variable retention forestry is used to promote multifunctional landscapes by maintaining pre-harvest legacy structures, such as patches of live standing trees, dead woody debris, and understory vegetation, throughout the harvest cycle [@Franklin2000; @Lindenmayer2012]. Retention accelerates the recovery of harvest blocks by increasing structural complexity relative to clear cuts, emulates natural disturbances like fire and insect outbreaks, and optimizes habitat for keystone or protected species [@galettoVariableRetentionHarvesting2019; @Lindermayer2002; @SerrouyaR.DEon2004]. By maintaining legacy structures, retention can improve the continuity of ecological processes and organisms across forest generations by maintaining habitat for species with low dispersal and/or small home ranges, enhancing ecological connectivity via residual stepping stones, accelerating stand recolonization by late successional species, moderating changes to microclimate, and maximizing niche availability through the maintenance of structural complexity [@Franklin2000; @FedrowitzKoricheva2014; @ChanMcleod2007; @BakerRead2011; @bakerMicroclimateSpaceTime2014; @BakerJordan2016; @heitheckerEdgerelatedGradientsMicroclimate2007; @TewsBrose2004]. 








Researchers have found that retention forestry can facilitate post-harvest ecological recovery [@FedrowitzKoricheva2014; @moriRetentionForestryMajor2014]. However, studies assessing the impacts of retention on bird communities have often relied on broad categorical harvest intensity metrics based on basal area or percent canopy coverage (e.g., 1-25%, 26-50%, 51-75%, 75-95%, 96-100% disturbance) to compare the short-term (<15 years) effects of harvesting. While these metrics can to predict post-harvest bird communities [@odsenBorealSongbirdsVariable2018; @priceLongtermResponseForest2020], they are often obtained from intensive fieldwork or from digital land cover maps that can be slow to update, require specialized knowledge to produce, and may not capture the full range of harvest intensities present in the landscape.






Remote sensing can provide detailed, up-to-date information on the magnitude and recovery of harvests, offering a promising alternative to categorical harvest intensity metrics, without intensive field-based vegetation surveys. Continuous measures of disturbance magnitude from remote sensing may reveal subtler relationships between forestry practices and bird communities. The Landsat program has collected a continuous 50-year archive of global land surface imagery, making it well-suited for monitoring long-term vegetation change [@wulderLandsatContinuityIssues2008; @cohenLandsatRoleEcological2004]. A growing body of image processing and change detection methods, along with the public availability of the full Landsat archive, are providing new ways to analyze optical time series data [@tewkesburyCriticalSynthesisRemotely2015; @zhuChangeDetectionUsing2017a; @gomezOpticalRemotelySensed2016]. For example, algorithms such as Landsat-based detection of Trends in Disturbance and Recovery (LandTrendr) [@Kennedy2018], Breaks For Additive Seasonal and Trend (BFAST) [@verbesseltDetectingTrendSeasonal2010], and the Continuous Change Detection and Classification (CCDC) [@zhuContinuousChangeDetection2014] can produce spectral recovery and disturbance metrics using Google Earth Engine.


In this study, we used bird point count data, acoustic monitoring tools, and spectral measures of harvest intensity to examine the impact of retention forestry on bird communities over time. Specifically, we employed an annual time series of Landsat Normalized Burn Ratio (NBR) to quantify the intensity and recovery of forest harvests and used these metrics to predict the functional and species diversity of birds within harvested areas. Our objectives were threefold: (1) quantify the influence of variable retention on bird communities along a gradient of recovery, (2) assess the efficacy of NBR as a measure of harvest intensity in predictive bird models, and (3) determine and compare the time required for bird communities to return to pre-harvest conditions following different harvest treatments.






# Methods

## Study area
Human and acoustic point counts were conducted within harvested forest areas spread across 413,161 km$^2$ of the Foothills and Boreal Forest Natural Regions of Alberta, Canada (Figure \@ref(fig:studyArea)). The Foothills Natural Region is situated at the eastern edge of the Rocky Mountains and is dominated by mixedwood forests comprising lodgepole pine (Pinus contorta), white spruce (Picea glauca), trembling aspen (Populus tremuloides), and balsam poplar (Populus balsamifera) at lower elevations. Lodgepole pine forests are typical of higher elevations [@Downing2006]. The Boreal Forest Natural Region, which spans across most of northern Alberta, comprises coniferous forests dominated by black spruce (*Picea mariana*), white spruce, and jack pine (*Pinus banksiana*), mixedwood forests with trembling aspen and balsam poplar, and shrubby black spruce fens [@Downing2006]. Despite the differences in vegetation and topography, there is a significant overlap of bird species between the foothills and boreal. In both regions, energy production and the harvesting of aspen and conifer are common. 

```{r studyArea, fig.cap= "Locations of point counts in the Boreal Central Mixedwood Natural Subregion of Northern Alberta.", out.width="100%",fig.pos="H"}
```

## Forest harvest data

We used two sources of data to identify forest harvest areas: the Common Attribute Schema for Forest Resource Inventories (CAS-FRI) [@Cumming2011a] and the Wall-to-Wall Human Footprint Inventory (HFI) [@abmi2019HFI]. CAS-FRI provides standardized 2 ha forest attributes derived from 1:10,000 to 1:40,000 aerial photography flown between 1987 and 2010 [@Cumming2011a]. HFI contains human footprint attributes interpreted from aerial and SPOT satellite imagery acquired between 1950 and 2019  [@abmi2019HFI]. We selected harvests with a range of CAS-FRI defined categorical harvest intensities (100-95, 95-75, 75-50, 50-25, and 25-0 percent). 



## Bird data 
Avian point counts were conducted in harvest areas ranging from 5 to 30 years post-harvest. We used bird detection data from databases managed by the Boreal Avian Modelling Project (BAM) [@BAM2018] and the University of Alberta's Bioacoustic Unit [@bioacousticunit2021]. Our study included data from 9700 individual point counts conducted in or within 800 meters of 1181 harvest blocks between 1995 and 2021. These point counts were conducted by human observers and autonomous recording units deployed by graduate students and field technicians. Surveys were conducted within sampling radii ranging from 50 to 150 meters and lasted between one and ten minutes. We included point counts conducted during the breeding season (May 16 to July 7) between sunrise and 10:00 h. Each location was surveyed between three and ten times approximately three days apart.  



For each point count location, we computed site-level bird species richness, Shannon diversity, and functional diversity indices based on the diet, foraging, and nesting traits of species [@EltonTraits2021] (Table \@ref(tab:responseDes)). We calculated functional divergence (FDiv), functional evenness (FEve), functional richness (FRic), and functional dispersion (FDis) using the *dbFD* function in the *FD* package in R  [@laliberteDistanceBasedFramework2010; @masonFunctionalRichnessFunctional2005; @villegerNewMultidimensionalFunctional2008; @R-FD]. Shannon diversity was calculated using the *diversity* function in the *vegan* package in R [@R-vegan]. Our analysis was limited to bird species with known breeding ranges in the Foothills and Boreal Forest Natural Regions of Alberta, Canada and that were observed at over three point count locations. 

(ref:laliberteDistanceBasedFramework2010) [@laliberteDistanceBasedFramework2010]

(ref:masonFunctionalRichnessFunctional2005) [@masonFunctionalRichnessFunctional2005]

(ref:villegerNewMultidimensionalFunctional2008) [@villegerNewMultidimensionalFunctional2008]

(ref:R-Vegan) [@R-vegan]

```{r responseDes}
```


## Model covariates

A time series' of NBR can be used to assess changes in forest structure and vegetation following disturbances. NBR is a spectral index that is calculated using near-infrared (NIR) and shortwave-infrared (SWIR) reflectance bandwidths [@keyLandscapeAssessmentGround2006]. The index can be used to differentiate between undisturbed forests and harvested or burned stands. NBR is calculated by subtracting SWIR from NIR and then dividing the result by the sum of SWIR and NIR ($NBR=\frac{NIR-SWIR}{NIR+SWIR}$)
A decrease in NBR indicates a loss of green vegetation, and the magnitude of loss can be used as a proxy for harvest intensity. When harvest areas regenerate, NBR increases with the return of green vegetation and forest structure [@hislopUsingLandsatSpectral2018; @veraverbekeEvaluatingSpectralIndices2011;@whiteConfirmationPostharvestSpectral2018]. 

To calculate harvest intensity and recovery metrics using NBR, we used methods developed by Hird et al. [-@hirdSatelliteTimeSeries2021]. First, we pre-processed harvest polygons using the *sf* package in R [@R-sf]. We repaired errors in polygon geometries, dissolved "doughnut shaped" polygons, buffered polygons by -30 m to minimize the influence of harvest edges on NBR estimates, and simplified polygons using a tolerance of 5 m. We uploaded the pre-processed harvest polygons as a shapefile asset to Google Earth Engine [@gorelickGoogleEarthEngine2017]. Next, we generated 30 m summer (June-September) composite NBR rasters from 1984 to 2021 using images from the Landsat 5 Thematic Mapper (bands 4 and 7), Landsat 7 Enhanced Thematic Mapper (bands 4 and 7), and Landsat 8 Operational Land Imager (bands 5 and 7) via Google Earth Engine's JavaScript API [@geologicalsurveyLandsat47Surface2018]. Snow, cloud, and cloud shadow pixels were masked and removed using the CFMask algorithm [@fogaCloudDetectionAlgorithm2017]. Finally, we applied the LandTrendr algorithm to NBR composites to generate the following spectral change metrics for each forest harvest area: the mean NBR value for the five years pre-harvest (NBR~pre-disturbance~), the lowest NBR post-harvest value (NBR~post-disturbance~), NBR spectral change ($\Delta$NBR=NBR~pre-disturbance~-NBR~post-disturbance~), and relative spectral change ($R\Delta$NBR=NBR~pre-disturbance~ - NBR~post-disturbance~ / ($\sqrt(|$NBR~pre-disturbance~$/1000|)$)[@MillerThode2007]

In addition to NBR recovery metrics, for each harvest we calculated the harvest area, perimeter-area ratio of harvest boundaries, and the Euclidean distance between point count locations and the nearest unharvested forest edge. These calculations were performed using the *sf* package in R [@R-sf]. The mean area of individual harvests was 115.40 ha (SD=110.23). We also calculated fractional land cover for the year of the point count, mean canopy height, and the standard deviation of canopy height within a 300 m circular buffer of point count locations using Google Earth Engine [@hermosillaLandCoverClassification2022; @lang2022high]. Model covariates and their definitions can be found in Table \@ref(tab:covDes).

(ref:MillerThode2007) [@MillerThode2007] 
(ref:lang2022high) [@lang2022high] 
(ref:hermosillaLandCoverClassification2022) [@hermosillaLandCoverClassification2022] 

```{r covDes}
```

## Analyses 

We employed mixed-effects regression models to investigate the impact of NBR harvest intensity indices on bird community metrics over time using the *glmer* package in R *lme4* [@batesFittingLinearMixedeffects2015]. We used a Poisson generalized linear mixed model with a log link to predict species richness. To predict Shannon diversity, functional divergence, functional evenness, functional richness, and functional dispersion, we employed Gamma generalized linear mixed models with a log link. To account for differences in bird detections resulting from varying sampling effort, we used the log of sampling effort (i.e., the number of survey days multiplied by point count durations) as an offset term.


We followed the same modelling process for each response variable. First, we fit a global model that included all potential predictor variables as fixed effects, and nested random effects for harvest ID and survey year. Second, we assessed the linearity of predictor-response relationships by fitting separate models using linear, quadratic, and cubic terms. Third, we addressed multicollinearity by calculating pairwise Pearson correlation coefficients and VIF scores for all predictors, and iteratively removed highly correlated predictors from the global model. We kept only predictors with low correlation (*r* < 0.5 and VIF < 3). Fourth, we used the 'dredge' function from the R package *MuMIn* to assess the performance of models containing combinations of the remaining predictors [@bartonMuMInMultimodelInference2020]. For each model we calculated pseudo-R2 as an indicator of model fit @nakagawaGeneralSimpleMethod2013]. The model with the lowest Akaike's Information Criterion (AIC) was selected as the top model [@burnhamModelSelectionMultimodel2002]. Where models had similar AIC values (differing by less than two) we chose the one with the highest pseudo-$R^2$. Finally, we calculated semi-partial $R^2$ values for predictor variables using the *r2beta* function from the *r2glmm* R package with standardized general variances.







# Results


Several fixed effects were common across the top models (Table \@ref(tab:topModels)). Time since harvest and RdNBR were applied to all the top models, and the interaction between these two variables were included in the top models for richness, functional richness, functional divergence, functional dispersion, and functional evenness. Standard deviation of canopy height and fractional water and shrub cover were common fixed effects that were top contributors to model performance across response variables. In contrast, tree species composition and harvest polygon metrics, such as the perimeter to area ratio and harvest area, were not predictive of bird communities. 
 
```{r topModels}
```





## Species diversity
Results show that species richness was negatively associated with the percentage of shrub cover, and positively associated with fractional water and herb cover, and the standard deviation of canopy height. The top model for species richness included the interaction between RdNBR and time since harvest as a fixed effect (Figure \@ref(fig:dotDiv)). The interaction led to an inverted U-shaped effect curve at levels of RdNBR <85% with maximum species richness occurring between 15 and 20 years post harvest. At high levels of harvest intensity <85%, maximum species richness occurred between 20 and 25 years post-harvest. After 20 years post harvest, the top model predicted that species richness in all harvest intensities would converge with the mean species richness of unharvested stands. The linear term for time since harvest was the strongest predictor of species richness (*b* = 0.197, SE = 0.023, *p* < 0.001, semi-partial $R^2$ = 0.037), followed by fractional cover of water  (*b* = 0.093, SE = 0.008, *p* < 0.001, semi-partial $R^2$ = 0.034), and the quadratic term for time since harvest (*b* = -0.083, SE = 0.016, *p* < 0.001, semi-partial $R^2$ = 0.014).

Shannon diversity index decreased with the percentage of exposed barren land and shrubs, and increased with the standard deviation of canopy height and the fractional cover of water. The response of Shannon diversity to time since harvest was nonlinear and followed an inverted U-shaped effect curve, with maximum Shannon diversity occurring after $\approx$ 20 years post-harvest for all harvest intensities. The top model for Shannon diversity included a quadratic term for RdNBR, resulting in a shallow U-shaped effect curve with a minimum harvest intensity of $\approx$ 45%. However, we found no significant interaction between time since harvest and RdNBR. After 20 years post-harvest, Shannon diversity in areas with harvest intensities below 75% was predicted to converge with the mean Shannon diversity of unharvested forests. The quadratic term for time since harvest was the strongest predictor of Shannon diversity (*b* = -0.046, SE = 0.015, *p* < 0.001, semi-partial $R^2$ = 0.099), followed by the proportion of water(*b* = 0.038, SE = 0.004, *p* < 0.001, semi-partial $R^2$ = 0.016), the linear term for time since harvest (*b* = 0.118, SE = 0.026, *p* < 0.001, semi-partial $R^2$ = 0.013), and the proportion of exposed barren land(*b* = , SE = , *p* < 0.001, semi-partial $R^2$ = 0.002). 


```{r dotDiv, fig.cap= "Richness and Shannon diversity estimates for all time periods and NBR derived harvest intensities with 95$\\%$  confidence intervals.", out.width="100%",fig.pos="H"}
```

## Functional diversity 

Functional richness had a negative response to shrub cover and a positive response to fractional cover of water and herbs. Our findings revealed a mild interaction between time since harvest and RdNBR, leading to an inverted U-shaped curve at harvest intensities greater than 60%. Peak functional richness observed between 15-20 years post-harvest (Figure \@ref(fig:dotFDiv). Functional richness approached the mean pre-harvest functional richness level between 15 and 20 years post-harvest across all harvest intensities. Among the predictors, fractional cover of water was the strongest predictor of functional richness (*b* = 0.123, SE = 0.016, *p* < 0.001, semi-partial $R^2$ = 0.024), followed by the linear term of time since harvest (*b* = 0.295, SE = 0.065, *p* < 0.001, semi-partial $R^2$ = 0.020), and fractional shrub cover (*b* = -0.125, SE = 0.028, *p* < 0.001, semi-partial $R^2$ = 0.014). 

Functional divergence decreased with time since harvest, shrub cover, and the distance to the nearest forested edge, and increased with RdNBR. The top model for functional divergence included the interaction between time since harvest and RdNBR. For harvest intensities below 50%, functional divergence converged with the mean functional divergence of unharvested stands after 20 years. Time since harvest was the strongest predictor of functional divergence (*b* = -0.021, SE = 0.007, *p* < 0.001, semi-partial $R^2$ = 0.038), followed by RdNBR (*b* = 0.028, SE = 0.006, *p* < 0.001, semi-partial $R^2$ = 0.038), and the distance to the nearest forested edge (*b* = -0.022, SE = 0.002, *p* < 0.001, semi-partial $R^2$ = 0.028) 



Functional dispersion was found to be negatively associated with the percentage of exposed barren land and the distance to the nearest forested edge, and positively associated with RdNBR. The top model for functional dispersion included a quadratic term for time since harvest, revealing a U-shaped response curve with the highest functional dispersion occurring between 10 and 15 years post-harvest. Functional dispersion approached the levels observed in unharvested stands between 20 and 25 years post-harvest. We found no significant interaction effect between time since harvest and RdNBR on functional dispersion. The linear term of time since harvest was the strongest predictor of functional dispersion (*b* = 0.070, SE = 0.017, *p* < 0.001, semi-partial $R^2$ = 0.017), followed by the quadratic term of time since harvest (*b* = -0.039, SE = -.010, *p* < 0.001, semi-partial $R^2$ = 0.017), and the distance to the nearest forested edge (*b* = -0.021, SE = 0.003, *p* < 0.001, semi-partial $R^2$ = 0.008) 



Functional evenness responded negatively to water, shrubs, and the standard deviation of canopy height. The top model for functional evenness included an interaction between time since harvest and RdNBR which led to an inverted U-shaped curve for harvest intensities above 50%, with maximum functional evenness observed between 15-20 years post-harvest. Functional evenness at all harvest intensities approached mean pre-harvest functional evenness between 20 and 25 years post-harvest. The interaction between the linear term of time since harvest and RdNBR was the strongest predictor of functional evenness (*b* = 0.017, SE = 0.007, *p* < 0.01, semi-partial $R^2$ = 0.010), followed by the interaction between the quadratic term of time since harvest and RdNBR (*b* = -0.008, SE = 0.003, *p* < 0.05, semi-partial $R^2$ = 0.005), and the standard deviation of canopy height(*b* = -0.007, SE = 0.004, *p* < 0.075, semi-partial $R^2$ = 0.001).
 
```{r dotFDiv, fig.cap= "Functional diversity estimates for all time periods and NBR derived harvest intensities with 95$\\%$  confidence intervals.", out.width="100%",fig.pos="H"}
```


# Discussion


Our results suggest that spectral measures of harvest intensity, post-harvest recovery time, and fractional land cover variables associated with low-lying vegetation and water are important drivers of variation in post-harvest bird communities. Our study suggests that harvest residuals can buffer the impacts of forest harvesting on birds over time. In the short term, i.e., less than five years after harvesting, we observed significant changes in bird community metrics for harvest intensities ranging from 35% to 100%. However, the rates of community recovery varied depending on the intensity and extent of harvesting. Our analysis revealed strong evidence that harvest residuals accelerated the recovery of bird species richness, functional evenness, and functional divergence. Additionally, across all harvest intensities, taxonomic richness, functional richness, functional dispersion, and functional evenness converged with levels of unharvested reference areas in less than 25 years.








Our study provides evidence that reducing harvest intensity (I.e., spectral change) can hasten the recovery of species diversity in post-harvest stands. We found that the response curves for both richness and Shannon diversity followed an inverted U-shape in relation to time since harvest, with immediate declines observed across all harvest intensities. However, we noted opposing trends between richness and Shannon diversity in response to harvest intensity. Specifically, between 1 and 20 years after harvest, increasing harvest intensity negatively impacted richness, but resulted in higher Shannon diversity. This finding may indicate increased evenness in the distribution of species abundances post-harvest [@hill1973diversity]. 

Retention promoted the recovery of species richness. Harvest intensities of less than 60% led to convergence of mean richness with un-harvested reference areas within 10 years. Shannon diversity at all harvest intensities converged with non-harvested areas within 10 years. However, between 10 and 20 years post-harvest, areas with harvest intensities >75% showed increased Shannon diversity beyond the mean Shannon diversity in un-harvested areas before eventually converging with unharvested areas 25 years post-harvest.


These findings are consistent with studies that have reported a decrease in species richness following harvests and higher levels of richness in areas with high retention compared to clear-cuts [@priceLongtermResponseForest2020; @Twedt2019; @odsenBorealSongbirdsVariable2018; @FedrowitzKoricheva2014]. This is likely due to harvest residuals increasing the horizontal and vertical vegetation heterogeneity of regenerating stands. Increases in vegetation heterogeneity can expand the niche space, leading to higher bird species richness [@TewsBrose2004; @culbert2013influence].

Our study revealed complex and contrasting trends between functional diversity indices and harvest intensity over time. Specifically, we found that functional richness exhibited a similar response to that of species richness, with an immediate decline in the size of the communities' functional trait space after harvest, followed by a rapid recovery and convergence with the mean functional richness of unharvested areas after approximately 20 years However, although our models suggested a mild negative effect of harvest intensity on functional richness, this effect was negligible and did not significantly contribute to the overall model performance. Therefore, our findings suggest that the size of the post-harvest functional trait space is similar across different harvest intensities.

 
Similar to the results of Leaver et al. [-@leaver2019response] and Edwards et al. [-@edwards2013impacts], we observed short-term declines in functional evenness (the regularity of the abundance distribution of species with different traits [@villegerNewMultidimensionalFunctional2008])  in areas with harvest intensities greater than 50%. This could be due to short-term increases in species dominance of cavity-nesting birds and birds that nest and forage in shrubs, as new recruits out-compete other species [@Schieck2006]. In areas with less than 50% harvest intensity, functional evenness did not deviate from the mean functional evenness of unharvested areas across all stages of recovery. This suggests that high amounts of harvest residuals can maintain the relative abundance of functional traits that would otherwise decline with greater harvest intensities.

Conversely, our models showed an increase in functional divergence post-harvest. Functional divergence refers to the extent to which the distribution of individual species abundances maximizes differences between functional traits [@masonFunctionalRichnessFunctional2005]. We found that harvest severity had a positive linear relationship with functional divergence, providing further evidence of the increased dominance of a few functionally distinct species at higher harvest intensities. While functional divergence decreased across all harvest intensities with time, it did not fully converge with the mean functional divergence of unharvested areas within 25 years post-harvest.

Functional dispersion, a measure of community heterogeneity that estimates the mean distance of species to the centroid of all species in the functional trait space [@laliberteDistanceBasedFramework2010], declined immediately post-harvest, followed by an inverted U-shaped response curve that peaked after 10-15 years and converged with the functional dispersion of unharvested areas after 20 years. Our models suggest that, following an initial decline, functional heterogeneity of communities increases relative to that of unharvested areas before converging after 20 years post-harvest. The increase in functional dispersion between 5 and 15 years post-harvest could be from the establishment of early seral specialists [@swanson2011forgotten].




The differential response of community metrics to harvest severity suggests that minimum retention recommendations may depend on the community indices used and the optimal timeline for recovery. Our study found that in the first 10 years after harvesting, harvest intensities of less than 50% reduced changes to species richness, functional richness, functional evenness, and functional dispersion. However, differences between harvest treatments shrank over time and for all harvest intensities, these metrics recovered to non-harvest levels within 25 years. With functional divergence and Shannon diversity, recovery took longer to reach baseline levels. For Shannon diversity, harvest intensities greater than 75% did not reach unharvested levels within 25 years. Among the metrics assessed, functional divergence was the slowest to recover, and harvest intensities greater than 50% were not predicted to reach baseline levels within 25 years post-harvest. Our findings align with a growing body of forestry research that shows that even small amounts of retention can reduce community change and accelerate post-harvest recovery [@Craig2009; @gustafssonTreeRetentionConservation2010; @Halpern2012].





Our research shows that spectral measures of disturbance derived from Landsat Time Series (LTS) data, particularly the differences in NBR between pre- and post-disturbance forests, can serve as a useful substitute for common harvest intensity metrics from ground-based stem volume and canopy cover measurements or photo-interpretation of canopy cover; metrics that are often temporally and spatially limited [@bartelsTrendsPostdisturbanceRecovery2016; @whiteConfirmationPostharvestSpectral2018]. For measuring harvest intensity, NBR has several advantages over other indices calculated from Landsat imagery such as the Normalized Difference Vegetation Index (NDVI), Tasseled Cap Greenness (TCG), and the Normalized Difference Moisture Index (NDMI) [@schultzPerformanceVegetationIndices2016]). First, the Short-Wave Infrared (SWIR) reflectance band used in NBR is sensitive to forest structure and vegetation moisture, making it well-suited for assessing the characteristics of regenerating stands [@hislopUsingLandsatSpectral2018]. Second, NBR has a slower saturation rate than other indices, allowing for more accurate measurement of long-term forest changes [@pickellForestRecoveryTrends2016]. And third, NBR can outperform other SWIR-based indices for forest disturbance monitoring [@cohenLandTrendrMultispectralEnsemble2018; @veraverbekeEvaluationPrePostfire2012]. 















While dNBR is a useful measure of harvest severity, it does not distinguish between the various silvicultural treatments present in the landscape. Our study included a variety of treatments, such as understory protection, structural retention, dispersed single-tree retention, and aggregated retention, but we did not differentiate between these treatments when estimating harvest severity. This is a common limitation in retention studies [@gustafssonTreeRetentionConservation2010; @Rosenvald2008a]. However, canopy closure and the availability of large trees are closely linked to percent retention and may compensate for the lack of differentiation between treatments [@VANDERWEL2007]. Nevertheless, distinguishing between treatments in future research could inform managers about the optimum density and spatial configuration of retained vegetation. For example, aggregated retention may create different habitat conditions than dispersed single tree retention, influencing tree mortality, edge habitat, vegetation heterogeneity, and the persistence of mature forest and early successional bird species [@CurzonKern2020; @VANDERWEL2007]. 

Retention forestry should aim to provide suitable conditions for a range of species representing different functional guilds [@gustafssonRetentionForestryMaintain2012]. Towards this, future research should target specific federally listed species at risk or those representative of important functional guilds. Furthermore, comparing natural disturbance regimes like fire with retention forestry is crucial to understanding whether the variability of post-harvest bird communities falls within natural bounds.




# Conclusion



To effectively conserve boreal birds and their habitats, it is crucial to understand the impact of forestry activities on bird populations. Our study used an annual time series of Landsat Normalized Burn Ratio to measure the intensity of forest harvesting and assess its impact on bird functional and taxonomic diversity within recovering harvested areas. While previous studies have used categorical harvest intensity metrics to assess the impacts of retention on bird communities, ours is the first to use continuous spectral measures of harvest intensity. Our findings indicate that retention forestry can mitigate the impacts of forest harvesting on bird communities. Also, including functional indices as response variables provides a more comprehensive understanding of community response compared to relying on taxonomic diversity alone [@mouillot2013functional]. Furthermore, we demonstrated the value of LandTrendr, a cloud-based change detection algorithm, as a tool for assessing harvest intensity and recovery in boreal forests. Landsat change metrics derived using LandTrendr are useful alternatives to those from traditional classified land cover maps. The findings show that methods incorporating novel remote sensing algorithms and community functional indices can reveal subtle relationships between forestry practices and bird communities over time. 


\pagebreak

# References {-}

<div id="refs"></div>

\pagebreak
