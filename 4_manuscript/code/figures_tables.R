library(knitr)
library(kableExtra)
library(xlsx)
library(tinytex)
library(rrtable)
library(htmltools)
library(remotes)
remotes::install_github("gorkang/html2latex")
library("html2latex")
library(lme4)
library(huxtable)
library(ggplot2)
library(dplyr)
library(ggstance)

## ---- covDes
cov_description<-xlsx::read.xlsx("../0_data/covariate_list.xlsx", "variables")
  kable(cov_description, position = "h!", col.names = c("Metric",
                                     "Source",
                                     "Description"),
  align = c("l","l","r"), escape=F, caption = 'Spatial covariates included in the analysis.',"latex", booktabs=TRUE, linesep="")%>% 
  column_spec(column = 3, width = "30em") %>% 
  kable_styling(latex_options="scale_down", font_size = 9)
  # kable_styling(font_size = 6, position = "center", full_width = T)


## ---- workflow
knitr::include_graphics("../3_output/figures/chapter1_workflow.pdf", dpi = 300, auto_pdf = TRUE)

## ---- studyArea
knitr::include_graphics("../3_output/maps/studyArea_inset.png", dpi = 300, auto_pdf = TRUE)

## ---- topModels
topModels<-xlsx::read.xlsx("figures_tables/top_models.xlsx", "Sheet2")
kable(topModels, position = "h!",col.names = c("Species",
                                               "Fixed effects",
                                               "$r^{2}m$","$r^{2}c$", "AUC"),
      align = c("l","l","r", "r", "r"), escape=F, caption = 'Top models for each species.',"latex", booktabs=TRUE, linesep="", digits=2)%>% 
  column_spec(2, width = "30em")%>% 
  kable_styling(latex_options="scale_down", font_size = 9)
# kable_styling(font_size = 6, position = "center", full_width = T)


# ## ---- AUC_lag
# knitr::include_graphics("../3_output/figures/timelag_stats/sppAll_AUC_lm.png", dpi = 300, auto_pdf = TRUE)
# 
# 
# ## ---- sppAll_FACO
# knitr::include_graphics("../3_output/figures/timelag_stats/sppAll_FACO_lm.png", dpi = 300, auto_pdf = TRUE)
# 
# 
# ## ---- spp_scatter
# knitr::include_graphics("../3_output/figures/timelag_stats/spp_scatter_00to15.png", dpi = 300, auto_pdf = TRUE)
