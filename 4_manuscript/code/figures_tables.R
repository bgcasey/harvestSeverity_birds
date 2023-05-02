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
library(jtools)
library(readxl)


## ---- covDes
# cov_description<-xlsx::read.xlsx("../3_output/tables/covariate_list.xlsx", "variables")
# kable(cov_description, position = "h!", col.names = c("Metric",
#                                                       "Source",
#                                                       "Description"),
#       align = c("l","l","r"), escape=F, caption = 'Spatial covariates included in the analysis.',"latex", booktabs=TRUE, linesep=c("", "", "", "", "\\addlinespace"))%>% 
#   # column_spec(column = 3, width = "50em") %>% 
#   kable_styling(latex_options="scale_down", font_size = 9)
# kable_styling(font_size = 6, position = "center", full_width = T)
cov_description<-readxl::read_xlsx("../0_data/covariate_list.xlsx", "variables")
# cov_description<-xlsx::read.xlsx("../0_data/covariate_list_1.xlsx", "variables")
kable(cov_description, position = "h!",
      align = c("l","l"), escape=F, caption = 'Spatial covariates included in the analysis.',"latex", booktabs=TRUE, linesep="")%>%
  kable_styling(latex_options="scale_down", font_size = 8)
# kable_styling(font_size = 6, position = "center", full_width = T)


## ---- responseDes
res_description<-readxl::read_xlsx("../0_data/response_variable_list.xlsx", "Sheet1")
# cov_description<-xlsx::read.xlsx("../0_data/covariate_list_1.xlsx", "variables")
kable(res_description, position = "h!",
      align = c("l","l"), escape=F, caption = 'Response variables included in the analysis.',"latex", booktabs=TRUE, linesep="")%>%
  column_spec(2, width = "30em")%>% 
  kable_styling(latex_options="scale_down", font_size = 8)
# kable_styling(font_size = 6, position = "center", full_width = T)



# ## ---- workflow
# knitr::include_graphics("../3_output/figures/chapter1_workflow.pdf", dpi = 300, auto_pdf = TRUE)

## ---- studyArea
knitr::include_graphics("../3_output/maps/studyArea_inset.png", dpi = 300, auto_pdf = TRUE)

## ---- dotFDiv
knitr::include_graphics("../3_output/figures/time_nbrDist_dot_fdiv_2023-04-25.png")

## ---- dotDiv
knitr::include_graphics("../3_output/figures/time_nbrDist_dot_div_2023-04-25.png", dpi = 300)



## ---- mPerf
load("../3_output/tables/models_performance_2023_05-01.rData")
m_perf<-m_perf%>%
  rename(R2c=R2_conditional, R2m=R2_marginal)%>%
  mutate(across(where(is.numeric), round, 2))
kable(m_perf, position = "h!",
      # align = c("l","l"), 
      escape=F, caption = 'Performance of top models.',"latex", booktabs=TRUE, linesep="")%>%
  # column_spec(2, width = "30em")%>% 
  kable_styling(latex_options="scale_down", font_size = 8)
# kable_styling(font_size = 6, position = "center", full_width = T)



## ---- topModels
topModels<-readxl::read_xlsx("../3_output/tables/top_models.xlsx", "Sheet2")%>%
  left_join(m_perf)%>%
  rename(R2c=R2_conditional, R2m=R2_marginal)%>%
  mutate(across(where(is.numeric), round, 2))%>%
  dplyr::select("Response variable", "Fixed effects", "AIC", "R2c", "R2m")%>%
  rename("$R^{2}m$"=R2m,"$R^{2}c$"=R2c)%>%
  mutate(across(where(is.numeric), round, 2))
kable(topModels, position = "h!",
      # align = c("l","l"), 
      escape=F, caption = 'The fixed effects and summary statistics for top models for each response variable.',"latex", booktabs=TRUE, linesep="")%>%
  column_spec(2, width = "30em")%>%
  kable_styling(latex_options="scale_down", font_size = 8)
# kable_styling(font_size = 6, position = "center", full_width = T)

## ---- modelStats
knitr::include_graphics("../3_output/tables/models_stats_2023-04-26.png", dpi = 300, auto_pdf = TRUE)







# ## ---- topModels
# topModels<-xlsx::read.xlsx("figures_tables/top_models.xlsx", "Sheet2")
# kable(topModels, position = "h!",col.names = c("Species",
#                                                "Fixed effects",
#                                                "$r^{2}m$","$r^{2}c$", "AUC"),
#       align = c("l","l","r", "r", "r"), escape=F, caption = 'Top models for each species.',"latex", booktabs=TRUE, linesep="", digits=2)%>% 
#   column_spec(2, width = "30em")%>% 
#   kable_styling(latex_options="scale_down", font_size = 9)
# # kable_styling(font_size = 6, position = "center", full_width = T)


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
