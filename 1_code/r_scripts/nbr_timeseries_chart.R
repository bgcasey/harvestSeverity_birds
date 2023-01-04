library(readr)

nbr<-read_csv("ee-chart.csv")


library(ggplot2)

p <- ggplot(nbr, aes(x=Year)) +
  geom_line(aes(y=Fitted), color="black") + 
  # geom_line(aes(y=Original), color="red") + 
  geom_point(aes(y=Original), size=.8) +
  ylab("Normalized Burn Ratio (NBR)")+
  theme_bw()+
  theme(legend.title=element_blank())+
  theme(legend.position="bottom") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")))+
  theme(axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
ggsave("nbr_2.png", plot=last_plot(), width=6, height = 4)





library("tidyverse")
df <- nbr %>%
  select(Year, Fitted, Original) %>%
  gather(key = "variable", value = "value", -Year)
head(df)
  

gg<-ggplot(df, aes(x = Year, y = value)) + 
  geom_line(aes(linetype = variable)) + 
  # scale_color_manual(values = c("darkred", "steelblue"))+
  ylab("Normalized Burn Ratio (NBR)")+
  theme_bw()+
  theme(legend.title=element_blank())+
  theme(legend.position="bottom") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")))+
  theme(axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
ggsave("nbr.png", plot=last_plot(), width=6, height = 4)
