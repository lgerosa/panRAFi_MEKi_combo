library(tidyverse)
library(GGally)
library(ggrepel)
library(ggforce)
library(ggpubr)
library(mrgsolve)
library(here)

theme.to.use <-
  theme_bw() +
  theme(strip.text.x = element_text(size=12),
        text=element_text(family="Arial"),
        strip.background = element_rect(fill="#c2c2c2"),
        axis.title.y=element_text(size=12, margin=margin(0, 10, 0, 0), face="bold"),
        axis.title.x=element_text(size=12, margin=margin(10, 0, 0, 0), face="bold"),
        axis.text.x = element_text(size=12, face="bold"),
        axis.text.y = element_text(size=12, face="bold"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        plot.title = element_text(size=12, face = "bold"),
        plot.subtitle = element_text(size=15),
        plot.caption = element_text(size=15),
        legend.position="bottom"
        )
theme_set(theme.to.use)


tgi_data <- read_csv(here("data","PK_variability" ,"clinical_tgi_data.csv")) %>% 
  mutate(cobi.level = factor(cobi.level, levels = c("Belva Mono",
                                                    "Belva + Cetux","Belva + Cobi 20mg QOD/TIW",
                                                    "Belva + Cobi 20mg QD","Belva + Cobi 40mg QD")))

set.seed(12345)
colors <- c("#C183C9","#B4D39D")
scale_color_cohort <- function(...){
  ggplot2:::manual_scale(
    'color', 
    values = setNames(colors, unique(tgi_data$mutation)), 
    ...
  )
}

scale_fill_cohort <- function(...){
  ggplot2:::manual_scale(
    'fill', 
    values = setNames(colors, unique(tgi_data$mutation)), 
    ...
  )
}


p1<- tgi_data %>% 
  filter(mutation == "BRAF+") %>% 
  filter(STUDYID %in% c("HM-RAFI-103","HM-RAFI-102")) %>% 
  ggplot(aes(cobi.level, LNKG)) + 
  geom_boxplot(aes(fill = cobi.level),alpha = 0.4, outlier.shape = NA, show.legend = F) +
  geom_point(aes(color = cobi.level), position = position_jitterdodge(0.1), show.legend = F) + 
  ggtitle("a) Tumor Growth Rate BRAF+") + 
  #scale_color_cohort() +
  #scale_fill_cohort() + 
  xlab("") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 0)) + 
  ylab("log(Growth Rate)") 

p2<- tgi_data %>% 
  filter(mutation == "NRAS+") %>% 
  filter(STUDYID %in% c("HM-RAFI-103","HM-RAFI-102")) %>% 
  ggplot(aes(cobi.level, LNKG)) + 
  geom_boxplot(aes(fill = cobi.level),alpha = 0.4, outlier.shape = NA, show.legend = F) +
  geom_point(aes(color = cobi.level), position = position_jitterdodge(0.1), show.legend = F) + 
  ggtitle("b) Tumor Growth Rate NRAS+") + 
  #  scale_color_cohort() +
  #scale_fill_cohort() + 
  xlab("") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 0)) + 
  ylab("log(Growth Rate)") 

p3 <- tgi_data %>% 
  filter(mutation == "BRAF+") %>% 
  filter(STUDYID %in% c("HM-RAFI-103","HM-RAFI-102")) %>% 
  ggplot(aes(cobi.level, LNKS)) + 
  geom_boxplot(aes(fill = cobi.level),alpha = 0.4, outlier.shape = NA, show.legend = F) +
  geom_point(aes(color = cobi.level), position = position_jitterdodge(0.1), show.legend = F) + 
  theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1, size = 12)) + 
  #  scale_color_cohort() +
  #scale_fill_cohort() + 
  ggtitle("c) Tumor Shrinkage Rate BRAF+") + 
  xlab("") + 
  ylab("log(Shrinkage Rate)") 


p4 <- tgi_data %>% 
  filter(mutation == "NRAS+") %>% 
  filter(STUDYID %in% c("HM-RAFI-103","HM-RAFI-102")) %>% 
  ggplot(aes(cobi.level, LNKS)) + 
  geom_boxplot(aes(fill = cobi.level),alpha = 0.4, outlier.shape = NA, show.legend = F) +
  geom_point(aes(color = cobi.level), position = position_jitterdodge(0.1), show.legend = F) + 
  theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1, size = 12)) + 
  #scale_color_cohort() +
  #scale_fill_cohort() + 
  ggtitle("d) Tumor Shrinkage Rate NRAS+") + 
  xlab("") + 
  ylab("Log(Shrinkage Rate)") 


psum_alt <- ggpubr::ggarrange(p1, p2,p3,p4, nrow = 2, ncol = 2)

psum_alt


p1<- tgi_data %>% 
  filter(mutation == "BRAF+") %>% 
  filter(STUDYID %in% c("HM-RAFI-103","HM-RAFI-102")) %>% 
  ggplot(aes(cobi.level, LNKG)) + 
  geom_boxplot(aes(fill = cobi.level),alpha = 0.4, outlier.shape = NA, show.legend = F) +
  geom_point(aes(color = cobi.level), position = position_jitterdodge(0.1), show.legend = F) + 
  ggtitle("a) Tumor Growth Rate BRAF+") + 
  #scale_color_cohort() +
  #scale_fill_cohort() + 
  xlab("") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 0)) + 
  ylab("log(Growth Rate)") 

p2 <- tgi_data %>% 
  filter(mutation == "BRAF+") %>% 
  filter(STUDYID %in% c("HM-RAFI-103","HM-RAFI-102")) %>% 
  ggplot(aes(cobi.level, LNKS)) + 
  geom_boxplot(aes(fill = cobi.level),alpha = 0.4, outlier.shape = NA, show.legend = F) +
  geom_point(aes(color = cobi.level), position = position_jitterdodge(0.1), show.legend = F) + 
  theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1, size = 0)) + 
  #  scale_color_cohort() +
  #scale_fill_cohort() + 
  ggtitle("b) Tumor Shrinkage Rate BRAF+") + 
  xlab("") + 
  ylab("log(Shrinkage Rate)") 


p3<- tgi_data %>% 
  filter(mutation == "NRAS+") %>% 
  filter(STUDYID %in% c("HM-RAFI-103","HM-RAFI-102")) %>% 
  ggplot(aes(cobi.level, LNKG)) + 
  geom_boxplot(aes(fill = cobi.level),alpha = 0.4, outlier.shape = NA, show.legend = F) +
  geom_point(aes(color = cobi.level), position = position_jitterdodge(0.1), show.legend = F) + 
  ggtitle("c) Tumor Growth Rate NRAS+") + 
  #  scale_color_cohort() +
  #scale_fill_cohort() + 
  xlab("") + 
  theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1, size = 12)) + 
  ylab("log(Growth Rate)") 

p4 <- tgi_data %>% 
  filter(mutation == "NRAS+") %>% 
  filter(STUDYID %in% c("HM-RAFI-103","HM-RAFI-102")) %>% 
  ggplot(aes(cobi.level, LNKS)) + 
  geom_boxplot(aes(fill = cobi.level),alpha = 0.4, outlier.shape = NA, show.legend = F) +
  geom_point(aes(color = cobi.level), position = position_jitterdodge(0.1), show.legend = F) + 
  theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1, size = 12)) + 
  #scale_color_cohort() +
  #scale_fill_cohort() + 
  ggtitle("d) Tumor Shrinkage Rate NRAS+") + 
  xlab("") + 
  ylab("Log(Shrinkage Rate)") 


psum <- ggpubr::ggarrange(p1, p2,p3,p4, nrow = 2, ncol = 2)

psum



### TGI Simulation

mod <- mread(here('data','pk_models', "tgi_model.cpp"))


tmp_summary_alt <- tgi_data %>% 
  filter(STUDYID %in% c("HM-RAFI-103","HM-RAFI-102")) %>% 
  group_by(mutation, cobi.level) %>% 
  summarise(TVKG = exp(mean(LNKG)),
            TVKS = exp(mean(LNKS)),
            TVBLS = 57) %>% 
  ungroup()

tmp_summary_sim <- tmp_summary_alt %>% 
  group_by(mutation, cobi.level) %>% 
  mutate(ID = cur_group_id()) %>% ungroup() %>% 
  mutate(TVBSL = 50 ) %>% 
  expand_grid(TIME = 0:52) %>% 
  arrange(ID, TIME)

set.seed(12345)
tgi_sim_out <- mod %>% 
  data_set(tmp_summary_sim) %>% 
  carry_out(ID,KG, KS, BSL) %>% 
  zero_re() %>% 
  mrgsim(end = 52, delta = 0.5, tad = TRUE, recover = c("mutation", "cobi.level")) %>% 
  as_tibble()



p5 <- tgi_sim_out %>% 
  ggplot(aes(TIME, IPRED)) + 
  geom_line(aes(color = cobi.level), size =1) + 
  facet_wrap(~mutation, nrow = 2) + 
  ylab("Simulatated Tumor Size (mm)") +
  ggtitle("e) Simulated Tumor Size") + 
  xlab("Time (wks)") +
  guides(color = guide_legend(nrow = 3)) + 
  labs(color = "")

p5

pfinal <- ggpubr::ggarrange(psum_alt, p5, nrow = 2)


ggsave(filename = here("figures", "clinical_tgi", "figure_7_other.png"),
       pfinal,
       width = 10,
       height = 14,
       units = "in",
       device = "png",
       dpi = 300
        )
