rm(list=ls())

library("tidyverse")
library("tidyr")
library("dplyr")
library("readxl")
library("ggplot2")
library("lme4")
library("lmerTest")
library("performance")
library("ggiraphExtra")
library("modelr")
library("ggpubr")
library("marginaleffects")
library("patchwork")
library("ggeffects")

#### read in overyielding data ####

over.out <- readRDS("../01_data_transect_2019/transect_over.RDS")

# models for ggpredict 

bas.0.linear <- lmer(bas_over ~ dist.T1.center + (1|plot), data = over.out %>% filter(depth == 0))
bas.5.linear <- lmer(bas_over ~ dist.T1.center + (1|plot), data = over.out %>% filter(depth == 5))
cmic.0.linear <- lmer(cmic_over ~ dist.T1.center + (1|plot), data = over.out %>% filter(depth == 0)) 
cmic.5.linear <- lmer(cmic_over ~ dist.T1.center + (1|plot), data = over.out %>% filter(depth == 5))

## H3 mean - Microbial respiration 

(mean.bas <- over.out %>% 
  ggplot() +
  geom_hline(yintercept = 0, linetype ='dashed', col = 'grey') +
  geom_pointrange(aes(x = mean(depth), y = mean(bas_over), 
                      ymin = over.out %>% 
                        filter(depth == 0) %>% 
                        pull(bas_over) %>% 
                        mean(), 
                      ymax = over.out %>% 
                        filter(depth == 5) %>% 
                        pull(bas_over) %>% 
                        mean())) +
  geom_point(data = over.out, aes(x = mean(depth), y = mean(bas_over)), 
             shape = 15, size = 6, color = "grey", fill = "grey") +
  theme_classic() +
  theme(axis.title = element_text(size = 18),
        axis.text=element_text(size=15),
        legend.title=element_text(size=18),
        legend.text=element_text(size=15),
        legend.position = "none",
        plot.subtitle = element_text(size = 16,
                                     hjust = 0.5,
                                     colour = "grey30")) +
  ylim(-3.5,5.8) +
  scale_x_continuous(breaks = 2.417582,
                      labels = "mean") +
  labs( x = "", 
        y = "Microbial respiration \n ",
        subtitle = "mean")
)

## H3 mean - Microbial biomass

(mean.cmic <- over.out %>% 
  ggplot() +
  geom_hline(yintercept = 0, linetype ='dashed', col = 'grey') +
  geom_pointrange(aes(x = mean(depth), y = mean(cmic_over), 
                      ymin = over.out %>% 
                        filter(depth == 0) %>% 
                        pull(cmic_over) %>% 
                        mean(), 
                      ymax = over.out %>% 
                        filter(depth == 5) %>% 
                        pull(cmic_over) %>% 
                        mean())) +
  geom_point(data = over.out, aes(x = mean(depth), y = mean(cmic_over)), 
             shape = 15, size = 6, color = "grey", fill = "grey") +
  theme_classic() +
  theme(axis.title = element_text(size = 18),
        axis.text=element_text(size=15),
        legend.title=element_text(size=18),
        legend.text=element_text(size=15),
        legend.position = "none",
        plot.subtitle = element_text(size = 16,
                                      hjust = 0.5,
                                      colour = "grey30")) +
  ylim(-450, 750) +
  scale_x_continuous(breaks = 2.417582,
                     labels = "mean") +
  labs( x = "", 
        y = "Microbial biomass \n ",
        subtitle = "mean") 
)

# Tukey Test - depth effect
tukey.mod <- aov(bas_over ~ depth  %>% as.factor, over.out)
tukey.cmic <- aov(cmic_over ~ depth %>% as.factor, over.out)


summary(tukey.mod)
summary(tukey.cmic)

TukeyHSD(tukey.mod, conf.level=.95)
TukeyHSD(tukey.cmic, conf.level=.95)

mod.bas <- lmer(bas_over ~ depth %>% as.factor + (1|plot), over.out)
summary(mod.bas)
anova(mod.bas)

mod.cmic <- lmer(cmic_over ~ depth %>% as.factor + (1|plot), over.out)
summary(mod.cmic)

## H3 depth - Microbial respiration 

(depth.bas <- over.out %>% 
   ggplot() +
   geom_hline(yintercept = 0, linetype ='dashed', col = 'grey') +
   geom_line(data = data.frame(x = c(0, 5),
                               y = c(5.3, 5.3)), aes(x = x, y = y)) +
   geom_line(data = data.frame(x = c(0, 0),
                               y = c(5.3, 5.1)), aes(x = x, y = y)) +
   geom_line(data = data.frame(x = c(5, 5),
                               y = c(5.3, 5.1)), aes(x = x, y = y)) +
   geom_pointrange(data = . %>% 
                     filter(depth == 0), aes(x = depth, y = mean(bas_over), ymin = min(bas_over), ymax = max(bas_over))) +
   geom_point(data = . %>% 
                filter(depth == 0), aes(x = depth, y = mean(bas_over)), 
              size = 6, color = "#c4bbee", fill = "#c4bbee") +
   geom_pointrange(data = . %>% 
                     filter(depth == 5), aes(x = depth, y = mean(bas_over), ymin = min(bas_over), ymax = max(bas_over))) +
   geom_point(data = . %>% 
                filter(depth == 5), aes(x = depth, y = mean(bas_over)), 
              shape = 17, size = 6, color = "#332288", fill = "#332288") +
   theme_classic() +
   ylim(-3.5,5.8) +
   scale_x_continuous(limits=c(-1,6),
                      breaks = c(0,5),
                      labels = c("0 - 5", "5 - 10")) +
   labs( x = "", 
         y = "Basal respiration \n",
         subtitle = "depth") +
   theme(axis.title.y = element_blank(), 
         axis.text.y = element_blank(), 
         axis.ticks.y = element_blank(), 
         axis.line.y = element_blank(),
         axis.title = element_text(size = 18),
         axis.text=element_text(size=15),
         legend.title=element_text(size=18),
         legend.text=element_text(size=15),
         legend.position = "none",
         plot.subtitle = element_text(size = 16,
                                      hjust = 0.5,
                                      colour = "grey30"))+
   annotate("text", x = 2.5, y = 5.45, label = "*", size = 7)
)

## H3 depth - Microbial biomass

(depth.cmic <- over.out %>% 
    ggplot() +
    geom_hline(yintercept = 0, linetype ='dashed', col = 'grey') +
    geom_line(data = data.frame(x = c(0, 5),
                                y = c(730, 730)), aes(x = x, y = y)) +
    geom_line(data = data.frame(x = c(0, 0),
                                y = c(730, 707)), aes(x = x, y = y)) +
    geom_line(data = data.frame(x = c(5, 5),
                                y = c(730, 707)), aes(x = x, y = y)) +
    geom_pointrange(data = . %>% 
                      filter(depth == 0), aes(x = depth, y = mean(cmic_over), ymin = min(cmic_over), ymax = max(cmic_over))) +
    geom_point(data = . %>% 
                 filter(depth == 0), aes(x = depth, y = mean(cmic_over)), 
               size = 6, color = "#FFB000", fill = "#FFB000") +
    geom_pointrange(data = . %>% 
                      filter(depth == 5), aes(x = depth, y = mean(cmic_over), ymin = min(cmic_over), ymax = max(cmic_over))) +
    geom_point(data = . %>% 
                 filter(depth == 5), aes(x = depth, y = mean(cmic_over)), 
               shape = 17, size = 6, color = "#805800", fill = "#805800") +
    theme_classic() +
    ylim(-450, 750) +
    scale_x_continuous(limits=c(-1,6),
                       breaks = c(0,5),
                       labels = c("0 - 5", "5 - 10")) +
    labs( x = "", 
          y = "Microbial biomass \n",
          subtitle = "depth") +
    theme(axis.title.y = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          axis.line.y = element_blank(),
          axis.title = element_text(size = 18),
          axis.text=element_text(size=15),
          legend.title=element_text(size=18),
          legend.text=element_text(size=15),
          legend.position = "none",
          plot.subtitle = element_text(size = 16,
                                       hjust = 0.5,
                                       colour = "grey30"))+
    annotate("text", x = 2.5, y = 750, label = "***", size = 7)
)


## H3 distance x depth - Microbial respiration

(over.basal <- over.out %>% 
    ggplot() +
    geom_hline(yintercept = 0, linetype ='dashed', col = 'grey') +
    geom_point(aes(x = dist.T1.center, y = bas_over, 
                   shape = as.factor(depth), color = as.factor(depth))) +
    geom_line(data = ggpredict(bas.0.linear, type = "fe") %>% as.data.frame(), 
              aes(y = dist.T1.center.predicted, x = dist.T1.center.x),
              color = "#c4bbee", linewidth = 1) +
    geom_ribbon(data = ggpredict(bas.0.linear, type = "fe") %>% as.data.frame(), 
                aes(ymin =  dist.T1.center.conf.low, ymax =  dist.T1.center.conf.high, x = dist.T1.center.x),
                fill = "#c4bbee", alpha = 0.1, linewidth = 1) +
    geom_line(data = ggpredict(bas.5.linear, type = "fe") %>% as.data.frame(), 
              aes(y = dist.T1.center.predicted, x = dist.T1.center.x),
              color = "#332288", linewidth = 1) +
    geom_ribbon(data = ggpredict(bas.5.linear, type = "fe") %>% as.data.frame(), 
                aes(ymin =  dist.T1.center.conf.low, ymax =  dist.T1.center.conf.high, x = dist.T1.center.x),
                fill = "#332288", alpha = 0.075, linewidth = 1) +
    theme_classic() +
    theme(axis.title.y = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          axis.line.y = element_blank(),
          axis.title = element_text(size = 18),
          axis.text=element_text(size=15),
          legend.title=element_text(size=18),
          legend.text=element_text(size=15),
          legend.position = "right",
          plot.subtitle = element_text(size = 16, 
                                       hjust = 0.5,
                                       colour = "grey30")) +
    ylim(-3,5) +
    labs( x = "", 
          y = "Basal respiration overyielding \n",
          subtitle = "distance x depth") +
    scale_x_continuous(breaks = c(-0.8, 0, 0.8), 
                       labels = c(expression(atop("tree", paste(italic("L. formosana")))),
                                  "center",
                                  expression(atop("tree", paste(italic("S. mukorossi")))))) +
    scale_color_manual(name = "Depth [cm]", 
                       values = c("#c4bbee", "#332288"), 
                       labels = c("0 - 5","5 - 10")) +
    scale_shape_manual(name = "Depth [cm]", 
                       values = c(16, 17), 
                       labels = c("0 - 5","5 - 10"))
)

## H3 distance x depth - Microbial biomass

(over.cmic <- over.out %>% 
    ggplot() +
    geom_hline(yintercept = 0, linetype ='dashed', col = 'grey') +
    geom_point(aes(x = dist.T1.center, y = cmic_over, shape = as.factor(depth), color = as.factor(depth))) +
    geom_line(data = ggpredict(cmic.0.linear, type = "fe") %>% as.data.frame(), 
              aes(y = dist.T1.center.predicted, x = dist.T1.center.x),
              color = "#FFB000", linewidth = 1) +
    geom_ribbon(data = ggpredict(cmic.0.linear, type = "fe") %>% as.data.frame(), 
                aes(ymin =  dist.T1.center.conf.low, ymax =  dist.T1.center.conf.high, x = dist.T1.center.x),
                fill = "#FFB000", alpha = 0.1, linewidth = 1) +
    geom_line(data = ggpredict(cmic.5.linear, type = "fe") %>% as.data.frame(), aes(y = dist.T1.center.predicted, x = dist.T1.center.x),
              color = "#805800", linewidth = 1) +
    geom_ribbon(data = ggpredict(cmic.5.linear, type = "fe") %>% as.data.frame(), 
                aes(ymin =  dist.T1.center.conf.low, ymax =  dist.T1.center.conf.high, x = dist.T1.center.x),
                fill = "#805800", alpha = 0.1, linewidth = 1) +
    theme_classic() +
    theme(axis.title.y = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          axis.line.y = element_blank(),
          axis.title = element_text(size = 18),
          axis.text=element_text(size=15),
          legend.title=element_text(size=18),
          legend.text=element_text(size=15),
          legend.position = "right",
          plot.subtitle = element_text(size = 16, 
                                       hjust = 0.5,
                                       colour = "grey30")) +
    ylim(-450, 750) +
    labs( x = "", 
          y = "Microbial biomass overyielding \n",
          subtitle = "distance x depth") +
    scale_x_continuous(breaks = c(-0.8, 0, 0.8), 
                       labels = c(expression(atop("tree", paste(italic("L. formosana")))),
                                  "center",
                                  expression(atop("tree", paste(italic("S. mukorossi")))))) +
    scale_color_manual(name = "Depth [cm]", 
                       values = c("#FFB000", "#805800"), 
                       labels = c("0 - 5","5 - 10")) +
    scale_shape_manual(name = "Depth [cm]", 
                       values = c(16, 17), 
                       labels = c("0 - 5","5 - 10"))
)


## putting plots together 

library(patchwork)
basal.resp <- mean.bas + theme(axis.title.y = element_blank()) + depth.bas + over.basal + 
  plot_layout(widths = c(0.25, 0.5, 1))

(basal.resp.out <- ggplot(data.frame(l = "Microbial respiration", x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), angle = 90, size = 8) + 
  theme_void() +
  coord_cartesian(clip = "off") + 
    ggplot(data.frame(l = "underyielding        overyielding", x = 1, y = 1)) +
    geom_text(aes(x, y, label = l), angle = 90, size = 6, color = "grey30") +
    theme_void() +
    coord_cartesian(clip = "off") +
    basal.resp +
  plot_layout(widths = c(1,1.1,25))  )
  


micro.biomass <- mean.cmic + theme(axis.title.y = element_blank()) + depth.cmic + over.cmic + 
  plot_layout(widths = c(0.25, 0.5, 1))

(micro.biomass.out <- ggplot(data.frame(l = "Microbial biomass", x = 1, y = 1)) +
    geom_text(aes(x, y, label = l), angle = 90, size = 8) + 
    theme_void() +
    coord_cartesian(clip = "off") + 
    ggplot(data.frame(l = "underyielding        overyielding", x = 1, y = 1)) +
    geom_text(aes(x, y, label = l), angle = 90, size = 6, color = "grey30") +
    theme_void() +
    coord_cartesian(clip = "off") +
    micro.biomass +
    plot_layout(widths = c(1,1.1,25)))



(overyielding <- micro.biomass.out / basal.resp.out +
    plot_annotation(tag_levels = list(c("", "", "A", "B", "C", "", "", "D", "E", "F"))))

ggsave("../03_output/Figure_3.png", overyielding, width = 23, height = 20, units = "cm", dpi = 600)