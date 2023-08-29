rm(list=ls())

library("tidyverse")
library("tidyr")
library("readxl")
library("lme4")
library("lmerTest")
library("performance")
library("ggiraphExtra")
library("modelr")
library("ggpubr")
library("patchwork")
library("ggeffects")


# read in data ----

transect2 <- read_xlsx("../01_data_transect_2019/transect_data.xlsx")

# create dataset with mono pairs only
transect <- transect2 %>% filter(species_mix != "Liqu-Sapi")

# predicting data for microbial respiration ---- 

bas.h1 <- lmer(formula = bas_res ~ close_distance * depth + (1|plot), data = transect)
summary(bas.h1)
performance::check_model(bas.h1)

modelr::data_grid(transect, .model = bas.h1)

d <- seq(from = 0, to = 10, by = 0.1)
d
t <- seq(from = 0, to = 90, by = 0.1)
t

newdata <- expand.grid(d,t)
newdata$plot <- as.character("H31")

colnames(newdata) <- c("depth2", "close_distance", "plot")

# adding 2.5 to make 0-5cm = 2.5 and 5-10cm = 7.5
transect$depth2 = transect$depth + 2.5

# model with adapted depth
bas.h1.depth <- lmer(formula = bas_res ~ close_distance * depth2 +(1|plot), data = transect)
predict(bas.h1.depth)
newdata$basres <- predict(bas.h1.depth, newdata = newdata)

min(newdata$basres)
max(newdata$basres)
mean(newdata$basres) #used as midpoint in graphic colour scheme


# Fig.2 - microbial respiration ----
basal_fig <- ggplot(data = newdata, aes(x=close_distance, y=depth2, z=basres, fill = basres)) +
  geom_raster() +
  scale_fill_gradient2(low = ("purple4"), 
                       high = ("darkorange"), 
                       midpoint = mean(newdata$basres),
                       breaks = c(1.9,2.4,2.9)) +
  scale_y_continuous(trans = "reverse", 
                     breaks = c(0, 5, 10), 
                     labels = c("0", "5", "10")) +
  scale_x_continuous(breaks = c(0, 30, 60, 90), 
                     labels = c("0", "30", "60", "90")) +
  coord_equal() +
  labs(x= "Distance to tree [cm]", 
       y="Depth [cm]", 
       fill = expression(atop(Basal~respiration,"["*µl~O[2]~h^{-1}~g^{-1}~dry~soil*"]")),
       subtitle = "distance    , depth *   , distance x depth") +
  theme(panel.background = element_blank(),
        plot.title=element_text(hjust=0),
        axis.title = element_text(size = 18),
        axis.ticks=element_blank(),
        axis.text=element_text(size=16),
        legend.title=element_text(size=18),
        legend.text=element_text(size=16),
        legend.position = "right",     axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), 
                                                                                      ends = "first")),
        plot.subtitle = element_text(size = 17, 
                                     hjust = 0.5,
                                     colour = "grey30"))


basal_fig

# predicting data for microbial biomass ----

cmic.h1 <- lmer(formula = cmic ~ close_distance * depth + (1|plot), data = transect)
summary(cmic.h1)
performance::check_model(cmic.h1)

modelr::data_grid(transect, .model = cmic.h1)

d <- seq(from = 0, to = 10, by = 0.1)
d
t <- seq(from = 0, to = 90, by = 0.1)
t

newdata1 <- expand.grid(d,t)

newdata1$plot <- as.character("H31")

colnames(newdata1) <- c("depth2", "close_distance", "plot")

# adding 2.5 to make 0-5cm = 2.5 and 5-10cm = 7.5
transect$depth2 = transect$depth + 2.5

# model with adapted depth
cmic.h1.depth <- lmer(formula = cmic ~ close_distance * depth2 + (1|plot), data = transect)
predict(cmic.h1.depth)
newdata1$cmic <- predict(cmic.h1.depth, newdata = newdata1)

min(newdata1$cmic)
max(newdata1$cmic)
mean(newdata1$cmic) #used as midpoint in graphic colour scheme

# Fig. 2 - microbial biomass ----
cmic_fig <- ggplot(data = newdata1, aes(x=close_distance, y=depth2, z=cmic, fill = cmic)) +
  geom_raster() +
  scale_fill_gradient2(low = ("royalblue4"), 
                       high = ("firebrick3"), 
                       midpoint = mean(newdata1$cmic), 
                       breaks = c(600, 425, 250)) +
  scale_y_continuous(trans = "reverse", 
                     breaks = c(0, 5, 10), 
                     labels = c("0", "5", "10")) +
  scale_x_continuous(breaks = c(0, 30, 60, 90), 
                     labels = c("0", "30", "60", "90")) +
  coord_equal() +
  labs(x= "Distance to tree [cm]", 
       y= "Depth [cm]", 
       fill = expression(atop(Microbial~biomass,"["*µg~Cmic~g^{-1}~dry~soil*"]")),
       subtitle = "distance ** , depth ***, distance x depth") +
  theme(panel.background = element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_text(size=16),
        axis.title = element_text(size = 18),
        legend.title=element_text(size=18),
        legend.text=element_text(size=16),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), 
                                                       ends = "first")),
        legend.position = "right",
        plot.subtitle = element_text(size = 17,
                                     hjust = 0.5,
                                     colour = "grey30"))
cmic_fig

# arrange figures together and save ---- 

dist_depth <- cmic_fig + theme(axis.title.y = element_blank()) +
      basal_fig + theme(axis.title.y = element_blank()) + 
    plot_layout(nrow = 2, guides = "collect")

dist_out <- ggplot(data.frame(l = cmic_fig$labels$y, x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), angle = 90, size = 6) + 
  theme_void() +
  coord_cartesian(clip = "off") +
  dist_depth +
  plot_layout(
    widths = c(1, 25)
  )

ggsave("Fig_2.png", dist_out, width = 30, height = 20, units = "cm", dpi = 600)