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

#### read in data ####

transect <- read_xlsx("../01_data_transect_2019/transect_data.xlsx")

# modeling for H1 ----

# create dataset with mono pairs only
transect <- transect2 %>% filter(species_mix != "Liqu-Sapi")

# microbial respiration
bas.h1 <- lmer(formula = bas_res ~ close_distance * depth  %>% as.factor + (1|plot), data = transect)
summary(bas.h1)
performance::check_model(bas.h1)

# microbial biomass
cmic.h1 <- lmer(formula = cmic ~ close_distance * depth  %>% as.factor + (1|plot), data = transect)
summary(cmic.h1)
performance::check_model(cmic.h1)


# modeling for H2 and H3 ----

# define distance to tree 1 and tree 2
transect.out <- transect %>% 
  group_by(tsp) %>% 
  mutate(dist.T1.center = (distance_T_1 / # pos 1, 2, 3, 4, 5
                             ((min(distance_T_1)+max(distance_T_1))/2))  #pos 3   
         -1
  ) %>% 
  mutate(dist.center = ifelse(position == 1 | position == 2 | position == 3, (distance_T_1 / # pos 1, 2, 3, 4, 5
                                                        ((min(distance_T_1)+max(distance_T_1))/2)),  #pos 3
                              (distance_T_2 / # pos 1, 2, 3, 4, 5
                                 ((min(distance_T_2)+max(distance_T_2))/2))
  )) %>% 
  mutate(dist.liqu = (distance_T_1 / # pos 1, 2, 3, 4, 5
                        ((min(distance_T_1)+max(distance_T_1))/2))  #pos 3   
  ) %>% 
  mutate(dist.sapi = (distance_T_2 / # pos 1, 2, 3, 4, 5
                        ((min(distance_T_2)+max(distance_T_2))/2))  #pos 3   
  )

# check point distributions

transect.out %>% 
  ggplot() +
  geom_point(aes(x = dist.T1.center, y = bas_res, color = as.factor(depth))) 
# this shows the center as point 0 and left and right the distances to the trees 
# (left = L.formosana as T1; right = S.mukorossi as T2)


transect.out %>% 
  ggplot() +
  geom_point(aes(x = dist.center, y = bas_res, color = as.factor(depth)))
# this shows the distance of the soil-cores to the center
# center is here at 0, 0
# this is needed to calculate the overyielding based on position
# using dist.T1.center includes negative values

# calculating overyielding ----

# model microbial respiration against centered dist only for mono specific pairs 
# here we use dist.center as dist.T1.center has negative values 

#sapindus, depth = 0
mod.sap.bas.0 <- transect.out %>% 
  filter(species_mix == "Sapi-Sapi", 
         depth == 0) %>% 
  lmer(bas_res ~ dist.center + (1|plot), data = .)

summary(mod.sap.bas.0)

#sapindus, depth = 5
mod.sap.bas.5 <- transect.out %>% 
  filter(species_mix == "Sapi-Sapi", 
         depth == 5) %>% 
  lmer(bas_res ~ dist.center + (1|plot), data = .)

summary(mod.sap.bas.5)

#liquidambar, depth = 0
mod.li.bas.0 <- transect.out %>% 
  filter(species_mix == "Liqu-Liqu", 
         depth == 0) %>% 
  lmer(bas_res ~ dist.center + (1|plot), data = .)

summary(mod.li.bas.0)

#liquidambar, depth = 5
mod.li.bas.5 <- transect.out %>% 
  filter(species_mix == "Liqu-Liqu", 
         depth == 5) %>% 
  lmer(bas_res ~ dist.center + (1|plot), data = .)

summary(mod.li.bas.5)


# model microbial biomass against centered dist only for mono specific pairs

#model cmic
#sapindus, depth = 0
mod.sap.cmic.0 <- transect.out %>% 
  filter(species_mix == "Sapi-Sapi", 
         depth == 0) %>% 
  lmer(cmic ~ dist.center + (1|plot), data = .)

summary(mod.sap.cmic.0)

#sapindus, depth = 5
mod.sap.cmic.5 <- transect.out %>% 
  filter(species_mix == "Sapi-Sapi", 
         depth == 5) %>% 
  lmer(cmic ~ dist.center + (1|plot), data = .)

summary(mod.sap.cmic.5)

#liquidambar, depth = 0
mod.li.cmic.0 <- transect.out %>% 
  filter(species_mix == "Liqu-Liqu", 
         depth == 0) %>% 
  lmer(cmic ~ dist.center + (1|plot), data = .)

summary(mod.li.cmic.0)

#liquidambar, depth = 5
mod.li.cmic.5 <- transect.out %>% 
  filter(species_mix == "Liqu-Liqu", 
         depth == 5) %>% 
  lmer(cmic ~ dist.center + (1|plot), data = .)

summary(mod.li.cmic.5)


#calculating overyielding in hetero-specific pairs based on mono-specific pair predictions
trans.over.0 <- transect.out %>% 
  filter(species_mix == "Liqu-Sapi", depth == 0) %>% 
  mutate(bas_over = bas_res - (predictions(mod.sap.bas.0, 
                                           newdata = datagrid(dist.center = dist.sapi))$predicted + 
                                 predictions(mod.li.bas.0, 
                                             newdata = datagrid(dist.center = dist.liqu))$predicted)/2) %>% 
  mutate(cmic_over = cmic - (predictions(mod.sap.cmic.0, 
                                         newdata = datagrid(dist.center = dist.sapi))$predicted + 
                               predictions(mod.li.cmic.0, 
                                           newdata = datagrid(dist.center = dist.liqu))$predicted)/2)

trans.over.5 <- transect.out %>% 
  filter(species_mix == "Liqu-Sapi", depth == 5) %>% 
  mutate(bas_over = bas_res - (predictions(mod.sap.bas.5, 
                                           newdata = datagrid(dist.center = dist.sapi))$predicted + 
                                 predictions(mod.li.bas.5, 
                                             newdata = datagrid(dist.center = dist.liqu))$predicted)/2) %>% 
  mutate(cmic_over = cmic - (predictions(mod.sap.cmic.5, 
                                         newdata = datagrid(dist.center = dist.sapi))$predicted + 
                               predictions(mod.li.cmic.5, 
                                           newdata = datagrid(dist.center = dist.liqu))$predicted)/2)

over.out <- rbind(trans.over.0, trans.over.5) 

# save prediction data (overyielding)
saveRDS(object = over.out, "../01_data_transect_2019/transect_over.RDS")

# models ----

# depth effect
bas.depth <- lmer(bas_over ~ depth  %>% as.factor + (1|plot), data = over.out)
cmic.depth <- lmer(cmic_over ~ depth  %>% as.factor + (1|plot), data = over.out)

summary(bas.depth)
summary(cmic.depth)

# distance and depth effect

bas.dist.depth <- lmer(bas_over ~ dist.T1.center * depth  %>% as.factor + (1|plot), data = over.out)
cmic.dist.depth <- lmer(cmic_over ~ dist.T1.center * depth  %>% as.factor + (1|plot), data = over.out) 

summary(dist.depth)
summary(dist.depth)

# tukey test

tukey.bas <- aov(bas_over ~ depth  %>% as.factor, over.out)
summary(tukey.bas)
TukeyHSD(tukey.bas, conf.level=.95)

tukey.cmic <- aov(cmic_over ~ depth %>% as.factor, over.out)
summary(tukey.cmic)
TukeyHSD(tukey.cmic, conf.level=.95)


