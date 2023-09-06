##                        ##
##    Rhizome Depletion   ##
##    Nonlinear Models    ##
##                        ##

## In this script we  fit and test nonlinear models
## to all the variables I measured, linear models may be used.
## Then I'll test treatment effects on parameters. 
## models will also be saved separately. 
## I'll load specific dataset for each variable, not sure if
## I should trust the all variable integration.

## Off we go!
# rm(list=ls(all=T))


## Load Packages ###

require(dplyr)
require(ggplot2)
require(emmeans)
require(cowplot)
require(kableExtra)
require(nlme)
require(plotly)
require(ggrepel)
require(lme4)


myTheme <- {theme(panel.background = element_rect(fill = NA, color = "black"),
                  panel.grid = element_blank(),
                  panel.spacing.x = unit(0, "cm"),
                  strip.background = element_rect(fill = NA, color = "black"),
                  strip.placement = "outside",
                  strip.text = element_text(size = 8, color = "black"),
                  legend.background = element_blank(),
                  legend.title = element_blank(),
                  legend.key = element_blank(),
                  axis.ticks.length = unit(.25, "cm"),
                  axis.ticks = element_line(lineend = "round"), 
                  text = element_text(size = 12), 
                  axis.text.x = element_text(size = 12, color = "black", 
                                             margin = margin(t = 7)),
                  axis.text.x.top = element_text(size = 12, color = "black", 
                                                 margin = margin(b = 7)),
                  axis.text.y = element_text(size = 12, color = "black", 
                                             margin = margin(r = 7))
)}
TrtCol <- c("dodgerblue3", "brown3")
## Load compiled Data
Allvar.df <- read.csv("./SSIII_dryad.csv") %>% 
  ## ^Download dataset from dryad 
  mutate(Phase = ifelse(HarvestDate >= "2021-04-25", "Growth", "Storage" ),
         HarvestDate = HarvestDate %>% as.Date,
         FieldRhiHarvest = "2021-01-13" %>% as.Date,
         ExperimentStart = "2021-03-03" %>% as.Date,
         PlantedRhi = "2021-03-17" %>% as.Date,
         Days = (HarvestDate - PlantedRhi) %>% as.numeric,
         
         RhiSample_g = ifelse(Phase == "Storage", RhiInitialFW_g, RhiSample_g),
         Caputd = ifelse(RhiID %in% ripID$RhiID, "dead", "alive")
  )

### a) Outlier removal
Allvar.df[with(Allvar.df, !is.na(RhiSample_g) & Sampling == "Grain" & Caputd == "alive" & RhiSample_g < 6), 
          "RhiSample_g"] <- NA
Allvar.df[with(Allvar.df, !is.na(RhiSample_g) & Sampling == "Final" & Caputd == "alive" & RhiSample_g < 6), 
          "RhiSample_g"] <- NA

Allvar.df[with(Allvar.df, !is.na(RhiSample_g) & Caputd == "alive" & RhiSample_g > 27), 
          "RhiSample_g"] <- NA

Allvar.df[with(Allvar.df, !is.na(RhiSample_g) & Sampling == "Final" & Caputd == "alive" & Root_g < 6), 
          "Root_g"] <- NA

Allvar.df[with(Allvar.df, !is.na(AboveGround_g) & HarvestDate == "2021-08-25" & Caputd == "alive" & AboveGround_g < 1), "AboveGround_g"] <- NA

### b) remove Duplicated observations
for(i in c("AboveGround_g", "RhiSample_g", "Root_g", "Starch_Rhi")) {Allvar.df[, i] <- ifelse(!duplicated(Allvar.df[, i]), Allvar.df[, i], NA)}

###                   ###
### ---- STORAGE ----  
###                   ###

## ---- I. Rhizome Respiration ----

### ___ i) Model fitting ----

#### Linear
Rrhi_lm_cnt <- lm(RhiResp_nmol.g.min/1000 ~ -1 + Days*Treatment,
                  # SSlogis(Days, Asym, xmid, scal),
                  # params = list(Asym + scal + xmid ~ Treatment),
                  data = RhiResp )
Rrhi_lm_cnt %>% anova
Rrhi_lm_cnt %>% summary
Rrhi_lm_cnt %>% emtrends(~Treatment, var = "Days")  %>% multcomp::cld(Letters = letters, reversed = TRUE)
Rrhi_lm_cnt %>% emtrends(~Treatment, var = "Days")  %>% pwpm

Rrhi_lm_dcrt <- lm(RhiResp_nmol.g.min/1000 ~ Treatment*factor(Days), 
                      ##^ So it's on umolCO2
                # SSlogis(Days, Asym, xmid, scal),
                # params = list(Asym + scal + xmid ~ Treatment),
                data = RhiResp )
Rrhi_lm_dcrt %>% anova

### ___ ii) Plot ----

Rrhi.trt.emm <-  emmeans(Rrhi_lm_cnt, ~Treatment) %>% multcomp::cld(Letters = letters, reversed = TRUE)
emmeans(Rrhi_lm_cnt, ~Treatment) %>% pwpm()

## Estimate Q10: Q10 = (R2/R1)^10/(T2 - T1)
(0.0538/0.0101)^(10/(25 - 5))

Rrhi.emm.pvals <- emmeans(Rrhi_lm_dcrt, pairwise ~ ~Treatment|Days) %>% .$contrasts %>% 
  as.data.frame() %>% select(Days, p.value) %>% cbind(Age = "Room")

Rrhi_emm <- 
  Rrhi_lm_dcrt %>% emmeans(~Treatment|Days) %>% multcomp::cld(Letters = letters, reversed = TRUE) %>% 
  left_join(Rrhi.emm.pvals) %>% 
  group_by(Days) %>% 
  mutate(Signif = ifelse(p.value < 0.05, "*", 
                         ifelse(p.value > 0.05 & p.value < 0.1, ".", NA)),
         ypos = max(emmean) + SE)

Rrhi_lm_dcrt %>% emmeans(~Treatment|Days) %>% pwpm

wi.trt.lttr <- Rrhi_lm_dcrt %>% emmeans(~Days|Treatment) %>% 
  multcomp::cld(Letters = letters, reversed = TRUE) %>% as.data.frame() %>% 
  group_by(Days) %>% 
  mutate(SE = ifelse(Days == 0, 0, SE),
         ypos = ifelse(emmean == max(emmean), emmean + SE + SE*.5, emmean - SE - SE*.5),
         .group = if(all(.group == " a")) {""} else {gsub(" ", "", .group)}) %>% 
  as.data.frame() 

Rrhi_plt <-  
  Rrhi_emm %>% as.data.frame() %>% 
  mutate(SE = ifelse(Days == 14, 0, SE)) %>% 
  ggplot(aes(x = Days, y = emmean, col = Treatment)) +
  scale_color_manual(values = TrtCol) +
  
  geom_point(position = position_dodge(0.6), show.legend = T, size = 2) +
  
  geom_errorbar(aes(ymin = emmean  - SE, ymax = emmean  + SE),
                width = .5, show.legend = F, position = position_dodge(0.6)) +
  
  ## Add Trt diff
  geom_text(aes(y = 0, label =Signif), 
            col = "grey35", show.legend = F, cex = 8) +
  
  ## Add change diff
  # geom_text(data = wi.trt.lttr, aes(y = ypos, label = .group), show.legend = F) +
  
  # geom_smooth(method = "lm",
  #             formula = y ~ poly(x, 3), se = F, lwd = .5) +
  

  ## Add horizontal line for Rrhi treatment average
  geom_segment(data = Rrhi.trt.emm,
               aes(y = emmean, yend = emmean, col = Treatment), x = -15, xend = 0, 
               lty = 2, show.legend = F) +
  
  ## Add arrow to indicate carb sampling
  # annotate("segment", #inherit.aes = T,
  #          x = c(-14, -8, 0), xend = c(-14, -8, 0), y = 0.02, yend = 0.003,
  #          arrow = arrow(length = unit(.25, units = "cm"))) +
  
  scale_x_continuous(name = "Days before planting", limits = c(-14.5, 0.5), breaks = c(-14, -7, 0), labels = c(14, 7, 0)) +
  scale_y_continuous(name = latex2exp::TeX('$R_{rhi}$ ($\\mu mol$ $\\CO_{2}$ $\\g^{-1} s^{-1})'),
                     limits = c(0, .075), expand = c(0,0)) +
  myTheme + theme(legend.position = c(.9, .9))
Rrhi_plt


## ---- II. Rhizome Starch (%) ----
Storage.starch.df <- Allvar.df %>%  
  filter(!(is.na(Starch_Rhi))) %>% 
  select(RhiID, Treatment, Sampling, HarvestDate, Days, Caputd, Phase, Starch_Rhi, Sucrose_Rhi, Glucose_Rhi) %>% 
  # filter(Root_g < 17) %>% 
  
  filter(Caputd == "alive") %>%
  filter(Days <= 0)

### ___ i) Model fitting ----

#### Linear
Storage.starchRhi_lm_cnt <- lm(Starch_Rhi ~ Days*Treatment,
                               # SSlogis(Days, Asym, xmid, scal),
                               # params = list(Asym + scal + xmid ~ Treatment),
                               data = Storage.starch.df )
Storage.starchRhi_lm_cnt %>% anova
Storage.starchRhi_lm_cnt %>% summary
Storage.starchRhi_lm_cnt %>% emtrends(~Treatment, var = "Days")  %>% multcomp::cld(Letters = letters, reversed = TRUE)

Storage.starchRhi_lm_cnt %>% emtrends(~Treatment, var = "Days") %>% test

Storage.starchRhi_lm_dcrt <- lm(Starch_Rhi ~ Treatment*factor(Days), 
                                data = Storage.starch.df )
Storage.starchRhi_lm_dcrt %>% anova

Storage.starchRhi.emm.pvals <- emmeans(Storage.starchRhi_lm_dcrt, pairwise ~ Treatment|Days) %>% .$contrasts %>% 
  as.data.frame() %>% select(Days, p.value) %>% cbind(Treatment = "Room")

emmeans(Storage.starchRhi_lm_dcrt,  ~Days|Treatment) %>% pwpm
 ## Rhi at Room lost 10.6/15.7 -1 = 32% of starch
emmeans(Storage.starchRhi_lm_dcrt,  ~Treatment|Days) %>% pwpm
  ## Rhi at Room had  10.6/14.5 - 1 = 27% less starch than cold

Storage.starchRhi.emm <- 
  Storage.starchRhi_lm_dcrt %>% emmeans(~Treatment|Days) %>% multcomp::cld(Letters = letters, reversed = TRUE) %>% 
  left_join(Storage.starchRhi.emm.pvals) %>% 
  group_by(Days) %>% 
  mutate(Signif = ifelse(p.value < 0.05, "*", 
                         ifelse(p.value > 0.05 & p.value < 0.1, ".", NA)),
         ypos = max(emmean) + SE)

Storage.starchRhi_lm_dcrt %>% emmeans(~Treatment|Days) %>% pwpm
## Day 161: 14.8/12.9 - 1 (P = 0.5232)
## Day 161: 4.03/2.32 - 1 (P = 0.404)

wi.trt.lttr <- Storage.starchRhi_lm_dcrt %>% emmeans(~Days|Treatment) %>% 
  multcomp::cld(Letters = letters, reversed = TRUE) %>% as.data.frame() %>% 
  group_by(Days) %>% 
  mutate(ypos = ifelse(emmean == max(emmean), emmean + SE + SE*.5, emmean - SE - SE*.5),
         .group = if(all(.group == " a")) {""} else {gsub(" ", "", .group)}) %>% 
  as.data.frame()

### ___ ii) Plot ----
Storage.starchRhi.plt <- 
  Storage.starchRhi.emm %>% as.data.frame() %>% 
  ggplot(aes(x = Days, y = emmean, col = Treatment, lty = Treatment)) +
  scale_color_manual(values = TrtCol) +
  scale_linetype_manual(values = c(2,1)) +
  
  geom_point(show.legend = T, size = 2) +
  
  geom_errorbar(aes(ymin = emmean  - SE, ymax = emmean  + SE),
                width = .5, show.legend = F, lty = 1) +
  
  ## Add change diff
  geom_text(data = wi.trt.lttr, aes(y = ypos, label = .group), show.legend = F) +
  
  geom_text(aes(y = 18, label =Signif), 
            col = "grey35", show.legend = F, cex = 8) +
  
  geom_smooth(method = "lm", se = F, lwd = .5, show.legend = F) +
  
  scale_x_continuous(name = "Days before planting", limits = c(-14.5, 0.5), breaks = c(-14, -7, 0), labels = c(14, 7, 0)) +
  scale_y_continuous(name = "Rhizome Starch (%)", limits = c(-0, 20), expand = c(0,0)) +
  myTheme + theme(legend.position = c(.9, .1))
Storage.starchRhi.plt
FigS1 <- cowplot::plot_grid(Rrhi_plt, Storage.starchRhi.plt, ncol = 1,
                   align = "hv", labels = "AUTO")

## ___ iii) End of storage ----
endStorage.starch.plt <- Storage.starchRhi.emm %>% as.data.frame() %>% 
  filter(Days == 0) %>% 
  ggplot(aes(x = Treatment, y = emmean, fill = Treatment)) +
  scale_fill_manual(values = TrtCol) +
  
  geom_col(show.legend = F, col = "black") +
  
  geom_errorbar(aes(ymin = emmean  - SE, ymax = emmean  + SE),
                width = .5, show.legend = F, lty = 1) +
  
  ## Add change diff

  geom_text(aes(y = 18, label =Signif), 
            col = "grey35", show.legend = F, cex = 8) +
  
  annotate("text", x = "Room", y = 12, label = "b") +
  annotate("text", x = "Cold", y = 16, label = "a") +
  
  scale_y_continuous(name = "Rhizome Starch (%)", limits = c(-0, 18), expand = c(0,0)) +
  myTheme + theme(legend.position = c(.9, .1))

Storage.starchRhi.plt
# FigS1.1 <- 
cowplot::plot_grid(Rrhi_plt, endStorage.starch.plt, ncol = 2,
                   align = "hv", labels = "AUTO", rel_widths = c(7, 3))

## ---- II. Rhizome Sucrose (%) ----
Storage.starch.df <- Allvar.df %>%  
  filter(!(is.na(Starch_Rhi))) %>% 
  select(RhiID, Treatment, Sampling, HarvestDate, Days, Caputd, Phase, Starch_Rhi, Sucrose_Rhi, Glucose_Rhi) %>% 
  # filter(Root_g < 17) %>% 
  
  filter(Caputd == "alive") %>%
  filter(Days <= 0)

### ___ i) Model fitting ----

#### Linear
Storage.SucroseRhi_lm_cnt <- lm(Sucrose_Rhi ~ Days*Treatment,
                               # SSlogis(Days, Asym, xmid, scal),
                               # params = list(Asym + scal + xmid ~ Treatment),
                               data = Storage.starch.df )
Storage.SucroseRhi_lm_cnt %>% anova
Storage.SucroseRhi_lm_cnt %>% summary
Storage.SucroseRhi_lm_cnt %>% emtrends(~Treatment, var = "Days")  %>% multcomp::cld(Letters = letters, reversed = TRUE)

Storage.SucroseRhi_lm_cnt %>% emtrends(~Treatment, var = "Days") %>% test

Storage.SucroseRhi_lm_dcrt <- lm(Sucrose_Rhi ~ Treatment*factor(Days), 
                                data = Storage.starch.df )
Storage.SucroseRhi_lm_dcrt %>% anova

Storage.SucroseRhi.emm.pvals <- emmeans(Storage.SucroseRhi_lm_dcrt, pairwise ~ Treatment|Days) %>% .$contrasts %>% 
  as.data.frame() %>% select(Days, p.value) %>% cbind(Treatment = "Room")

emmeans(Storage.SucroseRhi_lm_dcrt,  ~Days|Treatment) %>% pwpm
## Rhi at Room lost 1.61/2.84 = 57% of Sucrose
emmeans(Storage.SucroseRhi_lm_dcrt,  ~Treatment|Days) %>% pwpm
## Rhi at Room had  1.61/3.06 = 27% less Sucrose than cold

Storage.SucroseRhi.emm <- 
  Storage.SucroseRhi_lm_dcrt %>% emmeans(~Treatment|Days) %>% multcomp::cld(Letters = letters, reversed = TRUE) %>% 
  left_join(Storage.SucroseRhi.emm.pvals) %>% 
  group_by(Days) %>% 
  mutate(Signif = ifelse(p.value < 0.05, "*", 
                         ifelse(p.value > 0.05 & p.value < 0.1, ".", NA)),
         ypos = max(emmean) + SE)

Storage.SucroseRhi_lm_dcrt %>% emmeans(~Treatment|Days) %>% pwpm
## Day 161: 14.8/12.9 - 1 (P = 0.5232)
## Day 161: 4.03/2.32 - 1 (P = 0.404)

wi.trt.lttr <- Storage.SucroseRhi_lm_dcrt %>% emmeans(~Days|Treatment) %>% 
  multcomp::cld(Letters = letters, reversed = TRUE) %>% as.data.frame() %>% 
  group_by(Days) %>% 
  mutate(ypos = ifelse(emmean == max(emmean), emmean + SE + SE*.5, emmean - SE - SE*.5),
         .group = if(all(.group == " a")) {""} else {gsub(" ", "", .group)}) %>% 
  as.data.frame()

### ___ ii) Plot ----
Storage.SucroseRhi.plt <- 
  Storage.SucroseRhi.emm %>% as.data.frame() %>% 
  ggplot(aes(x = Days, y = emmean, col = Treatment, lty = Treatment)) +
  scale_color_manual(values = TrtCol) +
  scale_linetype_manual(values = c(2,1)) +
  
  geom_point(show.legend = T, size = 2) +
  
  geom_errorbar(aes(ymin = emmean  - SE, ymax = emmean  + SE),
                width = .5, show.legend = F, lty = 1) +
  
  ## Add change diff
  geom_text(data = wi.trt.lttr, aes(y = ypos, label = .group), show.legend = F) +
  
  geom_text(aes(y = 3.5, label =Signif), 
            col = "grey35", show.legend = F, cex = 8) +
  
  geom_smooth(method = "lm", se = F, lwd = .5, show.legend = F) +
  
  scale_x_continuous(name = "Days before planting", limits = c(-14.5, 0.5), breaks = c(-14, -7, 0), labels = c(14, 7, 0)) +
  scale_y_continuous(name = "Rhizome Sucrose (%)", limits = c(-0, 4), expand = c(0,0)) +
  myTheme + theme(legend.position = c(.9, .1))
Storage.SucroseRhi.plt
# FigS1 <- cowplot::plot_grid(Rrhi_plt, Storage.SucroseRhi.plt, ncol = 1,
#                             align = "hv", labels = "AUTO")

## ---- III. Rhizome Glucose (%) ----
Storage.starch.df <- Allvar.df %>%  
  filter(!(is.na(Starch_Rhi))) %>% 
  select(RhiID, Treatment, Sampling, HarvestDate, Days, Caputd, Phase, Starch_Rhi, Sucrose_Rhi, Glucose_Rhi) %>% 
  # filter(Root_g < 17) %>% 
  
  filter(Caputd == "alive") %>%
  filter(Days <= 0)

### ___ i) Model fitting ----

#### Linear
Storage.GlucoseRhi_lm_cnt <- lm(Glucose_Rhi ~ Days*Treatment,
                                # SSlogis(Days, Asym, xmid, scal),
                                # params = list(Asym + scal + xmid ~ Treatment),
                                data = Storage.starch.df )
Storage.GlucoseRhi_lm_cnt %>% anova
Storage.GlucoseRhi_lm_cnt %>% summary
Storage.GlucoseRhi_lm_cnt %>% emtrends(~Treatment, var = "Days")  %>% 
  multcomp::cld(Letters = letters, reversed = TRUE)

Storage.GlucoseRhi_lm_cnt %>% emtrends(~Treatment, var = "Days") %>% test

Storage.GlucoseRhi_lm_dcrt <- lm(Glucose_Rhi ~ Treatment*factor(Days), 
                                 data = Storage.starch.df )
Storage.GlucoseRhi_lm_dcrt %>% anova

Storage.GlucoseRhi.emm.pvals <- emmeans(Storage.GlucoseRhi_lm_dcrt, pairwise ~ Treatment|Days) %>% .$contrasts %>% 
  as.data.frame() %>% select(Days, p.value) %>% cbind(Treatment = "Room")

emmeans(Storage.GlucoseRhi_lm_dcrt,  ~Days|Treatment) %>% pwpm
## Rhi at Room lost 1- 2.00 /2.45  = 18% of Glucose
emmeans(Storage.GlucoseRhi_lm_dcrt,  ~Treatment|Days) %>% pwpm
## Rhi at Room had  1-2.00 /2.97  = 32% less Glucose than cold

Storage.GlucoseRhi.emm <- 
  Storage.GlucoseRhi_lm_dcrt %>% emmeans(~Treatment|Days) %>% multcomp::cld(Letters = letters, reversed = TRUE) %>% 
  left_join(Storage.GlucoseRhi.emm.pvals) %>% 
  group_by(Days) %>% 
  mutate(Signif = ifelse(p.value < 0.05, "*", 
                         ifelse(p.value > 0.05 & p.value < 0.1, ".", NA)),
         ypos = max(emmean) + SE)

Storage.GlucoseRhi_lm_dcrt %>% emmeans(~Treatment|Days) %>% pwpm
## Day 161: 14.8/12.9 - 1 (P = 0.5232)
## Day 161: 4.03/2.32 - 1 (P = 0.404)

wi.trt.lttr <- Storage.GlucoseRhi_lm_dcrt %>% emmeans(~Days|Treatment) %>% 
  multcomp::cld(Letters = letters, reversed = TRUE) %>% as.data.frame() %>% 
  group_by(Days) %>% 
  mutate(ypos = ifelse(emmean == max(emmean), emmean + SE + SE*.5, emmean - SE - SE*.5),
         .group = if(all(.group == " a")) {""} else {gsub(" ", "", .group)}) %>% 
  as.data.frame()

### ___ ii) Plot ----
Storage.GlucoseRhi.plt <- 
  Storage.GlucoseRhi.emm %>% as.data.frame() %>% 
  ggplot(aes(x = Days, y = emmean, col = Treatment)) +
  scale_color_manual(values = TrtCol) +

  geom_point(show.legend = T, size = 2) +
  
  geom_errorbar(aes(ymin = emmean  - SE, ymax = emmean  + SE),
                width = .5, show.legend = F, lty = 1) +
  
  ## Add change diff
  geom_text(data = wi.trt.lttr, aes(y = ypos, label = .group), show.legend = F) +
  
  geom_text(aes(y = 3.5, label =Signif), 
            col = "grey35", show.legend = F, cex = 8) +
  
  geom_smooth(method = "lm", se = F, lwd = .5, show.legend = F, lty = 2) +
  
  scale_x_continuous(name = "Days before planting", limits = c(-14.5, 0.5), breaks = c(-14, -7, 0), labels = c(14, 7, 0)) +
  scale_y_continuous(name = "Rhizome Glucose (%)", limits = c(-0, 4), expand = c(0,0)) +
  myTheme + theme(legend.position = c(.9, .1))
Storage.GlucoseRhi.plt

FigS1_new <- cowplot::plot_grid(Storage.starchRhi.plt + 
                                  scale_x_continuous(name = NULL, labels = NULL) ,  
                                Storage.SucroseRhi.plt + scale_x_continuous(name = NULL, labels = NULL),  
                                Storage.GlucoseRhi.plt , 
                                ncol = 1, align = "hv", labels = "AUTO")
FigS1_new

## ---- IV. CO2 respiration & Starch consumption ----

### ___ i) total CO2 respired (nmol.g.min/1000) ----
co2Resp.df <- Rrhi.trt.emm %>% as.data.frame %>% 
  select(Treatment, CO2rep_umol.g.min = emmean) %>% 
  merge(data.frame(Days = c(8, 14))) %>% 
  mutate(CO2resp_umol.g = round(CO2rep_umol.g.min*60*24*Days),
         CO2resp_mg.g = (CO2resp_umol.g/1000000)*(12 + 16*2)*1000, 
         ## ^/10000000: umol to mol; CO2molar mass: (12 + 16*2); *1000: g to mg (to match starch)
         CO2resp_mgC.g = CO2resp_mg.g * 12/(12 + 16*2) ## proportion of C by mass on CO2
         )

### ___ ii) starch consumed (mg starch/g) ----
###         (As difference between days)
starchConsmd.df <- Storage.starchRhi.emm %>% as.data.frame %>% 
  select(Days, Treatment, starch_Percent = emmean) %>% 
  # merge(data.frame(Length_Day = c(7, 14))) %>% 
  group_by(Treatment) %>% 
  mutate(Days = abs(Days + 14),
         starch_mg.g = starch_Percent*10,
         starchLoss_mg.g = c(NA, -diff(starch_mg.g)),
         starchLoss_mgC.g = starchLoss_mg.g*12*6/(12*6 + 1*10 + 16*5),
  ) %>% 
  arrange(Treatment)

### ___ iii) Plot ----
co2vStarch.df <- left_join(co2Resp.df, starchConsmd.df) %>% 
  select(Treatment, Days, CO2resp_mgC.g, starchLoss_mgC.g) %>% 
  mutate(percentRespired = 100*CO2resp_mgC.g/starchLoss_mgC.g)


# Any units, just to see trends
co2vStarch.df %>% 
  ggplot(aes(x = CO2resp_umol.g, y = starchLoss_mg.g))+
  geom_point(aes(col = Treatment, pch = Days %>% as.factor)) + 
  geom_smooth(data = . %>% filter(CO2resp_umol.g > 150),
              method = "lm", fullrange = T, se = F)

# mg of C per gram of biomass
co2vStarch.plt <- co2vStarch.df %>% 
  ggplot(aes(starchLoss_mgC.g, CO2resp_mgC.g)) +
  scale_color_manual(values = TrtCol) +
  
  geom_point(aes(col = Treatment), size = 3) +
  geom_abline(slope = 1, intercept = 0, col = "grey") +
  
  geom_smooth(method = "lm", se = F) +
  
  scale_y_continuous(name = "C respired as CO2 (mg C/ g)", limits = c(0, 15), expand = c(0,0)) +
  scale_x_continuous(name = "C consumed as starch loss (mg C/ g)", limits = c(-.15, 15), expand = c(0,0)) +
  myTheme + theme(legend.position = c(.1, .95))

co2vStarch.df %>% 
  lm(CO2resp_mgC.g ~ starchLoss_mgC.g, data = .) %>% summary


### ___ ii) starch consumed (mg starch/g) ----
###     (As integrated rate of consumtion)
StarchLoss.df <- Storage.starchRhi_lm_cnt %>% emtrends(~Treatment, var = "Days") %>% as.data.frame %>% 
  select(Treatment, starchLossRate_percent.day = Days.trend) %>% 
  merge(data.frame(Days = c(8, 14))) %>% 
  mutate(starchLossRate_mg.g.day = -starchLossRate_percent.day*10,
         starchLoss_mg.g = round(starchLossRate_mg.g.day*Days, 2),
         starchLoss_mgC.g = starchLoss_mg.g*12*6/(12*6 + 1*10 + 16*5)
  ) 

### ___ iii) Plot ----
co2vStarch.df <- left_join(co2Resp.df, StarchLoss.df) %>% 
  select(Treatment, Days, CO2resp_mgC.g, starchLoss_mgC.g) %>% 
  mutate(percentRespired = 100*CO2resp_mgC.g/starchLoss_mgC.g)

# Any units, just to see trends
co2vStarch.df %>% 
  ggplot(aes(x = starchLoss_mgC.g, y = CO2resp_mgC.g))+
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  # geom_smooth(method = "lm", se = F) +
  scale_y_continuous(name = "C respired as CO2 (mg C/ g)", limits = c(0, 25), expand = c(0,0)) +
  scale_x_continuous(name = "C consumed from starch (mg C/ g)", limits = c(0, 25), expand = c(0,0)) +
  myTheme

### New Fig S1

FigS1_new <- cowplot::plot_grid(Storage.starchRhi.plt + 
                                  scale_x_continuous(name = NULL, labels = NULL) ,  
                                Storage.SucroseRhi.plt + scale_x_continuous(name = NULL, labels = NULL),  
                                Storage.GlucoseRhi.plt , 
                                ncol = 1, align = "hv", labels = "AUTO")
FigS1_new



###                                 ###
### ---- GROWTH & DEVELOPMENT ----  
###                                 ###

## ---- I. Belowground Biomass ----
belowBM.df <- Allvar.df %>% select(RhiID, Treatment, HarvestDate, Sampling, Days, Caputd, Phase,RhiInitialFW_g, RhiSample_g, Root_g) %>% 
  tidyr::as_tibble() %>% 
  filter(Caputd == "alive") %>%
  filter(Days >= 0) %>%
  filter(is.na(Root_g) | Root_g < 17)  %>% 
  filter(is.na(Root_g) | !(Sampling == "Final" & (Root_g < 5  | Root_g > 15))) %>% 
  mutate(belowBM = RhiSample_g + ifelse(Days <= 0, 0, Root_g),
         Root_g = ifelse(Sampling == "14-Day", 0, Root_g),
         Rel.RhiBM = RhiSample_g/RhiInitialFW_g ) 
belowBM.df %>% dim#filter(is.na(Root_g))  
belowBM.df %>% tidyr::as_tibble() %>% filter(is.na(Root_g) | Root_g < 17) %>% dim

### ___ i) nonlinear fit ----

### test with nlme::gnls()
belowBM_null <- gnls(belowBM ~ Asym/(1+exp((xmid-Days)/scal)),
                   data = belowBM.df %>% filter(!is.na(belowBM)),
                   start = list(Asym = rep(30, 1), xmid = rep(100, 1), scal = rep(20, 1)))

datatest <- belowBM.df %>% filter(!is.na(belowBM)) %>% 
  mutate(Treatment = as.factor(Treatment))
belowBM_full <- gnls(belowBM ~ Asym/(1+exp((xmid-Days)/scal)),
                   params = list(Asym  + xmid + scal ~ Treatment),
                   data = belowBM.df %>% filter(!is.na(belowBM)) %>% 
                     mutate(Treatment = as.factor(Treatment)),
                   start = list(Asym = rep(20, 2), xmid = rep(80, 2), scal = rep(20, 2)))
anova(belowBM_null, belowBM_full)
anova(belowBM_full)

belowBM_full %>% plot

belowBM_full %>% emmeans(param = "Asym", specs = ~Treatment) %>% multcomp::cld(Letters = letters, reversed = TRUE) 
belowBM_full %>% emmeans(param = "Asym", specs = ~Treatment) %>% pwpm
belowBM_full %>% emmeans(param = "xmid", specs = ~Treatment) %>% pwpm
belowBM_full %>% emmeans(param = "xmid", specs = ~Treatment) %>% multcomp::cld(Letters = letters, reversed = TRUE) 
belowBM_full %>% emmeans(param = "scal", specs = ~Treatment) %>% multcomp::cld(Letters = letters, reversed = TRUE) 

belowBM_full %>% predict(newdata = expand.grid(Days = 1:180,
                                               Treatment = c("Cold", "Room")))
belowBM_full %>% attributes()
belowBM_full$contrasts
# DaysTrt.df %>% with(table(Treatment, Days))

### ___ ii) Plot ----
belowBM.lm <- belowBM.df %>% 
  lm(belowBM ~ Treatment*factor(Days), data = .)

belowBM.emm <- belowBM.lm %>% emmeans(~Treatment|Days) %>% multcomp::cld(Letters = letters, reversed = TRUE)
belowBM.lm %>% emmeans(~Treatment|Days) %>% pwpm()

wi.trt.lttr <- belowBM.lm %>% emmeans(~Days|Treatment) %>% 
  multcomp::cld(Letters = letters, reversed = TRUE) %>% as.data.frame() %>% 
  # arrange(Days, desc(emmean)) %>% 
  group_by(Days) %>% 
  mutate(ypos = ifelse(emmean == max(emmean), emmean + SE + SE*.5, emmean - SE - SE*.5),
         .group = if(all(.group == " a")) {""} else {gsub(" ", "", .group)}) %>% 
  as.data.frame()

belowBM.plt <-
  belowBM.emm %>% as.data.frame() %>% 
  ggplot(aes(x = Days, y = emmean, col = Treatment)) +
  scale_color_manual(values = TrtCol) +
  
  geom_point(position = position_dodge(0.6), show.legend = T, size = 2) +
  
  geom_errorbar(aes(ymin = emmean  - SE, ymax = emmean  + SE),
                width = 5, show.legend = F, position = position_dodge(0.6)) +
  
  # ## Add Trt diff
  # geom_text(data = . %>% filter(Signif == "*"), aes(y = 1, label =Signif), 
  #           col = "grey35", show.legend = F, cex = 8) +
  
  ## Add change diff
  geom_text(data = wi.trt.lttr, aes(y = ypos, label = .group), show.legend = F) +
  
  geom_smooth(method = "nls",
              formula = y ~ SSlogis(x, Asym, xmid, scal), se = F, lwd = .5) +
  
  # geom_smooth(method = "lm", data = . %>% filter(Treatment == 'Room'), se = F, lwd = .5) +
  
  scale_x_continuous(name = "Days since planting", limits = c(-2, 165)) +
  scale_y_continuous(name = "Belowground Biomass (g)", limits = c(0, 30), expand = c(0,0)) +
  myTheme + theme(legend.position = c(.9, .1))
belowBM.plt

## ---- I. Rhizome Biomass ----
# rhiBM <- Allvar.df %>% select(RhiID, Treatment, HarvestDate, Days, Caputd, Phase, RhiSample_g, Root_g) %>% 
#   filter(!(is.na(RhiSample_g))) %>% 
#   filter(!(Days == 161 & (Root_g < 5  | Root_g > 15))) %>% 
#   # filter(Days >= 0) %>% 
#   filter(Caputd == "alive") 

### ___ i) nonlinear fit ----

### test with nlme::gnls()
rhiBM_null <- gnls(RhiSample_g ~ Asym/(1+exp((xmid-Days)/scal)),
                        # params = list(Asym + scal + xmid ~ Treatment),
                        data = belowBM.df %>% filter(!is.na(RhiSample_g)),
                        start = list(Asym = rep(10, 1), xmid = rep(100, 1), scal = rep(20, 1)))

rhiBM_full <- gnls(RhiSample_g ~ Asym/(1+exp((xmid-Days)/scal)),
     params = list(Asym  + xmid + scal ~ Treatment),
     data = belowBM.df %>% filter(!is.na(RhiSample_g))%>% 
       mutate(Treatment = as.factor(Treatment)),
     start = list(Asym = rep(10, 2), xmid = rep(80, 2), scal = rep(20, 2)))
anova(rhiBM_full)
anova(rhiBM_null, rhiBM_full)

rhiBM_full %>% plot

rhiBM_full %>% emmeans(param = "Asym", specs = ~Treatment) %>% multcomp::cld(Letters = letters, reversed = TRUE) 
rhiBM_full %>% emmeans(param = "Asym", specs = ~Treatment) %>% pwpm
rhiBM_full %>% emmeans(param = "xmid", specs = ~Treatment) %>% pwpm
rhiBM_full %>% emmeans(param = "xmid", specs = ~Treatment) %>% multcomp::cld(Letters = letters, reversed = TRUE) 
rhiBM_full %>% emmeans(param = "scal", specs = ~Treatment) %>% multcomp::cld(Letters = letters, reversed = TRUE) 

### ___ ii) Plot ----
rhiBM.lm <- belowBM.df %>% 
  lm(RhiSample_g ~ Treatment*factor(Days), data = .)

rhiBM.emm <- rhiBM.lm %>% emmeans(~Treatment|Days) %>% multcomp::cld(Letters = letters, reversed = TRUE)
rhiBM.lm %>% emmeans(~Treatment|Days) %>% pwpm()

wi.trt.lttr <- rhiBM.lm %>% emmeans(~Days|Treatment) %>% 
  multcomp::cld(Letters = letters, reversed = TRUE) %>% as.data.frame() %>% 
  # arrange(Days, desc(emmean)) %>% 
  group_by(Days) %>% 
  mutate(ypos = ifelse(emmean == max(emmean), emmean + SE + SE*.5, emmean - SE - SE*.5),
         .group = if(all(.group == " a")) {""} else {gsub(" ", "", .group)}) %>% 
  as.data.frame()

rhiBM.plt <- rhiBM.emm %>% as.data.frame() %>% 
  ggplot(aes(x = Days, y = emmean, col = Treatment)) +
  scale_color_manual(values = TrtCol) +
  
  geom_point(position = position_dodge(0.6), show.legend = T, size = 2) +
  
  geom_errorbar(aes(ymin = emmean  - SE, ymax = emmean  + SE),
                width = 5, show.legend = F, position = position_dodge(0.6)) +
  
  # ## Add Trt diff
  # geom_text(data = . %>% filter(Signif == "*"), aes(y = 1, label =Signif), 
  #           col = "grey35", show.legend = F, cex = 8) +
  
  ## Add change diff
  geom_text(data = wi.trt.lttr, aes(y = ypos, label = .group), show.legend = F) +
  
  geom_smooth(method = "nls",
              formula = y ~ SSlogis(x, Asym, xmid, scal), se = F, lwd = .5) +
  # geom_smooth(method = "lm", data = . %>% filter(Treatment == 'Room'), se = F, lwd = .5) +
  
  scale_x_continuous(name = "Days since planting", limits = c(-2, 165)) +
  scale_y_continuous(name = "Rhizome Biomass (g)", limits = c(0, 18), expand = c(0,0)) +
  myTheme + theme(legend.position = c(.9, .1))
rhiBM.plt
### ___ ii) Plot Relative ----
rel.rhiBM.lm <- belowBM.df %>% 
  lm(Rel.RhiBM ~ Treatment*factor(Days), data = .)

rel.rhiBM.emm <- rel.rhiBM.lm %>% emmeans(~Treatment|Days) %>% multcomp::cld(Letters = letters, reversed = TRUE)
rel.rhiBM.lm %>% emmeans(~Treatment|Days) %>% pwpm()

wi.trt.lttr <- rel.rhiBM.lm %>% emmeans(~Days|Treatment) %>% 
  multcomp::cld(Letters = letters, reversed = TRUE) %>% as.data.frame() %>% 
  # arrange(Days, desc(emmean)) %>% 
  group_by(Days) %>% 
  mutate(ypos = ifelse(emmean == max(emmean), emmean + SE + SE*.5, emmean - SE - SE*.5),
         .group = if(all(.group == " a")) {""} else {gsub(" ", "", .group)}) %>% 
  as.data.frame()

rel.rhiBM.plt <- rel.rhiBM.emm %>% as.data.frame() %>% 
  ggplot(aes(x = Days, y = emmean, col = Treatment)) +
  scale_color_manual(values = TrtCol) +
  
  geom_point(position = position_dodge(0.6), show.legend = T, size = 2) +
  
  geom_errorbar(aes(ymin = emmean  - SE, ymax = emmean  + SE),
                width = 5, show.legend = F, position = position_dodge(0.6)) +
  
  # ## Add Trt diff
  # geom_text(data = . %>% filter(Signif == "*"), aes(y = 1, label =Signif), 
  #           col = "grey35", show.legend = F, cex = 8) +
  
  ## Add change diff
  geom_text(data = wi.trt.lttr, aes(y = ypos, label = .group), show.legend = F) +
  
  geom_smooth(method = "nls",
              formula = y ~ SSlogis(x, Asym, xmid, scal), se = F, lwd = .5) +
  # geom_smooth(method = "lm", data = . %>% filter(Treatment == 'Room'), se = F, lwd = .5) +
  
  scale_x_continuous(name = "Days since planting", limits = c(-2, 165)) +
  scale_y_continuous(name = "rel Rhizome Biomass (g)", limits = c(0, 10), expand = c(0,0)) +
  myTheme + theme(legend.position = c(.9, .1))
rel.rhiBM.plt

## ---- II. Root Biomass ----
# rootBM <- Allvar.df %>%  
#   mutate(Root_g = ifelse(Sampling == "14-Day", 0, Root_g)) %>% 
#   filter(!(is.na(Root_g))) %>% 
#   select(RhiID, Treatment, Sampling, HarvestDate, Days, Caputd, Phase, Root_g) %>% 
#   filter(Root_g < 17) %>% 
#   
#   filter(Caputd == "alive") %>%
#   filter(Days >= 0) %>% 
#   filter(!(Sampling == "Final" & (Root_g < 5  | Root_g > 15))) 
  
### ___ i) nonlinear fit ----
rootBM_null <- gnls(Root_g ~ Asym/(1+exp((xmid-Days)/scal)),
                     # params = list(Asym + scal + xmid ~ Treatment),
                    data = belowBM.df %>% filter(!is.na(Root_g)),
                    start = list(Asym = rep(10, 1), xmid = rep(100, 1), scal = rep(20, 1)))
anova(rootBM_null)

rootBM_full <- gnls(Root_g ~ Asym/(1+exp((xmid-Days)/scal)),
                     params = list(Asym  + xmid + scal ~ Treatment),
                    data = belowBM.df %>% filter(!is.na(Root_g))%>% 
                      mutate(Treatment = as.factor(Treatment)),
                    start = list(Asym = rep(10, 2), xmid = rep(80, 2), scal = rep(20, 2)))
anova(rootBM_full)
anova(rootBM_null, rootBM_full)

rootBM_full %>% plot

rootBM_full %>% emmeans(param = "Asym", specs = ~Treatment) %>% pwpm #multcomp::cld(Letters = letters, reversed = TRUE) 
  ## 1 - 12.32/9.12 = 35% (P = 0.0431)
rootBM_full %>% emmeans(param = "xmid", specs = ~Treatment) %>% multcomp::cld(Letters = letters, reversed = TRUE) 
rootBM_full %>% emmeans(param = "xmid", specs = ~Treatment) %>% pwpm #multcomp::cld(Letters = letters, reversed = TRUE) 
## 24 days (P = 0.0141)

rootBM_full %>% emmeans(param = "scal", specs = ~Treatment) %>% multcomp::cld(Letters = letters, reversed = TRUE) 

### ___ ii) Plot ----
rootBM.lm <- belowBM.df %>% 
  lm(Root_g ~ Treatment*factor(Days), data = .)

rootBM.emm <- rootBM.lm %>% emmeans(~Treatment|Days) %>% multcomp::cld(Letters = letters, reversed = TRUE) %>% 
  group_by(Days) %>% 
  mutate(Signif = ifelse(all(.group == " a"), "", "*"),
         SE = ifelse(Days == 0, 0, SE),
         ypos = min(emmean) - SE)
rootBM.lm %>% emmeans(~Treatment|Days) %>% pwpm
  ## Day 64: 6.34/3.48 = ~ 82% (P = 0.0293)
  ## Day 128: 13.04/7.25 = ~ 79% (P = 0.0004)
##

wi.trt.lttr <- rootBM.lm %>% emmeans(~Days|Treatment) %>% 
  multcomp::cld(Letters = letters, reversed = TRUE) %>% as.data.frame() %>% 
  group_by(Days) %>% 
  mutate(SE = ifelse(Days == 0, 0, SE),
         ypos = ifelse(emmean == max(emmean), emmean + SE + SE*.5, emmean - SE - SE*.5),
         ypos = ifelse(Days == 0, 1, ypos),
         .group = if(all(.group == " a")) {""} else {gsub(" ", "", .group)}) %>% 
  as.data.frame()

rootBM.plt <- rootBM.emm %>% as.data.frame() %>% 
  mutate(SE = ifelse(Days == 14, 0, SE)) %>% 
  ggplot(aes(x = Days, y = emmean, col = Treatment)) +
  scale_color_manual(values = TrtCol) +
  
  geom_point(position = position_dodge(0.6), show.legend = T, size = 2) +
  
  geom_errorbar(aes(ymin = emmean  - SE, ymax = emmean  + SE),
                width = 5, show.legend = F, position = position_dodge(0.6)) +
  
  ## Add Trt diff
  geom_text(data = . %>% filter(Signif == "*"), aes(y = 1, label =Signif),
            col = "grey35", show.legend = F, cex = 8) +
  
  ## Add change diff
  geom_text(data = wi.trt.lttr, aes(y = ypos, label = .group), show.legend = F) +
  
  geom_smooth(method = "nls",
              formula = y ~ SSlogis(x, Asym, xmid, scal), se = F, lwd = .5) +
  
  scale_x_continuous(name = "Days since planting", limits = c(0, 165)) +
  scale_y_continuous(name = "Root Biomass (g)", limits = c(-0, 15.5), expand = c(0,0)) +
  myTheme + theme(legend.position = c(.9, .1))
rootBM.plt

## ---- III. Leaf Area ----
LeafArea.df <- Allvar.df %>%  
  mutate(Leaf.Area_cm2 = ifelse(Sampling == "14-Day", 0, Leaf.Area_cm2)) %>% 
  filter(!(is.na(Leaf.Area_cm2))) %>% 
  select(RhiID, Treatment, Sampling, HarvestDate, Days, Caputd, Phase, Leaf.Area_cm2)

### ___ ii) Plot ----
LeafArea.lm <- LeafArea.df %>% filter(Caputd == "alive") %>% 
  lm(Leaf.Area_cm2 ~ Treatment*factor(Days), data = .)

LeafArea.emm <- LeafArea.lm %>% emmeans(~Treatment|Days) %>% multcomp::cld(Letters = letters, reversed = TRUE) %>% 
  group_by(Days) %>% 
  mutate(Signif = ifelse(all(.group == " a"), "", "*"),
         ypos = min(emmean - SE))
LeafArea.lm %>% emmeans(~Treatment|Days) %>% pwpm()
    ## 239/ 149 - 1 = 60% (P = 0.013)

wi.trt.lttr <- LeafArea.lm %>% emmeans(~Days|Treatment) %>% 
  multcomp::cld(Letters = letters, reversed = TRUE) %>% as.data.frame() %>% 
  group_by(Days) %>% 
  mutate(SE = ifelse(Days == 0, 0, SE),
         ypos = ifelse(emmean == max(emmean), emmean + SE + SE*.5, emmean - SE - SE*.5),
         ypos = ifelse(Days == 0, 1, ypos),
         .group = if(all(.group == " a")) {""} else {gsub(" ", "", .group)}) %>% 
  as.data.frame()
LeafArea.lm %>% emmeans(~Days|Treatment) %>% pwpm()

LeafArea.plt <- LeafArea.emm %>% as.data.frame() %>% 
  mutate(SE = ifelse(Days == 0, 0, SE)) %>% 
  ggplot(aes(x = Days, y = emmean, col = Treatment)) +
  scale_color_manual(values = TrtCol) +
  
  geom_point(position = position_dodge(0.6), show.legend = T, size = 2) +
  
  geom_errorbar(aes(ymin = emmean  - SE, ymax = emmean  + SE),
                width = 5, show.legend = F, position = position_dodge(0.6)) +
  
  ## Add Trt diff
  geom_text(data = . %>% filter(Signif == "*"), aes(y = 1, label =Signif), 
            col = "grey35", show.legend = F, cex = 8) +
  
  ## Add change diff
  geom_text(data = wi.trt.lttr, aes(y = ypos, label = .group), show.legend = F) +
  geom_line() +
  # geom_smooth(#data = . %>% filter(Days != 128),
  #             method = "nls",
  #             formula = y ~ SSgompertz(x, Asym, b2, b3), se = F, lwd = .5) +
  # 
  # geom_smooth(method = "lm", formula = y ~ poly(x,2), se = F, lwd = .5) +
  
  scale_x_continuous(name = "Days since planting", limits = c(0, 165)) +
  scale_y_continuous(name = "Aboveground Biomass (g)", limits = c(-10, 300), expand = c(0,0)) +
  myTheme + theme(legend.position = c(.9, .1))

### ___ iii) LeafArea ~ LeafCt ----
LeafArea.lm %>% emmeans(~Treatment|Days) %>% pwpm()
Allvar.df %>% filter(Leaf.Area_cm2 < 400)

tottt.lm <- Allvar.df %>% filter(Caputd == "alive") %>% 
  lm(totLeaves ~ Treatment*factor(Days), data = .)

tott.emm <- tottt.lm %>% emmeans(~Treatment|Days) %>% multcomp::cld(Letters = letters, reversed = TRUE) %>% 
  group_by(Days) %>% 
  mutate(Signif = ifelse(all(.group == " a"), "", "*"),
         ypos = min(emmean - SE))


Allvar.df %>% 
  filter(Days <= 105) %>% 
  filter(Leaf.Area_cm2 < 400 & totLeaves > 5) %>% 
  ggplot(aes(x = totLeaves, y = Leaf.Area_cm2, col = Treatment)) + 
  geom_point() +
  geom_smooth(method = "lm", se = F) + 
  scale_x_continuous(limits = c(0, 42)) +
  scale_y_continuous(limits = c(0, 300))


left_join(LeafArea.emm %>% select(Treatment:SE) %>% dplyr::rename(LeafArea = emmean, seLeafArea = SE),
          tott.emm %>% select(Treatment:SE) %>% dplyr::rename(LeafCt = emmean, seLeafCt = SE)) %>% 
  filter(LeafArea < 239.0) %>% 
  filter(Days <= 105) %>% 
  ggplot(aes(x = LeafCt, y = LeafArea)) + 
  geom_point() +
  geom_smooth(method = "lm", se = F) +

  scale_x_continuous(limits = c(0, 42))

## Check root + rhizome = belowBM
left_join(
  rhiBM.emm %>% as.data.frame() %>% select(Treatment, Days, emmean) %>% dplyr::rename(Rhizome = emmean),
  rootBM.emm %>% as.data.frame() %>% select(Treatment, Days, emmean) %>% dplyr::rename(Root = emmean)
) %>% 
  left_join(belowBM.emm %>% as.data.frame() %>% select(Treatment, Days, emmean) %>% dplyr::rename(Below = emmean)) %>% 
  
  mutate(Total = Rhizome + Root)



## ---- IV. Aboveground Biomass ----
AboveBM <- Allvar.df %>%  
  mutate(AboveGround_g = ifelse(Sampling == "14-Day", 0, AboveGround_g)) %>% 
  filter(!(is.na(AboveGround_g))) %>% 
  select(RhiID, Treatment, Sampling, HarvestDate, Days, Caputd, Phase, AboveGround_g)


### ___ i) nonlinear fit ----

AboveBM_null <- gnls(AboveGround_g ~ Asym/(1+exp((xmid-Days)/scal)),
                     # params = list(Asym + scal + xmid ~ Treatment),
                     data = AboveBM %>% filter(Caputd == "alive") ,
                     start = list(Asym = rep(10, 1), xmid = rep(100, 1), scal = rep(20, 1)))

AboveBM_full <- gnls(AboveGround_g ~ Asym/(1+exp((xmid-Days)/scal)),
                     params = list(Asym  + xmid + scal ~ Treatment),
                     data = AboveBM %>% filter(Caputd == "alive") %>% 
                       mutate(Treatment = as.factor(Treatment)) ,
                     start = list(Asym = rep(10, 2), xmid = rep(80, 2), scal = rep(20, 2)))
anova(AboveBM_full)
anova(AboveBM_null, AboveBM_full)

AboveBM_full %>% plot

AboveBM_full %>% emmeans(param = "xmid", specs = ~Treatment) %>% pwpm()
AboveBM_full %>% emmeans(param = "scal", specs = ~Treatment) %>% pwpm()
AboveBM_full %>% emmeans(param = "Asym", specs = ~Treatment) %>% multcomp::cld(Letters = letters, reversed = TRUE)
AboveBM_full %>% emmeans(param = "xmid", specs = ~Treatment) %>% multcomp::cld(Letters = letters, reversed = TRUE)
AboveBM_full %>% emmeans(param = "scal", specs = ~Treatment) %>% multcomp::cld(Letters = letters, reversed = TRUE)

### ___ ii) Plot ----
AboveBM.lm <- AboveBM %>% filter(Caputd == "alive") %>% 
  lm(AboveGround_g ~ Treatment*factor(Days), data = .)

AboveBM.emm <- AboveBM.lm %>% emmeans(~Treatment|Days) %>% multcomp::cld(Letters = letters, reversed = TRUE) %>% 
  group_by(Days) %>% 
  mutate(Signif = ifelse(all(.group == " a"), "", "*"),
         ypos = min(emmean - SE))
AboveBM.lm %>% emmeans(~Treatment|Days) %>% pwpm()


wi.trt.lttr <- AboveBM.lm %>% emmeans(~Days|Treatment) %>% 
  multcomp::cld(Letters = letters, reversed = TRUE) %>% as.data.frame() %>% 
  group_by(Days) %>% 
  mutate(SE = ifelse(Days == 0, 0, SE),
         ypos = ifelse(emmean == max(emmean), emmean + SE + SE*.5, emmean - SE - SE*.5),
         ypos = ifelse(Days == 0, 1, ypos),
         .group = if(all(.group == " a")) {""} else {gsub(" ", "", .group)}) %>% 
  as.data.frame()
AboveBM.lm %>% emmeans(~Days|Treatment) %>% pwpm()
AboveBM.lm %>% emmeans(~Treatment|Days) %>% pwpm()

AboveBM.plt <- AboveBM.emm %>% as.data.frame() %>% 
  mutate(SE = ifelse(Days == 0, 0, SE)) %>% 
  ggplot(aes(x = Days, y = emmean, col = Treatment)) +
  scale_color_manual(values = TrtCol) +
  
  geom_point(position = position_dodge(0.6), show.legend = T, size = 2) +
  
  geom_errorbar(aes(ymin = emmean  - SE, ymax = emmean  + SE),
                width = 5, show.legend = F, position = position_dodge(0.6)) +
  
  ## Add Trt diff
  geom_text(data = . %>% filter(Signif == "*"), aes(y = 1, label =Signif), 
            col = "grey35", show.legend = F, cex = 8) +
  
  ## Add change diff
  geom_text(data = wi.trt.lttr, aes(y = ypos, label = .group), show.legend = F) +
  
  geom_smooth(method = "nls",
              formula = y ~ SSlogis(x, Asym, xmid, scal), se = F, lwd = .5) +
  
  # geom_smooth(method = "lm", data = . %>% filter(Treatment == 'Room'), se = F, lwd = .5) +
  
  scale_x_continuous(name = "Days since planting", limits = c(0, 165)) +
  scale_y_continuous(name = "Aboveground Biomass (g)", limits = c(-0, 12), expand = c(0,0)) +
  myTheme + theme(legend.position = c(.9, .1))

AboveBM.lm %>% emmeans(~Treatment|Days) %>% pwpm()

## This below is just to learn the shoot to rhizome ratio
## as a way to know how much starch is per area and how much
## carbon could be relocated.
left_join(
  rhiBM.emm %>% as.data.frame() %>% select(Treatment, Days, emmean) %>% dplyr::rename(rhiBM = emmean),
  belowBM.emm %>% as.data.frame() %>% select(Treatment, Days, emmean) %>% dplyr::rename(belowBM = emmean)
) %>% left_join(
  AboveBM.emm %>% select(Treatment, Days, emmean) %>% dplyr::rename(aboveBM = emmean)
) %>% 
  # filter(Days != 0) %>% 
  mutate(AboveRhiRatio = aboveBM/rhiBM,
         AboveBelowRatio = aboveBM/belowBM) %>% 
  
  ggplot(aes(x = Days, y = AboveRhiRatio, col = Treatment)) +
  geom_point()+
  geom_line()

## ___iii) estimate of carbon relocation ----
aveargeRhiBM <- 6 #Mg/Ha
starchRhi <- .15 # percent (first observation)
extraStarch <- .3 # percent
CcontentGlu = 12*6/(12*6 + 1*10 + 16*5)
# CcontentGlu = 12*6/(12*6 + 1*12 + 16*6)
CcontentABbm <- .4 ## C content in above ground biomass
aveargeRhiBM*starchRhi*extraStarch*1000*CcontentGlu*(1/CcontentABbm) # surplus grams of C from starch /ha


## ---- IV. Tillering ----

### This dataset comes from "SourceSinkIII_DevelopmentCompilation.R" in the DevelopmentFolder
### the function below works like a charm to just run the need lines of another script
source2 <- function(file, start, end, ...) {
  file.lines <- scan(file, what = character(), skip = start-1, nlines = end - start + 1, sep='\n')
  file.lines.collapsed <- paste(file.lines, collapse='\n')
  source(textConnection(file.lines.collapsed), ...)
}

#source2("//es.msu.edu/prl/labs/walker/General Project Storage/Mauri TN/Data&Analyses/swgSourceSink/Experiment III/Development/SourceSinkIII_DevelopmentCompilation.R",
#        start = 40, end = 95)
          ## ^ This information should be already in the Dryad file

captdRhiID <- Allvar.df %>% filter(Caputd != "alive") %>% .$RhiID %>% unique
Development$Caputd <- "alive"
Development[Development$RhiID %in% captdRhiID, "Caputd"] <- "dead"
Development <- Development %>% filter(!is.na(NumberOfTiller)) %>% 
  mutate(Days = as.numeric(Days))

### ___ i) nonlinear fit ----

Tiller_null <- gnls(NumberOfTiller ~ Asym*(1 - exp(-exp(lrc)*(Days - c0))),
                     # params = list(Asym + scal + xmid ~ Treatment),
                     data = Development, # %>% filter(Caputd == "alive") ,
                     start = list(Asym = rep(4, 1), lrc = rep(-3, 1), c0 = rep(-11, 1)))
summary(Tiller_null)
  
Tiller_full <- gnls(NumberOfTiller ~ Asym*(1 - exp(-exp(lrc)*(Days - c0))),
                     params = list(Asym  + lrc + c0 ~ Treatment),
                     data = Development, # %>% filter(Caputd == "alive") ,
                    start = list(Asym = rep(4, 2), lrc = rep(-3, 2), c0 = rep(-11, 2)))

anova(Tiller_full)
anova(Tiller_null, Tiller_full)

Tiller_full %>% plot

Tiller_full %>% emmeans(param = "Asym", specs = ~Treatment) %>% multcomp::cld(Letters = letters, reversed = TRUE)
Tiller_full %>% emmeans(param = "lrc", specs = ~Treatment) %>% multcomp::cld(Letters = letters, reversed = TRUE)
Tiller_full %>% emmeans(param = "c0", specs = ~Treatment) %>% multcomp::cld(Letters = letters, reversed = TRUE)

### ___ ii) Plot ----
Tillers.lm <- Development %>% filter(Caputd == "alive") %>% 
  lm(NumberOfTiller ~ Treatment*factor(Days), data = .)

Tillers.emm <- Tillers.lm %>% emmeans(~Treatment|Days) %>% multcomp::cld(Letters = letters, reversed = TRUE) %>% 
  group_by(Days) %>% 
  mutate(Signif = ifelse(all(.group == " a"), "", "*"),
         ypos = min(emmean - SE))

Tillers.lm %>% emmeans(~Treatment|Days) %>% pwpm()
  ## 5.5/3.75-1 = .46 (P = 0.0092)
wi.trt.lttr <- Tillers.lm %>% emmeans(~Days|Treatment) %>% 
  multcomp::cld(Letters = letters, reversed = TRUE) %>% as.data.frame() %>% 
  group_by(Days) %>% 
  mutate(SE = ifelse(Days == 0, 0, SE),
         ypos = ifelse(emmean == max(emmean), emmean + SE + SE*.5, emmean - SE - SE*.5),
         ypos = ifelse(Days == 0, 1, ypos),
         .group = if(all(.group == " a")) {""} else {gsub(" ", "", .group)}) %>% 
  as.data.frame()
Tillers.lm %>% emmeans(~Days|Treatment) %>% pwpm()
    ## Room - Day 8 v. Day 121: 5.50/1.83 = 3.0 (P < 0.0001)
    ## Cold - Day 8 v. Day 121: 3.75/1.92 = 

Tillers.plt <- 
  Tillers.emm %>% as.data.frame() %>% 
  mutate(SE = ifelse(Days == 0, 0, SE)) %>% 
  ggplot(aes(x = Days, y = emmean, col = Treatment)) +
  scale_color_manual(values = TrtCol) +
  
  geom_point(position = position_dodge(0.6), show.legend = T, size = 2) +
  
  geom_errorbar(aes(ymin = emmean  - SE, ymax = emmean  + SE),
                width = 5, show.legend = F, position = position_dodge(0.6)) +
  
  ## Add Trt diff
  geom_text(data = . %>% filter(Signif == "*"), aes(y = .2, label =Signif), 
            col = "grey35", show.legend = F, cex = 8) +
  
  ## Add change diff
  # geom_text(data = wi.trt.lttr, aes(y = ypos, label = .group), show.legend = F) +
  
  geom_smooth(method = "nls",
              formula = y ~ SSasympOff(x, Asym, lrc, c0), se = F, lwd = .5,
              fullrange = T) +
  
  # geom_smooth(method = "lm", data = . %>% filter(Treatment == 'Room'), se = F, lwd = .5) +
  
  
  scale_x_continuous(name = "Days since planting", limits = c(6, 165)) +
  scale_y_continuous(name = "Tiller Development (count/pot)", limits = c(-0, 6.5), expand = c(0,0)) +
  myTheme + theme(legend.position = c(.9, .1))

## ---- V. Leaf Development ----
Leaf.Flower.df <- read.csv("//es.msu.edu/prl/labs/walker/General Project Storage/Mauri TN/Data&Analyses/swgSourceSink/Experiment III/Development/SSIII_DevelopmentCompilation.csv") %>% 
  mutate(Days = (as.Date(SamplingDate) - as.Date(PlantedRhi)) %>%  as.numeric ) %>% 
  filter(!is.na(totLeaves))
         
captdRhiID <- Allvar.df %>% filter(Caputd != "alive") %>% .$RhiID %>% unique
Leaf.Flower.df$Caputd <- "alive"
Leaf.Flower.df[Leaf.Flower.df$RhiID %in% captdRhiID, "Caputd"] <- "dead"

orderedRhiID <- Leaf.Flower.df %>% .$RhiID %>% table %>% as.data.frame() %>% arrange(Freq) %>% .$. %>% as.character() %>%  as.numeric

Leaf.Flower.df %>% 
  with(., table(RhiID, SamplingDate, Treatment)) %>% as.data.frame() %>%
  filter(Freq != 0) %>% 
  
  mutate(RhiID = factor(RhiID, levels = orderedRhiID)) %>%
  arrange(RhiID) %>%
  # str
  ggplot(aes(x = SamplingDate %>% as.Date(), y = RhiID, col = Treatment)) +
  scale_color_manual(values = TrtCol) +
  facet_grid(Treatment ~., scales = "free_y" ) +
  geom_point(show.legend = F) +
  geom_line(show.legend = F) +
  myTheme

### ___ i) nonlinear fit ----
Leaf_null <- gnls(totLeaves ~ SSlogis(Days, Asym, xmid, scal),
                    # params = list(Asym + scal + xmid ~ Treatment),
                    data = Leaf.Flower.df , # %>% filter(Caputd == "alive")
                  )
summary(Leaf_null)

Leaf_full <- gnls(totLeaves ~ Asym/(1+exp((xmid-Days)/scal)),
                    params = list(Asym + scal + xmid ~ Treatment),
                    data = Leaf.Flower.df %>% 
                    mutate(Treatment = as.factor(Treatment)), # %>% filter(Caputd == "alive") ,
                    start = list(Asym = rep(22, 2), xmid = rep(50, 2), c0 = rep(23, 2)))

anova(Tiller_full)
anova(Leaf_null, Leaf_full)

Leaf_full %>% plot

Leaf_full %>% emmeans(param = "Asym", specs = ~Treatment) %>% multcomp::cld(Letters = letters, reversed = TRUE)
Leaf_full %>% emmeans(param = "xmid", specs = ~Treatment) %>% multcomp::cld(Letters = letters, reversed = TRUE)
Leaf_full %>% emmeans(param = "scal", specs = ~Treatment) %>% multcomp::cld(Letters = letters, reversed = TRUE)

### ___ ii) Plot ----
Leaf.lm <- Leaf.Flower.df %>% #filter(Caputd == "alive") %>% 
  lm(totLeaves ~ Treatment*factor(Days), data = .)

Leaf.emm <- Leaf.lm %>% emmeans(~Treatment|Days) %>% multcomp::cld(Letters = letters, reversed = TRUE) %>% 
  group_by(Days) %>% 
  mutate(Signif = ifelse(all(.group == " a"), "", "*"),
         ypos = min(emmean - SE))
Leaf.lm %>% emmeans(~Treatment|Days) %>% pwpm
  ## Day 143: 28.2/18.6 = 51% (P = 0.0119)
wi.trt.lttr <- Leaf.lm %>% emmeans(~Days|Treatment) %>% 
  multcomp::cld(Letters = letters, reversed = TRUE) %>% as.data.frame() %>% 
  group_by(Days) %>% 
  mutate(SE = ifelse(Days == 0, 0, SE),
         ypos = ifelse(emmean == max(emmean), emmean + SE + SE*.5, emmean - SE - SE*.5),
         ypos = ifelse(Days == 0, 1, ypos),
         .group = if(all(.group == " a")) {""} else {gsub(" ", "", .group)}) %>% 
  as.data.frame()

Leaf.plt <- Leaf.emm %>% as.data.frame() %>% 
  mutate(SE = ifelse(Days == 0, 0, SE)) %>% 
  ggplot(aes(x = Days, y = emmean, col = Treatment)) +
  scale_color_manual(values = TrtCol) +
  
  geom_point(position = position_dodge(0.6), show.legend = T, size = 2) +
  
  geom_errorbar(aes(ymin = emmean  - SE, ymax = emmean  + SE),
                width = 5, show.legend = F, position = position_dodge(0.6)) +
  
  ## Add Trt diff
  geom_text(data = . %>% filter(Signif == "*"), aes(y = .2, label =Signif), 
            col = "grey35", show.legend = F, cex = 8) +
  
  ## Add change diff
  geom_text(data = wi.trt.lttr, aes(y = ypos, label = .group), show.legend = F) +
  
  geom_smooth(method = "nls",
              formula = y ~ SSlogis(x, Asym, xmid, scal), se = F, lwd = .5,
              fullrange = T) +
  
  # geom_smooth(method = "lm", data = . %>% filter(Treatment == 'Room'), se = F, lwd = .5) +
  
  
  scale_x_continuous(name = "Days since planting", limits = c(6, 165)) +
  scale_y_continuous(name = "Leaf Development (count/pot)", limits = c(-0, 33), expand = c(0,0)) +
  myTheme + theme(legend.position = c(.9, .1))

## ---- VI. Leaf/Tiller ----
#Leaf.Flower.df <- read.csv("//es.msu.edu/prl/labs/walker/General Project Storage/Mauri TN/Data&Analyses/swgSourceSink/Experiment III/Development/SSIII_DevelopmentCompilation.csv") %>% 
#  mutate(Days = (as.Date(SamplingDate) - as.Date(PlantedRhi)) %>%  as.numeric ) %>% 
#  filter(!is.na(totLeaves))
      ## ^ This information is in the dryad dataset

captdRhiID <- Allvar.df %>% filter(Caputd != "alive") %>% .$RhiID %>% unique
Leaf.Flower.df$Caputd <- "alive"
Leaf.Flower.df[Leaf.Flower.df$RhiID %in% captdRhiID, "Caputd"] <- "dead"

Leaf.Flower.df %>% 
  filter(!(Days > 50 & avgLeaves < 1.2 | 
             Days > 75 & avgLeaves < 2.5 |
             Days > 120 & avgLeaves < 4.2)) %>% 
  
  ggplot(aes(x = Days, y = avgLeaves , col = Treatment)) +
  # facet_grid(~ Caputd) +
  scale_color_manual(values = TrtCol) +
  # geom_errorbar(aes(ymin = avgTillerCount - seTillerCount, ymax = avgTillerCount + seTillerCount), lwd = 1,
  #               width = 5, show.legend = F, position = "identity") +
  geom_smooth(method = "lm", lwd = 1, fullrange = F, se = F) +
  geom_jitter(size = 2, width = 2, show.legend = T, alpha = .2) +
  myTheme + theme(legend.position = c(.1, .9))

### ___ i) linear fit ----
leaf.till.lm <- lm(avgLeaves ~ Days*Treatment, 
                   data = Leaf.Flower.df %>% filter(!(Days > 50 & avgLeaves < 1.2 | 
                                                        Days > 75 & avgLeaves < 2.5 |
                                                        Days > 120 & avgLeaves < 4.2)) )
anova(leaf.till.lm)
emtrends(leaf.till.lm, specs = "Treatment", var = "Days") %>% 
  multcomp::cld(Letters = letters, reversed = TRUE)

### ___ ii) Plot ----
LeafTill.lm <- Leaf.Flower.df %>% #filter(Caputd == "alive") %>% 
  lm(avgLeaves ~ Treatment*factor(Days), data = .)

LeafTill.emm <- LeafTill.lm %>% emmeans(~Treatment|Days) %>% multcomp::cld(Letters = letters, reversed = TRUE) %>% 
  group_by(Days) %>% 
  mutate(Signif = ifelse(all(.group == " a"), "", "*"),
         ypos = min(emmean - SE))

wi.trt.lttr <- LeafTill.lm %>% emmeans(~Days|Treatment) %>% 
  multcomp::cld(Letters = letters, reversed = TRUE) %>% as.data.frame() %>% 
  group_by(Days) %>% 
  mutate(SE = ifelse(Days == 0, 0, SE),
         ypos = ifelse(emmean == max(emmean), emmean + SE + SE*.5, emmean - SE - SE*.5),
         ypos = ifelse(Days == 0, 1, ypos),
         .group = if(all(.group == " a")) {""} else {gsub(" ", "", .group)}) %>% 
  as.data.frame()

LeafTill.plt <- LeafTill.emm %>% as.data.frame() %>% 
  mutate(SE = ifelse(Days == 0, 0, SE)) %>% 
  ggplot(aes(x = Days, y = emmean, col = Treatment)) +
  scale_color_manual(values = TrtCol) +
  
  geom_point(position = position_dodge(0.6), show.legend = T, size = 2) +
  
  geom_errorbar(aes(ymin = emmean  - SE, ymax = emmean  + SE),
                width = 5, show.legend = F, position = position_dodge(0.6)) +
  
  ## Add Trt diff
  geom_text(data = . %>% filter(Signif == "*"), aes(y = .2, label =Signif), 
            col = "grey35", show.legend = F, cex = 8) +
  
  ## Add change diff
  # geom_text(data = wi.trt.lttr, aes(y = ypos, label = .group), show.legend = F) +
  
  geom_smooth(method = "nls",
              formula = y ~ SSlogis(x, Asym, xmid, scal), se = F, lwd = .5,
              fullrange = T) +
  
  # geom_smooth(method = "lm", data = . %>% filter(Treatment == 'Room'), se = F, lwd = .5) +
  
  
  scale_x_continuous(name = "Days since planting", limits = c(6, 165)) +
  scale_y_continuous(name = "leaves per tiller (leaf/tiller)", limits = c(-0, 7.5), expand = c(0,0)) +
  myTheme + theme(legend.position = c(.9, .1))
LeafTill.plt


### VII. Flowering

# Remove all samplingdates that are zero

# Leaf.Flower.df %>% filter(avgFlower > 0)
Flor.lm <- Leaf.Flower.df %>% filter(avgFlower > 0) %>% 
  lm(avgFlower  ~ Treatment, data = .)
#       Cold   Room
# Cold [6.11] 0.1395
# Room   2.67 [3.44]
  
# Flor.emm <- Flor.lm %>% emmeans(~Treatment) %>% multcomp::cld(Letters = letters, reversed = TRUE) %>% 
Flor.lm %>% emmeans(~Treatment) 
Flor.lm %>% emmeans(~Treatment) %>% pwpm()

wi.trt.lttr <- LeafTill.lm %>% emmeans(~Days|Treatment) %>% 
  multcomp::cld(Letters = letters, reversed = TRUE) %>% as.data.frame() %>% 
  group_by(Days) %>% 
  mutate(SE = ifelse(Days == 0, 0, SE),
         ypos = ifelse(emmean == max(emmean), emmean + SE + SE*.5, emmean - SE - SE*.5),
         ypos = ifelse(Days == 0, 1, ypos),
         .group = if(all(.group == " a")) {""} else {gsub(" ", "", .group)}) %>% 
  as.data.frame()

###                           ###
### ---- PHOTOSYNTHESIS ----  
###                           ###

## ---- I. Net CO2 Assimilation ----
#Photos.df <- read.csv("//es.msu.edu/prl/labs/walker/General Project Storage/Mauri TN/Data&Analyses/swgSourceSink/Experiment III/SourceActivity/SSIII_SurveyCompilation.csv") %>% 
#  mutate(Phase = "Growth",
#         SamplingDate = as.Date(SamplingDate),
#         FieldRhiHarvest = "2021-01-13" %>% as.Date,
#         ExperimentStart = "2021-03-03" %>% as.Date,
#         PlantedRhi = "2021-03-17" %>% as.Date,
#         Days = (SamplingDate - ExperimentStart) %>% as.numeric
#  ) %>%  
#  dplyr::rename(RhiID = PotNumber,
#                Treatment = Storage) %>% 
#  select(RhiID, Treatment, SamplingDate, Days, Phase, A, PhiPS2, gsw) 
         ## ^ This data is in the dryad dataset. 

### ___ i) Model fitting ----

#### Linear
Anet_lm_cnt <- lm(A ~ Days*Treatment,
                  # SSlogis(Days, Asym, xmid, scal),
                  # params = list(Asym + scal + xmid ~ Treatment),
                  data = Photos.df )
Anet_lm_cnt %>% anova
Anet_lm_cnt %>% summary
Anet_lm_cnt %>% emtrends(~Treatment, var = "Days")  %>% multcomp::cld(Letters = letters, reversed = TRUE)
Anet_lm_cnt %>% emtrends(~Treatment, var = "Days")  %>% pwpm

A_lm_dcrt <- lm(A ~ Treatment*factor(Days), 
                # SSlogis(Days, Asym, xmid, scal),
                # params = list(Asym + scal + xmid ~ Treatment),
                data = Photos.df )
A_lm_dcrt %>% anova

### x) Try segmented approach
Anet.sgmtd.df <- Photos.df %>% mutate(trtRoom = ifelse(Treatment == "Room", Days, 0),
                                      trtCold = ifelse(Treatment == "Cold", Days, 0))

Anet.sgmtd.lm <- lm(A ~ Treatment + trtRoom + trtCold , ## the - is to make the left slope zero
                    data = Anet.sgmtd.df)

#### Add Break-Point analysis
Anet.sgmtd.lm.Seg <- segmented(Anet.sgmtd.lm, seg.Z= ~ trtRoom + trtCold, 
                               psi = list(trtRoom = c(60), trtCold = c(100)))

Anet.sgmtd.lm.Seg %>% slope
summary(Anet.sgmtd.lm.Seg)
confint(Anet.sgmtd.lm.Seg)


data.frame(Treatment = Anet.sgmtd.df$Treatment,
           Days = Anet.sgmtd.df$Days, 
           A = Anet.sgmtd.df$A,
           Pred = fitted(Anet.sgmtd.lm.Seg)) %>% 
  ggplot(aes(x = Days, y = Pred, col = Treatment)) +
  geom_point(aes(y = A))+
  geom_line()

### ___ ii) Plot ----
Anet.emm.pvals <- emmeans(A_lm_dcrt, pairwise ~ ~Treatment|Days) %>% .$contrasts %>% 
  as.data.frame() %>% select(Days, p.value) %>% cbind(Age = "Room")

A_emm <- 
  A_lm_dcrt %>% emmeans(~Treatment|Days) %>% multcomp::cld(Letters = letters, reversed = TRUE) %>% 
  left_join(Anet.emm.pvals) %>% 
  group_by(Days) %>% 
  mutate(Signif = ifelse(p.value < 0.05, "*", 
                         ifelse(p.value > 0.05 & p.value < 0.1, ".", NA)),
         ypos = max(emmean) + SE)

A_lm_dcrt %>% emmeans(~Treatment|Days) %>% pwpm

wi.trt.lttr <- A_lm_dcrt %>% emmeans(~Days|Treatment) %>% 
  multcomp::cld(Letters = letters, reversed = TRUE) %>% as.data.frame() %>% 
  group_by(Days) %>% 
  mutate(SE = ifelse(Days == 0, 0, SE),
         ypos = ifelse(emmean == max(emmean), emmean + SE + SE*.5, emmean - SE - SE*.5),
         .group = if(all(.group == " a")) {""} else {gsub(" ", "", .group)}) %>% 
  as.data.frame() 

A_plt <-  A_emm %>% as.data.frame() %>% 
  mutate(SE = ifelse(Days == 14, 0, SE)) %>% 
  ggplot(aes(x = Days, y = emmean, col = Treatment)) +
  scale_color_manual(values = TrtCol) +
  
  geom_point(position = position_dodge(0.6), show.legend = T, size = 2) +
  
  geom_errorbar(aes(ymin = emmean  - SE, ymax = emmean  + SE),
                width = 5, show.legend = F, position = position_dodge(0.6)) +
  
  ## Add Trt diff
  geom_text(aes(y = 25, label =Signif), 
            col = "grey35", show.legend = F, cex = 8) +
  
  ## Add change diff
  # geom_text(data = wi.trt.lttr, aes(y = ypos, label = .group), show.legend = F) +
  
  # geom_smooth(method = "lm",
  #             formula = y ~ poly(x, 3), se = F, lwd = .5) +
  
  # geom_vline(xintercept = c(80)) +
  
  geom_smooth(method = "lm",  se = F, lwd = .5) +
  
  scale_x_continuous(name = "Days since planting", limits = c(0, 165)) +
  scale_y_continuous(name = "Anet", limits = c(-3, 27)) +
  myTheme + theme(legend.position = c(.9, .9))
A_plt

A_lm1 <- lm(A ~ poly(Days,1), data = Photos.df )
A_lm3 <- lm(A ~ Days*Treatment + Days2*Treatment + Days3*Treatment, 
            data = Photos.df %>% mutate(Days2 = Days^2, Days3 = Days^3) )

anova(A_lm1, A_lm3)

### ___ iii) Accumulated Anet ----
A_emm %>% 
  group_by(Treatment) %>% 
  arrange(Days) %>% 
  mutate(lengthA = cumsum(emmean*0.035185),  ## correlation between Anet in (umol m-2 s-1) and A'net (mol m-2 day-1)
         widthDay = c(diff(Days), NA) %>% cumsum,
         accAnet = lengthA*widthDay/2
  ) %>% 
  
  ggplot(aes(x = Days, y = accAnet, col = Treatment)) +
  scale_color_manual(values = TrtCol) +
  
  geom_point(show.legend = T, size = 2) +
  
  # geom_errorbar(aes(ymin = emmean  - SE, ymax = emmean  + SE),
  #               width = 5, show.legend = F, position = position_dodge(0.6)) +
  
  ## Add Trt diff
  # geom_text(data = . %>% filter(Signif == "*"), aes(y = 25, label =Signif), 
  #           col = "grey35", show.legend = F, cex = 8) +
  
  ## Add change diff
  # geom_text(data = wi.trt.lttr, aes(y = ypos, label = .group), show.legend = F) +
  
geom_smooth(method = "lm",
            formula = y ~ poly(x, 2), se = F, lwd = .5) +
  
  # geom_smooth(method = "lm", data = . %>% filter(Treatment == 'Room'), se = F, lwd = .5) 
  
  scale_x_continuous(name = "Days since planting", limits = c(6, 165)) +
  scale_y_continuous(name = "Accumulated Net CO2 Assimilation (mol CO2 m-2)", limits = c(-15, 400)) +
  myTheme + theme(legend.position = c(.9, .1))

## ---- II. gsw ----

### ___ i) Model fitting ----

#### Linear
gsw_lm_cnt <- lm(gsw ~ Days*Treatment,
                 # SSlogis(Days, Asym, xmid, scal),
                 # params = list(Asym + scal + xmid ~ Treatment),
                 data = Photos.df %>% filter(gsw > 0 & gsw < .3))
gsw_lm_cnt %>% anova
gsw_lm_cnt %>% summary
gsw_lm_cnt %>% emtrends(~Treatment, var = "Days")  %>% multcomp::cld(Letters = letters, reversed = TRUE)
gsw_lm_cnt %>% emtrends(~Treatment, var = "Days")  %>% pwpm

gsw.lm_dcrt <- lm(gsw ~ Treatment*factor(Days), 
                  # SSlogis(Days, Asym, xmid, scal),
                  # params = list(Asym + scal + xmid ~ Treatment),
                  data = Photos.df %>% filter(gsw > 0 & gsw < .3))
gsw.lm_dcrt %>% anova

### ___ ii) Plot ----
gsw.emm <- 
  gsw.lm_dcrt %>% emmeans(~Treatment|Days) %>% multcomp::cld(Letters = letters, reversed = TRUE) %>% 
  group_by(Days) %>% 
  mutate(Signif = ifelse(all(.group == " a"), "", "*"),
         ypos = min(emmean - SE))

wi.trt.lttr <- gsw.lm_dcrt %>% emmeans(~Days|Treatment) %>% 
  multcomp::cld(Letters = letters, reversed = TRUE) %>% as.data.frame() %>% 
  group_by(Days) %>% 
  mutate(SE = ifelse(Days == 0, 0, SE),
         ypos = ifelse(emmean == max(emmean), emmean + SE + SE*.5, emmean - SE - SE*.5),
         .group = if(all(.group == " a")) {""} else {gsub(" ", "", .group)}) %>% 
  as.data.frame() %>% 
  left_join(
    data.frame(Days = Photos.df %>% arrange(Days) %>% .$Days %>% unique,
               Treatment = rep(c("Cold", "Room"), each = 20),
               ypos2 = c(seq(0, -8, length.out =20),
                         seq(-2, -10, length.out =20))
    )
  ) 

gsw.plt <- 
  gsw.emm %>% as.data.frame() %>% 
  mutate(SE = ifelse(Days == 14, 0, SE)) %>% 
  ggplot(aes(x = Days, y = emmean, col = Treatment)) +
  scale_color_manual(values = TrtCol) +
  
  geom_point(position = position_dodge(0.6), show.legend = T, size = 2) +
  
  geom_errorbar(aes(ymin = emmean  - SE, ymax = emmean  + SE),
                width = 5, show.legend = F, position = position_dodge(0.6)) +
  
  ## Add Trt diff
  geom_text(data = . %>% filter(Signif == "*"), aes(y = 0, label =Signif), 
            col = "grey35", show.legend = F, cex = 8) +
  
  ## Add change diff
  # geom_text(data = wi.trt.lttr, aes(y = ypos, label = .group), show.legend = F) +
  
  # geom_smooth(method = "lm",
  #             formula = y ~ poly(x, 3), se = F, lwd = .5) +
  
  geom_smooth(method = "lm",  se = F, lwd = .5) +
  
  scale_x_continuous(name = "Days since planting", limits = c(0, 165)) +
  scale_y_continuous(name = "gsw", limits = c(0, .15)) +
  myTheme + theme(legend.position = c(.9, .9))
gsw.plt

## ---- III. PWUE ----

### ___ i) Model fitting ----

#### Linear
PWUE.lm_cnt <- lm(A/gsw ~ Days*Treatment,
                  # SSlogis(Days, Asym, xmid, scal),
                  # params = list(Asym + scal + xmid ~ Treatment),
                  data = Photos.df %>% filter(gsw > 0 & gsw < .3))
PWUE.lm_cnt %>% anova
PWUE.lm_cnt %>% summary
PWUE.lm_cnt %>% emtrends(~Treatment, var = "Days")  %>% multcomp::cld(Letters = letters, reversed = TRUE)
PWUE.lm_cnt %>% emtrends(~Treatment, var = "Days")  %>% pwpm

PWUE.lm_dcrt <- lm(A/gsw ~ Treatment*factor(Days), 
                   # SSlogis(Days, Asym, xmid, scal),
                   # params = list(Asym + scal + xmid ~ Treatment),
                   data = Photos.df %>% filter(gsw > 0 & gsw < .3))
PWUE.lm_dcrt %>% anova

### ___ ii) Plot ----
PWUE.emm <- 
  PWUE.lm_dcrt %>% emmeans(~Treatment|Days) %>% multcomp::cld(Letters = letters, reversed = TRUE) %>% 
  group_by(Days) %>% 
  mutate(Signif = ifelse(all(.group == " a"), "", "*"),
         ypos = min(emmean - SE))

wi.trt.lttr <- PWUE.lm_dcrt %>% emmeans(~Days|Treatment) %>% 
  multcomp::cld(Letters = letters, reversed = TRUE) %>% as.data.frame() %>% 
  group_by(Days) %>% 
  mutate(SE = ifelse(Days == 0, 0, SE),
         ypos = ifelse(emmean == max(emmean), emmean + SE + SE*.5, emmean - SE - SE*.5),
         .group = if(all(.group == " a")) {""} else {gsub(" ", "", .group)}) %>% 
  as.data.frame() %>% 
  left_join(
    data.frame(Days = Photos.df %>% arrange(Days) %>% .$Days %>% unique,
               Treatment = rep(c("Cold", "Room"), each = 20),
               ypos2 = c(seq(0, -8, length.out =20),
                         seq(-2, -10, length.out =20))
    )
  ) 

# PWUE.plt <-

PWUE.emm %>% as.data.frame() %>% 
  mutate(SE = ifelse(Days == 14, 0, SE)) %>% 
  ggplot(aes(x = Days, y = emmean, col = Treatment)) +
  scale_color_manual(values = TrtCol) +
  
  geom_point(position = position_dodge(0.6), show.legend = T, size = 2) +
  
  geom_errorbar(aes(ymin = emmean  - SE, ymax = emmean  + SE),
                width = 5, show.legend = F, position = position_dodge(0.6)) +
  
  # ## Add Trt diff
  # geom_text(data = . %>% filter(Signif == "*"), aes(y = 0, label =Signif), 
  #           col = "grey35", show.legend = F, cex = 8) +
  
  ## Add change diff
  # geom_text(data = wi.trt.lttr, aes(y = ypos, label = .group), show.legend = F) +
  
  # geom_smooth(method = "lm",
  #             formula = y ~ poly(x, 3), se = F, lwd = .5) +
  
geom_smooth(method = "lm",  se = F, lwd = .5) +
  
  scale_x_continuous(name = "Days since planting", limits = c(14 , 150)) +
  # scale_y_continuous(name = "gsw", limits = c(0, .3)) +
  myTheme + theme(legend.position = c(.9, .9))
PWUE.plt

## ---- IV. phiPSII ----

### ___ i) Model fitting ----

#### Linear
PhiPS2.lm_cnt <- lm(PhiPS2 ~ Days*Treatment,
                    # SSlogis(Days, Asym, xmid, scal),
                    # params = list(Asym + scal + xmid ~ Treatment),
                    data = Photos.df %>% filter(PhiPS2 > 0))
PhiPS2.lm_cnt %>% anova
PhiPS2.lm_cnt %>% summary
PhiPS2.lm_cnt %>% emtrends(~Treatment, var = "Days")  %>% multcomp::cld(Letters = letters, reversed = TRUE)
PhiPS2.lm_cnt %>% emtrends(~Treatment, var = "Days")  %>% pwpm

PhiPS2.lm_dcrt <- lm(PhiPS2 ~ Treatment*factor(Days),
                     data = Photos.df %>% filter(gsw > 0 & gsw < .3))

PhiPS2.lm_dcrt %>% anova

### ___ ii) Plot ----
PhiPS2.emm <- 
  PhiPS2.lm_dcrt %>% emmeans(~Treatment|Days) %>% multcomp::cld(Letters = letters, reversed = TRUE) %>% 
  group_by(Days) %>% 
  mutate(Signif = ifelse(all(.group == " a"), "", "*"),
         ypos = min(emmean - SE))

wi.trt.lttr <- PhiPS2.lm_dcrt %>% emmeans(~Days|Treatment) %>% 
  multcomp::cld(Letters = letters, reversed = TRUE) %>% as.data.frame() %>% 
  group_by(Days) %>% 
  mutate(SE = ifelse(Days == 0, 0, SE),
         ypos = ifelse(emmean == max(emmean), emmean + SE + SE*.5, emmean - SE - SE*.5),
         .group = if(all(.group == " a")) {""} else {gsub(" ", "", .group)}) %>% 
  as.data.frame() 

PhiPS2.plt <- PhiPS2.emm %>% as.data.frame() %>% 
  mutate(SE = ifelse(Days == 14, 0, SE)) %>% 
  ggplot(aes(x = Days, y = emmean, col = Treatment)) +
  scale_color_manual(values = TrtCol) +
  
  geom_point(position = position_dodge(0.6), show.legend = T, size = 2) +
  
  geom_errorbar(aes(ymin = emmean  - SE, ymax = emmean  + SE),
                width = 5, show.legend = F, position = position_dodge(0.6)) +
  
  ## Add Trt diff
  geom_text(data = . %>% filter(Signif == "*"), aes(y = 0, label =Signif), 
            col = "grey35", show.legend = F, cex = 8) +
  
  ## Add change diff
  # geom_text(data = wi.trt.lttr, aes(y = ypos, label = .group), show.legend = F) +
  
  # geom_smooth(method = "lm",
  #             formula = y ~ poly(x, 3), se = F, lwd = .5) +
  
  geom_smooth(method = "lm",  se = F, lwd = .5) +
  
  scale_x_continuous(name = "Days since planting", limits = c(0, 165)) +
  scale_y_continuous(name = "phiPSII", limits = c(0, .6)) +
  myTheme + theme(legend.position = c(.9, .9))
PhiPS2.plt



###                         ###
### ---- CARBOHYDRATES ----  
###                         ###

## ---- I. Rhizome Starch (%) ----
carbRhi.df <- Allvar.df %>%  
  filter(!(is.na(Starch_Rhi))) %>% 
  select(RhiID, Treatment, Sampling, HarvestDate, Days, Caputd, Phase, Starch_Rhi, Sucrose_Rhi, Glucose_Rhi) %>% 
  # filter(Root_g < 17) %>% 
  
  filter(Caputd == "alive") %>%
  filter(Days > 0)

### ___ i) Model fitting ----

#### Linear
starchRhi_lm_cnt <- lm(Starch_Rhi ~ Days*Treatment,
                       # SSlogis(Days, Asym, xmid, scal),
                       # params = list(Asym + scal + xmid ~ Treatment),
                       data = carbRhi.df )
starchRhi_lm_cnt %>% anova
starchRhi_lm_cnt %>% summary
starchRhi_lm_cnt %>% emtrends(~Treatment, var = "Days")  %>% multcomp::cld(Letters = letters, reversed = TRUE)

starchRhi_lm_dcrt <- lm(Starch_Rhi ~ Treatment*factor(Days), 
                       # SSlogis(Days, Asym, xmid, scal),
                       # params = list(Asym + scal + xmid ~ Treatment),
                       data = carbRhi.df )
starchRhi_lm_dcrt %>% anova
starchRhi.emm <- starchRhi_lm_dcrt %>% emmeans(~Treatment|Days) %>% multcomp::cld(Letters = letters, reversed = TRUE)

starchRhi_lm_dcrt %>% emmeans(~Treatment|Days) %>% pwpm
  ## Day 161: 14.8/12.9 - 1 (P = 0.5232)
  ## Day 161: 4.03/2.32 - 1 (P = 0.404)

wi.trt.lttr <- starchRhi_lm_dcrt %>% emmeans(~Days|Treatment) %>% 
  multcomp::cld(Letters = letters, reversed = TRUE) %>% as.data.frame() %>% 
  group_by(Days) %>% 
  mutate(SE = ifelse(Days == 0, 0, SE),
         ypos = ifelse(emmean == max(emmean), emmean + SE + SE*.5, emmean - SE - SE*.5),
         ypos = ifelse(Days == 0, 1, ypos),
         .group = if(all(.group == " a")) {""} else {gsub(" ", "", .group)}) %>% 
  as.data.frame()

### ___ ii) Plot ----
starchRhi.plt <- starchRhi.emm %>% as.data.frame() %>% 
  ggplot(aes(x = Days, y = emmean, col = Treatment)) +
  scale_color_manual(values = TrtCol) +
  
  geom_point(position = position_dodge(0.6), show.legend = T, size = 2) +
  
  geom_errorbar(aes(ymin = emmean  - SE, ymax = emmean  + SE),
                width = 5, show.legend = F, position = position_dodge(0.6)) +
  
  ## Add change diff
  geom_text(data = wi.trt.lttr, aes(y = ypos, label = .group), show.legend = F) +
  
  geom_smooth(method = "nls", data = . %>% filter(Treatment == 'Cold'),
              formula = y ~ Asym/(1+exp((xmid-x)/scal)) + c , 
              method.args = list(start = list(Asym = 12, xmid = 120, scal = 2, c = 4)), se = F, lwd = .5) +
  
  geom_smooth(method = "lm", data = . %>% filter(Treatment == 'Room'), se = F, lwd = .5) +
  
  scale_x_continuous(name = "Days since planting", limits = c(0, 165)) +
  scale_y_continuous(name = "Rhizome Starch (%)", limits = c(-0, 18), expand = c(0,0)) +
  myTheme + theme(legend.position = c(.9, .1))

#### nonlinear - not to be included
starchRhi_cold_c <- gnls(Starch_Rhi ~ Asym/(1+exp((xmid-Days)/scal)) + c,
                       # SSlogis(Days, Asym, xmid, scal),
                    # params = list(Asym + scal + xmid ~ Treatment),
                    data = carbRhi.df %>% filter(Treatment == "Cold"),
                    start = list(Asym = 12, xmid = 120, scal = 2, c = 4))
anova(starchRhi_cold_c)
summary(starchRhi_cold_c)

starchRhi_cold_Noc <- gnls(Starch_Rhi ~ SSlogis(Days, Asym, xmid, scal),
                         # SSlogis(Days, Asym, xmid, scal),
                         # params = list(Asym + scal + xmid ~ Treatment),
                         data = carbRhi.df %>% filter(Treatment == "Cold")
                         )
anova(starchRhi_cold_Noc)
summary(starchRhi_cold_Noc)

anova(starchRhi_cold_Noc, starchRhi_cold_c)

starchRhi_room_logis <- gnls(Starch_Rhi ~ SSlogis(Days, Asym, xmid, scal),
                           # SSlogis(Days, Asym, xmid, scal),
                           # params = list(Asym + scal + xmid ~ Treatment),
                           data = carbRhi.df %>% filter(Treatment == "Room")
)
summary(starchRhi_room_logis)

# starchRhi_room_logis_c <- gnls(Starch_Rhi ~ Asym/(1+exp((xmid-Days)/scal)) + c,
#                                # SSlogis(Days, Asym, xmid, scal),
#                                # params = list(Asym + scal + xmid ~ Treatment),
#                                data = carbRhi.df %>% filter(Treatment == "Room"),
#                                start = list(Asym = 17, xmid = 125, scal = 45, c = 2))

# starchRhi_room_asymOff <- gnls(Starch_Rhi ~ SSasympOff(Days, Asym, lrc, c0), 
#                              # SSlogis(Days, Asym, xmid, scal),
#                              # params = list(Asym + scal + xmid ~ Treatment),
#                              data = carbRhi.df %>% filter(Treatment == "Room")
# )
# 
# starchRhi_lm <- lm(Starch_Rhi ~ Days*Treatment,
#                              # SSlogis(Days, Asym, xmid, scal),
#                              # params = list(Asym + scal + xmid ~ Treatment),
#                              data = carbRhi.df %>% filter(Treatment == "Room") )
starchRhi_room_lm <-
  lm(Starch_Rhi ~ Days,
     data = carbRhi.df %>% filter(Treatment == "Room") )

starchRhi_full %>% plot
starchRhi_lm %>% anova

## ---- II. Rhizome Starch (g) ----
carbRhi.df <- Allvar.df %>%  
  filter(!(is.na(Starch_Rhi))) %>% 
  mutate(Starch_Rhi_g = Starch_Rhi*RhiSample_g) %>% 
  select(RhiID, Treatment, Sampling, HarvestDate, Days, Caputd, Phase, 
         RhiSample_g, Starch_Rhi, Sucrose_Rhi, Glucose_Rhi,
         Starch_Rhi_g) %>% 
  # filter(Root_g < 17) %>% 
  
  filter(!(is.na(Starch_Rhi_g))) %>%
  filter(Caputd == "alive" & HarvestDate >= "2021-04-25") %>%
  filter(Days > 14)

### ___ i) nonlinear fit ----
# starchRhi_g_cold_c <- gnls(Starch_Rhi_g ~ Asym/(1+exp((xmid-Days)/scal)) + c ,
#                            start = list(Asym = 240, xmid = 126, scal = 25, c = 2),
#                            data = carbRhi.df %>% filter(Treatment == "Cold"))
# summary(starchRhi_g_cold_c)
# anova(starchRhi_g_cold_c)
## ^ no better fit
# 
# starchRhi_g_cold_Noc <- gnls(Starch_Rhi_g ~ SSlogis(Days, Asym, xmid, scal),
#                            # SSlogis(Days, Asym, xmid, scal),
#                            # params = list(Asym + scal + xmid ~ Treatment),
#                            data = carbRhi.df %>% filter(Treatment == "Cold")
# )
# anova(starchRhi_cold_Noc)
# summary(starchRhi_cold_Noc)
# 
# anova(starchRhi_cold_Noc, starchRhi_g_cold_c)
# 
# starchRhi_g_room_Noc <- gnls(Starch_Rhi_g ~ SSlogis(Days, Asym, xmid, scal), 
#                            # SSlogis(Days, Asym, xmid, scal),
#                            # params = list(Asym + scal + xmid ~ Treatment),
#                            data = carbRhi.df %>% filter(Treatment == "Room")
# )
# summary(starchRhi_g_room_Noc)

starchRhi_g_null <- gnls(Starch_Rhi_g ~ SSlogis(Days, Asym, xmid, scal),
                         data = carbRhi.df)


starchRhi_g_full <- gnls(Starch_Rhi_g ~ Asym/(1+exp((xmid-Days)/scal)) ,
                           params = list(Asym  ~ Treatment,   xmid + scal ~ 1),
                           start = list(Asym = c(236.30925, -50.04012), xmid = rep(120, 1), scal = rep(25, 1)),
                           data = carbRhi.df)
summary(starchRhi_g_full)

anova(starchRhi_g_null, starchRhi_g_full)

starchRhi_g_full %>% emmeans(param = "Asym", specs = ~Treatment) %>% multcomp::cld(Letters = letters, reversed = TRUE) 
starchRhi_g_full %>% emmeans(param = "Asym", specs = ~Treatment) %>% pwpm

starchRhi_g_full %>% emmeans(param = "xmid", specs = ~1) %>% multcomp::cld(Letters = letters, reversed = TRUE) 
starchRhi_g_full %>% emmeans(param = "scal", specs = ~1) %>% multcomp::cld(Letters = letters, reversed = TRUE) 

starchRhi_g_lm <- lm(Starch_Rhi_g ~ Days*Treatment,
                   data = carbRhi.df )
starchRhi_g_lm %>% anova

AIC(starchRhi_g_lm, starchRhi_g_full)
BIC(starchRhi_g_lm, starchRhi_g_full)


### ___ ii) Plot ----
starchRhi_g_lm_dcrt <- lm(Starch_Rhi_g ~ factor(Days)*Treatment,
                     data = carbRhi.df )
starchRhi_g_lm_dcrt %>% anova
starchRhi_g.emm <- starchRhi_g_lm_dcrt %>% emmeans(~Treatment|Days) %>% multcomp::cld(Letters = letters, reversed = TRUE)
starchRhi_g_lm_dcrt %>% emmeans(~Treatment|Days) %>% pwpm
    ## Day 161: 210/162 -1 = 29% (P = 0.322)

wi.trt.lttr <- starchRhi_g_lm_dcrt %>% emmeans(~Days|Treatment) %>% 
  multcomp::cld(Letters = letters, reversed = TRUE) %>% as.data.frame() %>% 
  group_by(Days) %>% 
  mutate(SE = ifelse(Days == 0, 0, SE),
         ypos = ifelse(emmean == max(emmean), emmean + SE + SE*.5, emmean - SE - SE*.5),
         ypos = ifelse(Days == 0, 1, ypos),
         .group = if(all(.group == " a")) {""} else {gsub(" ", "", .group)}) %>% 
  as.data.frame()

starchRhi_g.plt <- starchRhi_g.emm %>% as.data.frame() %>% 
  ggplot(aes(x = Days, y = emmean, col = Treatment)) +
  scale_color_manual(values = TrtCol) +
  
  geom_point(position = position_dodge(0.6), show.legend = T, size = 2) +
  
  geom_errorbar(aes(ymin = emmean  - SE, ymax = emmean  + SE),
                width = 5, show.legend = F, position = position_dodge(0.6)) +
  
  ## Add change diff
  geom_text(data = wi.trt.lttr, aes(y = ypos, label = .group), show.legend = F) +
  
  geom_smooth(method = "nls",
              formula = y ~ SSlogis(x, Asym, xmid, scal), se = F, lwd = .5) +
  
  # geom_smooth(method = "lm", data = . %>% filter(Treatment == 'Room'), se = F, lwd = .5) +
  
  scale_x_continuous(name = "Days since planting", limits = c(0, 165)) +
  scale_y_continuous(name = "Rhizome Starch (g)", limits = c(-15, 270), expand = c(0,0)) +
  myTheme + theme(legend.position = c(.9, .1))

## ---- III. Rhizome Sucrose (%) ----
carbRhi.df <- Allvar.df %>%  
  filter(!(is.na(Starch_Rhi))) %>% 
  select(RhiID, Treatment, Sampling, HarvestDate, Days, Caputd, Phase, Starch_Rhi, Sucrose_Rhi, Glucose_Rhi) %>% 
  # filter(Root_g < 17) %>% 
  
  filter(Caputd == "alive") %>%
  filter(Days > 0)

### ___ i) Model fitting ----

#### Linear
sucroseRhi_lm_cnt <- lm(Sucrose_Rhi ~ Days*Treatment,
                       # SSlogis(Days, Asym, xmid, scal),
                       # params = list(Asym + scal + xmid ~ Treatment),
                       data = carbRhi.df )
sucroseRhi_lm_cnt %>% anova
sucroseRhi_lm_cnt %>% summary
sucroseRhi_lm_cnt %>% emtrends(~Treatment, var = "Days")  %>% multcomp::cld(Letters = letters, reversed = TRUE)

sucroseRhi_lm_dcrt <- lm(Sucrose_Rhi ~ Treatment*factor(Days), 
                        # SSlogis(Days, Asym, xmid, scal),
                        # params = list(Asym + scal + xmid ~ Treatment),
                        data = carbRhi.df )
sucroseRhi_lm_dcrt %>% anova
sucroseRhi.emm <- sucroseRhi_lm_dcrt %>% emmeans(~Treatment|Days) %>% multcomp::cld(Letters = letters, reversed = TRUE)

sucroseRhi_lm_dcrt %>% emmeans(~Treatment|Days) %>% pwpm
## Day 161: 14.8/12.9 - 1 (P = 0.5232)
## Day 161: 4.03/2.32 - 1 (P = 0.404)

wi.trt.lttr <- sucroseRhi_lm_dcrt %>% emmeans(~Days|Treatment) %>% 
  multcomp::cld(Letters = letters, reversed = TRUE) %>% as.data.frame() %>% 
  group_by(Days) %>% 
  mutate(SE = ifelse(Days == 0, 0, SE),
         ypos = ifelse(emmean == max(emmean), emmean + SE + SE*.5, emmean - SE - SE*.5),
         ypos = ifelse(Days == 0, 1, ypos),
         .group = if(all(.group == " a")) {""} else {gsub(" ", "", .group)}) %>% 
  as.data.frame()

### ___ ii) Plot ----
sucroseRhi.plt <- sucroseRhi.emm %>% as.data.frame() %>% 
  ggplot(aes(x = Days, y = emmean, col = Treatment)) +
  scale_color_manual(values = TrtCol) +
  
  geom_point(position = position_dodge(0.6), show.legend = T, size = 2) +
  
  geom_errorbar(aes(ymin = emmean  - SE, ymax = emmean  + SE),
                width = 5, show.legend = F, position = position_dodge(0.6)) +
  
  ## Add change diff
  geom_text(data = wi.trt.lttr, aes(y = ypos, label = .group), show.legend = F) +
  geom_line() +
  # geom_smooth(method = "nls", data = . %>% filter(Treatment == 'Cold'),
  #             formula = y ~ Asym/(1+exp((xmid-x)/scal)) + c , 
  #             method.args = list(start = list(Asym = 12, xmid = 120, scal = 2, c = 4)), se = F, lwd = .5) +
  # 
  # geom_smooth(method = "lm", se = F, lwd = .5) +
  
  scale_x_continuous(name = "Days since planting", limits = c(0, 165)) +
  scale_y_continuous(name = "Rhizome sucrose (%)", limits = c(-0, 5), expand = c(0,0)) +
  myTheme + theme(legend.position = c(.9, .1))
sucroseRhi.plt

## ---- IV. Rhizome Glucose (%) ----
carbRhi.df <- Allvar.df %>%  
  filter(!(is.na(Starch_Rhi))) %>% 
  select(RhiID, Treatment, Sampling, HarvestDate, Days, Caputd, Phase, Starch_Rhi, Sucrose_Rhi, Glucose_Rhi) %>% 
  # filter(Root_g < 17) %>% 
  
  filter(Caputd == "alive") %>%
  filter(Days >= 0)

### ___ i) Model fitting ----

#### Linear
glucoseRhi_lm_cnt <- lm(Glucose_Rhi ~ Days*Treatment,
                        # SSlogis(Days, Asym, xmid, scal),
                        # params = list(Asym + scal + xmid ~ Treatment),
                        data = carbRhi.df )
glucoseRhi_lm_cnt %>% anova
glucoseRhi_lm_cnt %>% summary
glucoseRhi_lm_cnt %>% emtrends(~Treatment, var = "Days")  %>% multcomp::cld(Letters = letters, reversed = TRUE)

glucoseRhi_lm_dcrt <- lm(Glucose_Rhi ~ Treatment*factor(Days), 
                         # SSlogis(Days, Asym, xmid, scal),
                         # params = list(Asym + scal + xmid ~ Treatment),
                         data = carbRhi.df )
glucoseRhi_lm_dcrt %>% anova
glucoseRhi.emm <- glucoseRhi_lm_dcrt %>% emmeans(~Treatment|Days) %>% multcomp::cld(Letters = letters, reversed = TRUE)

glucoseRhi_lm_dcrt %>% emmeans(~Treatment|Days) %>% pwpm
## Day 161: 14.8/12.9 - 1 (P = 0.5232)
## Day 161: 4.03/2.32 - 1 (P = 0.404)

wi.trt.lttr <- glucoseRhi_lm_dcrt %>% emmeans(~Days|Treatment) %>% 
  multcomp::cld(Letters = letters, reversed = TRUE) %>% as.data.frame() %>% 
  group_by(Days) %>% 
  mutate(SE = ifelse(Days == 0, 0, SE),
         ypos = ifelse(emmean == max(emmean), emmean + SE + SE*.5, emmean - SE - SE*.5),
         ypos = ifelse(Days == 0, 1, ypos),
         .group = if(all(.group == " a")) {""} else {gsub(" ", "", .group)}) %>% 
  as.data.frame()

### ___ ii) Plot ----
glucoseRhi.plt <- glucoseRhi.emm %>% as.data.frame() %>% 
  ggplot(aes(x = Days, y = emmean, col = Treatment)) +
  scale_color_manual(values = TrtCol) +
  
  geom_point(position = position_dodge(0.6), show.legend = T, size = 2) +
  
  geom_errorbar(aes(ymin = emmean  - SE, ymax = emmean  + SE),
                width = 5, show.legend = F, position = position_dodge(0.6)) +
  
  ## Add change diff
  geom_text(data = wi.trt.lttr, aes(y = ypos, label = .group), show.legend = F) +
  geom_line() +
  # geom_smooth(method = "nls", data = . %>% filter(Treatment == 'Cold'),
  #             formula = y ~ Asym/(1+exp((xmid-x)/scal)) + c , 
  #             method.args = list(start = list(Asym = 12, xmid = 120, scal = 2, c = 4)), se = F, lwd = .5) +
  # 
  # geom_smooth(method = "lm", se = F, lwd = .5) +
  
  scale_x_continuous(name = "Days since planting", limits = c(-2, 165)) +
  scale_y_continuous(name = "Rhizome glucose (%)", limits = c(-0, 5), expand = c(0,0)) +
  myTheme + theme(legend.position = c(.9, .1))
glucoseRhi.plt

FigS5_new <- cowplot::plot_grid(starchRhi.plt + scale_x_continuous(name = NULL, labels = NULL),
                                sucroseRhi.plt + scale_x_continuous(name = NULL, labels = NULL),  
                                glucoseRhi.plt, ncol = 1,
                                align = "hv", labels = "AUTO")
FigS5_new

###                     ###
### ---- GROWTH RATES ----  
###                     ###


## i) Predict values for each day

### Create new dataset for predict()
DaysTrt.df <- expand.grid(Days = 1:180,
                          Treatment = c("Cold", "Room"))

### Estimate predictions
AllPreds.df <- cbind(DaysTrt.df, 
                     belowBM = predict(belowBM_full, newdata = DaysTrt.df),
                     rhiBM = predict(rhiBM_full, newdata = DaysTrt.df),
                     rootBM = predict(rootBM_full, newdata = DaysTrt.df),
                     aboveBM = predict(AboveBM_full, newdata = DaysTrt.df),
                     tiller = predict(Tiller_full, newdata = DaysTrt.df),
                     LeafDev = predict(Leaf_full, newdata = DaysTrt.df),
                     Anet = predict(Anet_lm_cnt, newdata = DaysTrt.df),
                     Starch.p = c(predict(starchRhi_cold_c, newdata = DaysTrt.df %>% filter(Treatment == "Cold")),
                                  predict(starchRhi_room_lm, newdata = DaysTrt.df %>% filter(Treatment == "Room")))
) 

### Long format makes it easier to process
AllPreds.lng.df <- AllPreds.df %>% 
  tidyr::pivot_longer(cols = belowBM:Starch.p, 
                      names_to = "Variable", values_to = "Value")

## ii) Estimate growth and development rates
AllRates.lng.df <- 
  AllPreds.lng.df %>% 
  group_by(Variable, Treatment) %>% 
  mutate(
    Rate = c(0, diff(Value)/1),
      ## This step estimates the rate between consecutive dates (divide by 1)
    
    # CIValue = ifelse(Value >= max(Value)*.1 & Value <= max(Value)*.9, Value, NA),
    # CIRate = c(0, diff(CIValue)),
  )

## iii) Estimate start and end of the stage
StartEnd.Day.df <-
  AllPreds.lng.df %>% 
  filter(!(Variable %in% c("Anet", "belowBM"))) %>% 
  group_by(Variable, Treatment) %>% 
  summarize(Start.Day = min(Days[Value >= max(Value)*.15]),
            End.Day = max(Days[Value <= max(Value)*.85])) %>% 
  mutate(Variable = factor(Variable, levels = c("aboveBM", "rhiBM", "rootBM", "tiller", "LeafDev", "Starch.p")),
         yPos = ifelse(Treatment == "Room", -.1, -.2))

StartEnd.Day.df.lng <-
  StartEnd.Day.df%>% 
  select(-yPos) %>% 
  tidyr::pivot_longer(cols = Start.Day:End.Day)

StartEnd.Day.df.lng %>%
  ggplot(aes(x = Variable, y= value, col = Treatment, fill = Treatment)) +
  # facet_grid(~Treatment) +
  scale_color_manual(values = TrtCol) +
  scale_fill_manual(values = TrtCol) +

  geom_line(position = position_dodge(.5)) +  
  coord_flip() +
  myTheme

  
## iv) plot
Rates.plt <-
  AllRates.lng.df %>% 
  filter(Rate > 0) %>% 
  filter(!(Variable %in% c("Anet", "belowBM"))) %>% 
  
  mutate(Variable = factor(Variable, levels = c("aboveBM", "rhiBM", "rootBM", "tiller", "LeafDev", "Starch.p"))) %>%
  group_by(Variable) %>% 
  mutate(relRate = Rate/max(Rate)) %>% 
  
  ggplot(aes(x = Days, col = Treatment, fill = Treatment)) +
  scale_color_manual(values = TrtCol) +
  scale_fill_manual(values = TrtCol) +
  
  facet_grid(Variable ~., scales = "free_y", switch = "y") +
  # geom_polygon(aes(y = relRate), alpha = 1, show.legend = F) +
  geom_line(aes( y = relRate), lwd = 1) +
  
  geom_hline(yintercept = 0, lwd = .5) +
  
  ## Add rectangle with Start & End  
  geom_rect(data = StartEnd.Day.df, inherit.aes = F, show.legend = F,
            aes(ymin = yPos %>% as.numeric  - .05, ymax = yPos %>% as.numeric + .05, 
                xmin = Start.Day, xmax = End.Day,
                col = Treatment, fill = Treatment), col = "black"
            ) +

  scale_x_continuous(expand = c(0,0), breaks = c(0, 50, 100, 150), limits = c(0, 180))+
  scale_y_continuous(name = NULL, breaks = NULL, limits = c(-.30, 1.1), expand = c(0,0), ) +
  myTheme + theme(strip.background = element_blank(),
                  strip.text = element_text(size = 14),
                  # panel.background = element_blank(),
                  # plot.background = element_rect(colour = "black"),
                  legend.position = c(.9, .95),
                  # axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')
                  )


  
  Rates.plt

###### bottstraping
require(nlraa)

fm1.P.dm <- predict2_nls(AboveBM_full, interval = "conf")
## Reproducing book results
predict2_nls(starchRhi_cold_c, 
             interval = "conf", 
             newdata = data.frame(Treatment = "Cold", Days = 180))


## difference between trt over time
diff.ch_AllRates.lng.df <- AllRates.lng.df %>% 
  filter(Rate > 0) %>% 
  filter(!(Variable %in% c("Anet", "belowBM"))) %>% 
  
  mutate(Variable = factor(Variable, levels = c("aboveBM", "rhiBM", "rootBM", "tiller", "LeafDev", "Starch.p"))) %>%
  select(-Value) %>% 
  tidyr::pivot_wider(values_from = Rate, names_from = Treatment) %>% 
  
  mutate(chTrt = Room/Cold,
         diffTrt = Room - Cold) 

### relative importance of rates by treatment and day
stackedRates.plt <-
  AllRates.lng.df %>% 
  filter(Rate > 0) %>% 
  filter(!(Variable %in% c("Anet", "belowBM"))) %>% 
  
  # mutate(Variable = factor(Variable, levels = c("aboveBM", "rhiBM", "rootBM", "tiller", "LeafDev", "Starch.p"))) %>%
  group_by(Days, Treatment) %>% 
  mutate(relAlloc = Rate/sum(Rate)) %>% 
  
  # filter(Treatment == "Cold") %>% 
  
  ggplot(aes(x = Days, y = Rate, fill = Variable)) +
  facet_grid(Treatment ~ .) +
  scale_fill_manual(values = c('#7FC97F','#BEAED4','#FDC086','#FFFF99', "#386CB0", '#F0027F')) +
  geom_area(data = . %>% filter(Variable != "Starch.p"),
            position = 'stack', col = "Black") +
  
  scale_x_continuous(name = "Days since planting",  expand = c(0,0),  
                     breaks = c(0, 50, 100, 150), limits = c(0, 180))+
  scale_y_continuous(name = "Relative sink activity", expand = c(0,0), ) +
    
  geom_line(data = . %>% filter(Variable == "Starch.p"), lwd = 1, col = "#386CB0") +
  
  myTheme + theme(legend.position = "top")

lineAllRates.plt <-
  AllRates.lng.df %>% 
  filter(Rate > 0) %>% 
  filter(!(Variable %in% c("Anet", "belowBM"))) %>% 
  
  # mutate(Variable = factor(Variable, levels = c("aboveBM", "rhiBM", "rootBM", "tiller", "LeafDev", "Starch.p"))) %>%
  group_by(Days, Treatment) %>% 
  mutate(relAlloc = Rate/sum(Rate)) %>% 
  
  # filter(Treatment == "Cold") %>% 
  
  ggplot(aes(x = Days, y = Rate, col = Variable)) +
  facet_grid(Treatment ~ .) +
  scale_color_manual(values = c('#7FC97F','#BEAED4','#FDC086','#FFFF99',"#386CB0", '#F0027F')) +
  # geom_area(data = . %>% filter(Variable != "Starch.p"),
  #           position = 'stack', col = "Black") +
  
  scale_x_continuous(name = "Days since planting",  expand = c(0,0),  
                     breaks = c(0, 50, 100, 150), limits = c(0, 180))+
  scale_y_continuous(name = "Relative sink activity", expand = c(0,0), ) +
  
  geom_line(lwd = 2) +
  
  myTheme + theme(legend.position = "top")

diff.ch_AllRates.lng.df %>% 
  ggplot(aes(x = Days, y = chTrt, fill = Variable)) +
  facet_grid(Variable~., scales = "free_y") +
  scale_fill_manual(values = c('#7FC97F','#BEAED4','#FDC086','#FFFF99','#386CB0','#F0027F')) +
  geom_area(position = 'stack', col = "Black") +
  
  geom_hline(yintercept = 1 )+
  
  scale_x_continuous(name = "Days since planting", 
                     breaks = c(0, 50, 100, 150), limits = c(0, 180))+
  scale_y_continuous(name = "Room - Cold", expand = c(0,0), ) +
  
  myTheme

  
