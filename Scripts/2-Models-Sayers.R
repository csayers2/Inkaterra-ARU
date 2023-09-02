
# Chris Sayers
# updated April 9, 2023

# script designed for data visualizations and model fitting

library(tidyverse)
library(ggplot2)
library(stringr)
library(ggpubr)
library(mgcv)
library(MuMIn)
library(ggview)

# SPECIES RICHNESS --------------------------------------------------------

SR.Total <- read.csv("Outputs/SR.Total")
SR.Window.60 <- read.csv("Outputs/SR.Window.60") %>% 
  mutate(Minute = Time.Window/60,
         fTime.Window = as.factor(Time.Window),
         Site = as.factor(Site),
         Day = as.factor(Day),
         Hab1 = as.factor(Hab1),
         SiteDay = as.factor(str_c(Site, Day))) 

# OUTLIERS & NORMALITY OF RESPONSE VARIABLES -----------------------------------
ggdensity(SR.Window.60$SR, xlab = "Species Richness")
gghistogram(SR.Window.60$SR, xlab = "Species Richness")
ggqqplot(SR.Window.60$SR, ylab = "Species Richness")
shapiro.test(SR.Window.60$SR) # W = 0.97624, p-value = 2.612e-12, not normal

ggplot(data = SR.Window.60) +
  geom_point(mapping = aes(x = Minute, y = SR, color = Day)) +
  geom_smooth(mapping = aes(x = Minute, y = SR, color = Day)) + 
  facet_grid(Day ~ Site)

# SPECIES RICHNESS MODEL -------------------------------------------------------

# biologically relevant GAM structure
SRgam.60 <- mgcv::gamm(SR ~ s(Minute, by = SiteDay) + s(Site, bs = "re") + s(Day, bs = "re"),
                       family = "poisson",
                       method = "REML",
                       data = SR.Window.60)

# checking for autocorrelation issues
par(mfrow=c(1,1))
performance::check_singularity(SRgam.60$gam)
acf(SR.Window.60$SR) # raw data is autocorrelated
acf(resid(SRgam.60$lme, type = "normalized")) # we have mild autocorrelation in the residuals
pacf(resid(SRgam.60$lme, type = "normalized"))

# adding in first order autoregressive covariance structure (AR1) to account for
# residual autocorrelation
SRgam.60.ar1 <- mgcv::gamm(SR ~ s(Minute, by = SiteDay) + s(Site, bs = "re") + s(Day, bs = "re"),
                             correlation = corAR1(form = ~ Minute | SiteDay),
                             family = "poisson",
                             method = "REML",
                             data = SR.Window.60)
# another reasonable model structure is to include site and day as FEs

# checking for autocorrelation issues
performance::check_singularity(SRgam.60.ar1$gam)
acf(resid(SRgam.60.ar1$lme, type = "normalized")) # visually, we do not have autocorrelated data
pacf(resid(SRgam.60.ar1$lme, type = "normalized")) # visually, we do not have autocorrelated data

# checking model diagnostics
par(mfrow=c(2,2))
gam.check(SRgam.60.ar1$gam) # k values are too small, but we can't change them, residuals look great
concurvity(SRgam.60.ar1$gam, full = TRUE) # no issues with concurvity
concurvity(SRgam.60.ar1$gam, full = FALSE) # no issues with concurvity

summary(SRgam.60.ar1$gam)
anova.gam(SRgam.60.ar1$gam)
# visualizing partial effects
plot(SRgam.60.ar1$gam, shade = TRUE, shift = coef(SRgam.60.ar1$gam)[1],
     trans = exp, pages = 1, all.terms = TRUE, rug = FALSE)

SR.anova <- aov(SR ~ SiteDay, data = SR.Window.60)
SR.tukey <- TukeyHSD(SR.anova)

library(agricolae)
SR.HSD <- HSD.test(SR.anova, trt = 'group')




# creating a second model to smooth minute by day effects
SRgam.60.ar1.day <- mgcv::gamm(SR ~ s(Minute, by = Day) + s(Site, bs = "re") + s(Day, bs = "re"),
                               correlation = corAR1(form = ~ Minute | SiteDay),
                               family = "poisson",
                               method = "REML",
                               data = SR.Window.60)
summary(SRgam.60.ar1.day$gam)
anova.gam(SRgam.60.ar1.day$gam)
# visualizing partial effects
plot(SRgam.60.ar1.day$gam, shade = TRUE, shift = coef(SRgam.60.ar1$gam)[1],
     trans = exp, pages = 1, all.terms = TRUE, rug = FALSE)


# MODEL PREDICTION --------------------------------------------------------

# predicting from the model
predicted.sr <- data.frame(predict(SRgam.60.ar1$gam, type = "response", se.fit = TRUE))
predicted.sr.day <- data.frame(predict(SRgam.60.ar1.day$gam, type = "response", se.fit = TRUE)) %>% 
  rename(fit.day = fit, se.fit.day = se.fit)
SR.Window.60 <- cbind(SR.Window.60, predicted.sr, predicted.sr.day) %>% 
  # converting site and day names to meaningful numbers
  mutate(Day = if_else(Day == "A", "Day 1",
                       if_else(Day == "B", "Day 2", "Day 3")),
         Site = if_else(Site == 1, "Site A",
                        if_else(Site == 2, "Site B", 
                                if_else(Site == 4, "Site C",
                                        if_else(Site == 5, "Site D", 
                                                if_else(Site == 6, "Site E", "Site F"))))))

write.csv(SR.Window.60, "Outputs/SR.Window.60.pred")

# predictions across days and sites
SR.siteday.plot <- ggplot(data = SR.Window.60) +
  geom_point(mapping = aes(x = Minute, y = SR, color = Site)) +
  geom_line(mapping = aes(x = Minute, y = fit + se.fit), linetype = 2) +
  geom_line(mapping = aes(x = Minute, y = fit), size = 1.1) +
  geom_line(mapping = aes(x = Minute, y = fit - se.fit), linetype = 2) +
  facet_grid(Site ~ Day) +
  theme_bw(base_size = 16) +
  labs(x = "Minute", y = "Species Richness") +
  scale_color_brewer(palette = "Dark2") +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(hjust = 0.5),
        axis.text.y = element_text(hjust = 1),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold"),
        strip.text.y = element_text(face = "bold"),
        panel.spacing = unit(1, "lines"),
        aspect.ratio = 0.8)

ggview(device = "jpeg", units = "in", dpi = 1200, width = 6, height = 12)

ggsave("Figures/Fig3A-SR.jpg", dpi = 1200, width = 6, height = 12)

# smoothing predictions by day
SR.day.plot <- ggplot(data = SR.Window.60) +
  #geom_point(mapping = aes(x = Minute, y = SR, color = Site)) +
  geom_line(mapping = aes(x = Minute, y = fit, color = Site)) +
  geom_smooth(mapping = aes(x = Minute, y = fit.day + se.fit.day), se = FALSE, color = "black", linetype = 2) +
  geom_smooth(mapping = aes(x = Minute, y = fit.day), color = "black", size = 1.1, se = FALSE) +
  geom_smooth(mapping = aes(x = Minute, y = fit.day - se.fit.day), se = FALSE, color = "black", linetype = 2) +
  facet_grid(~ Day) +
  theme_bw(base_size = 16) +
  labs(x = "Minute", y = "Species Richness") +
  scale_color_brewer(palette = "Dark2") +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(hjust = 0.5),
        axis.text.y = element_text(hjust = 1),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold"),
        strip.text.y = element_text(face = "bold"),
        panel.spacing = unit(1, "lines"),
        aspect.ratio = 0.8)

ggview(device = "jpeg", units = "in", dpi = 1200, width = 12, height = 4)

ggsave("Figures/Fig3B-SR.jpg", dpi = 1200, width = 12, height = 4)


# TOTAL VOCAL PREVALENCE --------------------------------------------------------

TVP.Window.60 <- read.csv("Outputs/TVP.Window.60") %>% 
  mutate(Minute = Time.Window/60,
         Site = as.factor(Site),
         Day = as.factor(Day),
         Hab1 = as.factor(Hab1),
         SiteDay = as.factor(str_c(Site, Day))) 

# OUTLIERS & NORMALITY OF RESPONSE VARIABLES -----------------------------------
ggdensity(TVP.Window.60$TVP, xlab = "Total Vocal Prevalence")
gghistogram(TVP.Window.60$TVP, xlab = "Total Vocal Prevalence")
ggqqplot(TVP.Window.60$TVP, ylab = "Total Vocal Prevalence")
shapiro.test(TVP.Window.60$TVP) # W = 0.98492, p-value = 3.958e-09, not normal

ggplot(data = TVP.Window.60) +
  geom_point(mapping = aes(x = Minute, y = TVP, color = Day)) +
  geom_smooth(mapping = aes(x = Minute, y = TVP, color = Day)) + 
  facet_grid(Day ~ Site)

# TOTAL VOCAL PREVALENCE MODEL -------------------------------------------------

TVPgam.60 <- mgcv::gamm(TVP ~ s(Minute, by = SiteDay) + s(Site, bs = "re") + s(Day, bs = "re"),
                        family = "poisson",
                        method = "REML",
                        data = TVP.Window.60)

# checking for autocorrelation issues
par(mfrow=c(1,1))
performance::check_singularity(TVPgam.60$gam)
acf(resid(TVPgam.60$gam)) # we have strong autocorrelation in our residuals
acf(resid(TVPgam.60$lme, type = "normalized")) # we have strong autocorrelation in our residuals
pacf(resid(TVPgam.60$lme, type = "normalized"))

# adding in first order autoregressive covariance structure (AR1) to account for
# residual autocorrelation
TVPgam.60.ar1 <- mgcv::gamm(TVP ~ s(Minute, by = SiteDay) + s(Day, bs = "re") + s(Site, bs = "re"),
                        correlation = corAR1(form = ~ Minute | SiteDay),
                        family = "poisson",
                        method = "REML",
                        data = TVP.Window.60)

# checking for autocorrelation issues
performance::check_singularity(TVPgam.60.ar1$gam)
acf(resid(TVPgam.60.ar1$lme, type = "normalized")) # visually, we do not have autocorrelated residuals
pacf(resid(TVPgam.60.ar1$lme, type = "normalized"))

# checking model diagnostics
par(mfrow=c(2,2))
gam.check(TVPgam.60.ar1$gam) # k values are too small, but we can't change them, residuals look great
concurvity(TVPgam.60.ar1$gam, full = TRUE) # no issues with concurvity
concurvity(TVPgam.60.ar1$gam, full = FALSE) # no issues with concurvity

summary(TVPgam.60.ar1$gam)
anova.gam(TVPgam.60.ar1$gam)
# visualizing partial effects
plot(TVPgam.60.ar1$gam, shade = TRUE, shift = coef(TVPgam.60.ar1$gam)[1],
     trans = exp, pages = 1, all.terms = TRUE, rug = FALSE)

# creating a second model to smooth minute by day effects
TVPgam.60.ar1.day <- mgcv::gamm(TVP ~ s(Minute, by = Day) + s(Site, bs = "re") + s(Day, bs = "re"),
                               correlation = corAR1(form = ~ Minute | SiteDay),
                               family = "poisson",
                               method = "REML",
                               data = TVP.Window.60)
summary(TVPgam.60.ar1.day$gam)
anova.gam(TVPgam.60.ar1.day$gam)
# visualizing partial effects
plot(TVPgam.60.ar1.day$gam, shade = TRUE, shift = coef(TVPgam.60.ar1$gam)[1],
     trans = exp, pages = 1, all.terms = TRUE, rug = FALSE)


# MODEL PREDICTION --------------------------------------------------------

# predicting from the model
predicted.tvp <- data.frame(predict(TVPgam.60.ar1$gam, type = "response", se.fit = TRUE))
predicted.tvp.day <- data.frame(predict(TVPgam.60.ar1.day$gam, type = "response", se.fit = TRUE)) %>% 
  rename(fit.day = fit, se.fit.day = se.fit)
TVP.Window.60 <- cbind(TVP.Window.60, predicted.tvp, predicted.tvp.day) %>% 
  # converting site and day names to meaningful numbers
  mutate(Day = if_else(Day == "A", "Day 1",
                       if_else(Day == "B", "Day 2", "Day 3")),
         Site = if_else(Site == 1, "Site A",
                        if_else(Site == 2, "Site B", 
                                if_else(Site == 4, "Site C",
                                        if_else(Site == 5, "Site D", 
                                                if_else(Site == 6, "Site E", "Site F"))))))

write.csv(TVP.Window.60, "Outputs/TVP.Window.60.pred")

# predictions across days and sites
ggplot(data = TVP.Window.60) +
  geom_point(mapping = aes(x = Minute, y = TVP, color = Site)) +
  geom_line(mapping = aes(x = Minute, y = fit + se.fit), linetype = 2) +
  geom_line(mapping = aes(x = Minute, y = fit), size = 1.1) +
  geom_line(mapping = aes(x = Minute, y = fit - se.fit), linetype = 2) +
  facet_grid(Site ~ Day) +
  theme_bw(base_size = 16) +
  labs(x = "Minute", y = "Total Vocal Prevalence") +
  scale_color_brewer(palette = "Dark2") +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(hjust = 0.5),
        axis.text.y = element_text(hjust = 1),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold"),
        strip.text.y = element_text(face = "bold"),
        panel.spacing = unit(1, "lines"),
        aspect.ratio = 0.8)

ggview(device = "jpeg", units = "in", dpi = 1200, width = 6, height = 12)

ggsave("Figures/Fig3A-TVP.jpg", dpi = 1200, width = 6, height = 12)

# smoothing predictions by day
ggplot(data = TVP.Window.60) +
  #geom_point(mapping = aes(x = Minute, y = TVP, color = Site)) +
  geom_line(mapping = aes(x = Minute, y = fit, color = Site)) +
  geom_smooth(mapping = aes(x = Minute, y = fit.day + se.fit.day), se = FALSE, color = "black", linetype = 2) +
  geom_smooth(mapping = aes(x = Minute, y = fit.day), color = "black", se = FALSE, size = 1.1) +
  geom_smooth(mapping = aes(x = Minute, y = fit.day - se.fit.day), se = FALSE, color = "black", linetype = 2) +
  facet_grid(~ Day) +
  theme_bw(base_size = 16) +
  labs(x = "Minute", y = "Total Vocal Prevalence") +
  scale_color_brewer(palette = "Dark2") +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(hjust = 0.5),
        axis.text.y = element_text(hjust = 1),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold"),
        strip.text.y = element_text(face = "bold"),
        panel.spacing = unit(1, "lines"),
        aspect.ratio = 0.8)

ggview(device = "jpeg", units = "in", dpi = 1200, width = 12, height = 4)

ggsave("Figures/Fig3B-TVP.jpg", dpi = 1200, width = 12, height = 4)







# SPECIES VOCAL PREVALENCE -----------------------------------------------------

VP.Window.60 <- read.csv("Outputs/VP.Window.60") %>% 
  mutate(Minute = Time.Window/60,
         Site = as.factor(Site),
         Day = as.factor(Day),
         Hab1 = as.factor(Hab1),
         SiteDay = as.factor(str_c(Site, Day))) 

# OUTLIERS & NORMALITY OF RESPONSE VARIABLES -----------------------------------
ggdensity(VP.Window.60$VP, xlab = "Vocal Prevalence")
gghistogram(VP.Window.60$VP, xlab = "Vocal Prevalence")
ggqqplot(VP.Window.60$VP, xlab = "Vocal Prevalence")
# going to need to model this as a binomial distribution

Most.present.spp <- VP.Window.60 %>% 
  group_by(Species) %>% 
  summarize(sum = sum(VP)) %>% 
  view()
# Hauxwell's Thrush, Black-faced Antthrush, Thrush-like Wren

# SPECIES VOCAL PREVALENCE MODEL -----------------------------------------------

# Hauxwell's Thrush 
VP.Window.60.HATH <- VP.Window.60 %>% 
  filter(Species == "HATH")

gghistogram(VP.Window.60.HATH$VP, xlab = "VP")
ggqqplot(VP.Window.60.HATH$VP, ylab = "VP")

VPgam.60.ar1.HATH <- mgcv::gamm(VP ~ s(Minute) + s(Day, bs = "re") + s(Site, bs = "re"),
                            correlation = corAR1(form = ~ Minute | SiteDay),
                            family = "poisson",
                            method = "REML",
                            data = VP.Window.60.HATH)

# checking for autocorrelation issues
performance::check_singularity(VPgam.60.ar1.HATH$gam)
acf(resid(VPgam.60.ar1.HATH$lme, type = "normalized")) # visually, we do not have autocorrelated residuals
pacf(resid(VPgam.60.ar1.HATH$lme, type = "normalized"))

# checking model diagnostics
par(mfrow=c(2,2))
gam.check(VPgam.60.ar1.HATH$gam) # k values are too small, but we can't change them, residuals look great
concurvity(VPgam.60.ar1.HATH$gam, full = TRUE) # no issues with concurvity
concurvity(VPgam.60.ar1.HATH$gam, full = FALSE) # no issues with concurvity

summary(VPgam.60.ar1.HATH$gam)
anova.gam(VPgam.60.ar1.HATH$gam)
# visualizing partial effects
plot(VPgam.60.ar1.HATH$gam, shade = TRUE, shift = coef(VPgam.60.ar1.HATH$gam)[1],
     trans = exp, pages = 1, all.terms = TRUE, rug = FALSE)


# Black-faced Antthrush
VP.Window.60.BFAT <- VP.Window.60 %>% 
  filter(Species == "BFAT")

gghistogram(VP.Window.60.BFAT$VP, xlab = "VP")
ggqqplot(VP.Window.60.BFAT$VP, ylab = "VP")

VPgam.60.ar1.BFAT <- mgcv::gamm(VP ~ s(Minute) + s(Day, bs = "re") + s(Site, bs = "re"),
                                correlation = corAR1(form = ~ Minute | SiteDay),
                                family = "poisson",
                                method = "REML",
                                data = VP.Window.60.BFAT)

# checking for autocorrelation issues
performance::check_singularity(VPgam.60.ar1.BFAT$gam)
acf(resid(VPgam.60.ar1.BFAT$lme, type = "normalized")) # visually, we do not have autocorrelated residuals
pacf(resid(VPgam.60.ar1.BFAT$lme, type = "normalized"))

# checking model diagnostics
par(mfrow=c(2,2))
gam.check(VPgam.60.ar1.BFAT$gam) # k values are too small, but we can't change them, residuals look great
concurvity(VPgam.60.ar1.BFAT$gam, full = TRUE) # no issues with concurvity
concurvity(VPgam.60.ar1.BFAT$gam, full = FALSE) # no issues with concurvity

summary(VPgam.60.ar1.BFAT$gam)
anova.gam(VPgam.60.ar1.BFAT$gam)
# visualizing partial effects
plot(VPgam.60.ar1.BFAT$gam, shade = TRUE, shift = coef(VPgam.60.ar1.BFAT$gam)[1],
     trans = exp, pages = 1, all.terms = TRUE, rug = FALSE)


#Thrush-like Wren
VP.Window.60.TLWR <- VP.Window.60 %>% 
  filter(Species == "TLWR")

gghistogram(VP.Window.60.TLWR$VP, xlab = "VP")
ggqqplot(VP.Window.60.TLWR$VP, ylab = "VP")

VPgam.60.ar1.TLWR <- mgcv::gamm(VP ~ s(Minute) + s(Day, bs = "re") + s(Site, bs = "re"),
                                correlation = corAR1(form = ~ Minute | SiteDay),
                                family = "poisson",
                                method = "REML",
                                data = VP.Window.60.TLWR)

# checking for autocorrelation issues
performance::check_singularity(VPgam.60.ar1.TLWR$gam)
acf(resid(VPgam.60.ar1.TLWR$lme, type = "normalized")) # visually, we do not have autocorrelated residuals
pacf(resid(VPgam.60.ar1.TLWR$lme, type = "normalized"))

# checking model diagnostics
par(mfrow=c(2,2))
gam.check(VPgam.60.ar1.TLWR$gam) # k values are too small, but we can't change them, residuals look great
concurvity(VPgam.60.ar1.TLWR$gam, full = TRUE) # no issues with concurvity
concurvity(VPgam.60.ar1.TLWR$gam, full = FALSE) # no issues with concurvity

summary(VPgam.60.ar1.TLWR$gam)
anova.gam(VPgam.60.ar1.TLWR$gam)
# visualizing partial effects
plot(VPgam.60.ar1.TLWR$gam, shade = TRUE, shift = coef(VPgam.60.ar1.TLWR$gam)[1],
     trans = exp, pages = 1, all.terms = TRUE, rug = FALSE)


# Piratic Flycatcher
VP.Window.60.PIFL <- VP.Window.60 %>% 
  filter(Species == "PIFL")

gghistogram(VP.Window.60.PIFL$VP, xlab = "VP")
ggqqplot(VP.Window.60.PIFL$VP, ylab = "VP")

VPgam.60.ar1.PIFL <- mgcv::gamm(VP ~ s(Minute) + s(Day, bs = "re") + s(Site, bs = "re"),
                                correlation = corAR1(form = ~ Minute | SiteDay),
                                family = "poisson",
                                method = "REML",
                                data = VP.Window.60.PIFL)

# checking for autocorrelation issues
performance::check_singularity(VPgam.60.ar1.PIFL$gam)
acf(resid(VPgam.60.ar1.PIFL$lme, type = "normalized")) # visually, we do not have autocorrelated residuals
pacf(resid(VPgam.60.ar1.PIFL$lme, type = "normalized"))

# checking model diagnostics
par(mfrow=c(2,2))
gam.check(VPgam.60.ar1.PIFL$gam) # k values are too small, but we can't change them, residuals look great
concurvity(VPgam.60.ar1.PIFL$gam, full = TRUE) # no issues with concurvity
concurvity(VPgam.60.ar1.PIFL$gam, full = FALSE) # no issues with concurvity

summary(VPgam.60.ar1.PIFL$gam)
anova.gam(VPgam.60.ar1.PIFL$gam)
# visualizing partial effects
plot(VPgam.60.ar1.PIFL$gam, shade = TRUE, shift = coef(VPgam.60.ar1.PIFL$gam)[1],
     trans = exp, pages = 1, all.terms = TRUE, rug = FALSE)


# Buff-throated Woodcreeper
VP.Window.60.BTWO <- VP.Window.60 %>% 
  filter(Species == "BTWO")

gghistogram(VP.Window.60.BTWO$VP, xlab = "VP")
ggqqplot(VP.Window.60.BTWO$VP, ylab = "VP")

VPgam.60.ar1.BTWO <- mgcv::gamm(VP ~ s(Minute) + s(Day, bs = "re") + s(Site, bs = "re"),
                                correlation = corAR1(form = ~ Minute | SiteDay),
                                family = "poisson",
                                method = "REML",
                                data = VP.Window.60.BTWO)

# checking for autocorrelation issues
performance::check_singularity(VPgam.60.ar1.BTWO$gam)
acf(resid(VPgam.60.ar1.BTWO$lme, type = "normalized")) # visually, we do not have autocorrelated residuals
pacf(resid(VPgam.60.ar1.BTWO$lme, type = "normalized"))

# checking model diagnostics
par(mfrow=c(2,2))
gam.check(VPgam.60.ar1.BTWO$gam) # k values are too small, but we can't change them, residuals look great
concurvity(VPgam.60.ar1.BTWO$gam, full = TRUE) # no issues with concurvity
concurvity(VPgam.60.ar1.BTWO$gam, full = FALSE) # no issues with concurvity

summary(VPgam.60.ar1.BTWO$gam)
anova.gam(VPgam.60.ar1.BTWO$gam)
# visualizing partial effects
plot(VPgam.60.ar1.BTWO$gam, shade = TRUE, shift = coef(VPgam.60.ar1.BTWO$gam)[1],
     trans = exp, pages = 1, all.terms = TRUE, rug = FALSE)


# Amazonian Motmot
VP.Window.60.AMMO <- VP.Window.60 %>% 
  filter(Species == "AMMO")

gghistogram(VP.Window.60.AMMO$VP, xlab = "VP")
ggqqplot(VP.Window.60.AMMO$VP, ylab = "VP")

VPgam.60.ar1.AMMO <- mgcv::gamm(VP ~ s(Minute) + s(Day, bs = "re") + s(Site, bs = "re"),
                                correlation = corAR1(form = ~ Minute | SiteDay),
                                family = "poisson",
                                method = "REML",
                                data = VP.Window.60.AMMO)

# checking for autocorrelation issues
performance::check_singularity(VPgam.60.ar1.AMMO$gam)
acf(resid(VPgam.60.ar1.AMMO$lme, type = "normalized")) # visually, we do not have autocorrelated residuals
pacf(resid(VPgam.60.ar1.AMMO$lme, type = "normalized"))

# checking model diagnostics
par(mfrow=c(2,2))
gam.check(VPgam.60.ar1.AMMO$gam) # k values are too small, but we can't change them, residuals look great
concurvity(VPgam.60.ar1.AMMO$gam, full = TRUE) # no issues with concurvity
concurvity(VPgam.60.ar1.AMMO$gam, full = FALSE) # no issues with concurvity

summary(VPgam.60.ar1.AMMO$gam)
anova.gam(VPgam.60.ar1.AMMO$gam)
# visualizing partial effects
plot(VPgam.60.ar1.AMMO$gam, shade = TRUE, shift = coef(VPgam.60.ar1.AMMO$gam)[1],
     trans = exp, pages = 1, all.terms = TRUE, rug = FALSE)


# Buff-breasted Wren
VP.Window.60.BBWR <- VP.Window.60 %>% 
  filter(Species == "BBWR")

gghistogram(VP.Window.60.BBWR$VP, xlab = "VP")
ggqqplot(VP.Window.60.BBWR$VP, ylab = "VP")

VPgam.60.ar1.BBWR <- mgcv::gamm(VP ~ s(Minute) + s(Day, bs = "re") + s(Site, bs = "re"),
                                correlation = corAR1(form = ~ Minute | SiteDay),
                                family = "poisson",
                                method = "REML",
                                data = VP.Window.60.BBWR)

# checking for autocorrelation issues
performance::check_singularity(VPgam.60.ar1.BBWR$gam)
acf(resid(VPgam.60.ar1.BBWR$lme, type = "normalized")) # visually, we do not have autocorrelated residuals
pacf(resid(VPgam.60.ar1.BBWR$lme, type = "normalized"))

# checking model diagnostics
par(mfrow=c(2,2))
gam.check(VPgam.60.ar1.BBWR$gam) # k values are too small, but we can't change them, residuals look great
concurvity(VPgam.60.ar1.BBWR$gam, full = TRUE) # no issues with concurvity
concurvity(VPgam.60.ar1.BBWR$gam, full = FALSE) # no issues with concurvity

summary(VPgam.60.ar1.BBWR$gam)
anova.gam(VPgam.60.ar1.BBWR$gam)
# visualizing partial effects
plot(VPgam.60.ar1.BBWR$gam, shade = TRUE, shift = coef(VPgam.60.ar1.BBWR$gam)[1],
     trans = exp, pages = 1, all.terms = TRUE, rug = FALSE)


# Amazonian Barred Woodcreeper
VP.Window.60.AMBW <- VP.Window.60 %>% 
  filter(Species == "AMBW")

gghistogram(VP.Window.60.AMBW$VP, xlab = "VP")
ggqqplot(VP.Window.60.AMBW$VP, ylab = "VP")

VPgam.60.ar1.AMBW <- mgcv::gamm(VP ~ s(Minute) + s(Day, bs = "re") + s(Site, bs = "re"),
                                correlation = corAR1(form = ~ Minute | SiteDay),
                                family = "poisson",
                                method = "REML",
                                data = VP.Window.60.AMBW)

# checking for autocorrelation issues
performance::check_singularity(VPgam.60.ar1.AMBW$gam)
acf(resid(VPgam.60.ar1.AMBW$lme, type = "normalized")) # visually, we do not have autocorrelated residuals
pacf(resid(VPgam.60.ar1.AMBW$lme, type = "normalized"))

# checking model diagnostics
par(mfrow=c(2,2))
gam.check(VPgam.60.ar1.AMBW$gam) # k values are too small, but we can't change them, residuals look great
concurvity(VPgam.60.ar1.AMBW$gam, full = TRUE) # no issues with concurvity
concurvity(VPgam.60.ar1.AMBW$gam, full = FALSE) # no issues with concurvity

summary(VPgam.60.ar1.AMBW$gam)
anova.gam(VPgam.60.ar1.AMBW$gam)
# visualizing partial effects
plot(VPgam.60.ar1.AMBW$gam, shade = TRUE, shift = coef(VPgam.60.ar1.AMBW$gam)[1],
     trans = exp, pages = 1, all.terms = TRUE, rug = FALSE)



# MODEL PREDICTION --------------------------------------------------------

# predicting from the model
predicted.VP.HATH <- data.frame(predict(VPgam.60.ar1.HATH$gam, type = "response", se.fit = TRUE))
predicted.VP.HATH <- cbind(VP.Window.60.HATH, predicted.VP.HATH)

predicted.VP.BFAT <- data.frame(predict(VPgam.60.ar1.BFAT$gam, type = "response", se.fit = TRUE))
predicted.VP.BFAT <- cbind(VP.Window.60.BFAT, predicted.VP.BFAT)

#predicted.VP.TLWR <- data.frame(predict(VPgam.60.ar1.TLWR$gam, type = "response", se.fit = TRUE))
#predicted.VP.TLWR <- cbind(VP.Window.60.TLWR, predicted.VP.TLWR)

predicted.VP.PIFL <- data.frame(predict(VPgam.60.ar1.PIFL$gam, type = "response", se.fit = TRUE))
predicted.VP.PIFL <- cbind(VP.Window.60.PIFL, predicted.VP.PIFL)

predicted.VP.AMMO <- data.frame(predict(VPgam.60.ar1.AMMO$gam, type = "response", se.fit = TRUE))
predicted.VP.AMMO <- cbind(VP.Window.60.AMMO, predicted.VP.AMMO)

predicted.VP.BTWO <- data.frame(predict(VPgam.60.ar1.BTWO$gam, type = "response", se.fit = TRUE))
predicted.VP.BTWO <- cbind(VP.Window.60.BTWO, predicted.VP.BTWO)

#predicted.VP.BBWR <- data.frame(predict(VPgam.60.ar1.BBWR$gam, type = "response", se.fit = TRUE))
#predicted.VP.BBWR <- cbind(VP.Window.60.BBWR, predicted.VP.BBWR)

predicted.VP.AMBW <- data.frame(predict(VPgam.60.ar1.AMBW$gam, type = "response", se.fit = TRUE))
predicted.VP.AMBW <- cbind(VP.Window.60.AMBW, predicted.VP.AMBW)

predicted.VP.spp <- rbind(predicted.VP.HATH, predicted.VP.BFAT, predicted.VP.PIFL,
                          predicted.VP.BTWO, predicted.VP.AMMO, predicted.VP.AMBW) %>% 
  # converting species codes to names
  mutate(Species = if_else(Species == "AMMO", "Amazonian Motmot", Species),
         Species = if_else(Species == "AMBW", "Amazonian Barred Woodcreeper", Species),
         Species = if_else(Species == "BFAT", "Black-faced Antthrush", Species),
         Species = if_else(Species == "BTWO", "Buff-throated Woodcreeper", Species),
         Species = if_else(Species == "HATH", "Hauxwell's Thrush", Species),
         Species = if_else(Species == "PIFL", "Piratic Flycatcher", Species))

ggplot(data = predicted.VP.spp) +
  #geom_point(mapping = aes(x = Minute, y = VP)) +
  #geom_line(mapping = aes(x = Minute, y = fit)) +
  geom_smooth(mapping = aes(x = Minute, y = fit + se.fit), se = FALSE, color = "black", linetype = 2) +
  geom_smooth(mapping = aes(x = Minute, y = fit), se = FALSE, color = "black", size = 1.1) +
  geom_smooth(mapping = aes(x = Minute, y = fit - se.fit), se = FALSE, color = "black", linetype = 2) +
  facet_wrap(~ Species, ncol = 3) +
  theme_bw(base_size = 16) +
  labs(x = "Minute", y = "Vocal Prevalence") +
  scale_color_brewer(palette = "Dark2") +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(hjust = 0.5),
        axis.text.y = element_text(hjust = 1),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold"),
        strip.text.y = element_text(face = "bold"),
        panel.spacing = unit(1, "lines"),
        aspect.ratio = 1)

ggview(device = "jpeg", units = "in", dpi = 1200, width = 10, height = 7.5)

ggsave("Figures/Fig5-VP.png", dpi = 1200, width = 10, height = 7.5)

