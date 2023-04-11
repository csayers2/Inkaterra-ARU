
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

# smoothing predictions by day
ggplot(data = SR.Window.60) +
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

# predictions across days and sites
ggplot(data = SR.Window.60) +
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

ggsave("Figures/Fig1-SR.jpg", dpi = 1200, width = 6, height = 12)


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

ggsave("Figures/Fig2-TVP.jpg", dpi = 1200, width = 6, height = 12)


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

VPgam.60.ar1.HATH <- mgcv::gamm(VP ~ s(Minute, by = SiteDay) + s(Day, bs = "re") + s(Site, bs = "re"),
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

VPgam.60.ar1.BFAT <- mgcv::gamm(cbind(VP, VA) ~ s(Minute, by = Day) + s(Day, bs = "re") + s(Site, bs = "re"),
                                correlation = corAR1(form = ~ Minute | SiteDay),
                                family = binomial(link = "logit"),
                                method = "REML",
                                data = VP.Window.60.BFAT)

# checking for autocorrelation issues
performance::check_singularity(VP.Window.60.BFAT$gam)
acf(resid(VP.Window.60.BFAT$lme, type = "normalized")) # visually, we do not have autocorrelated residuals
pacf(resid(VP.Window.60.BFAT$lme, type = "normalized"))

# checking model diagnostics
par(mfrow=c(2,2))
gam.check(VP.Window.60.BFAT$gam) # k values are too small, but we can't change them, residuals look great
concurvity(VP.Window.60.BFAT$gam, full = TRUE) # no issues with concurvity
concurvity(VP.Window.60.BFAT$gam, full = FALSE) # no issues with concurvity

summary(VP.Window.60.BFAT$gam)
anova.gam(VP.Window.60.BFAT$gam)
# visualizing partial effects
plot(VP.Window.60.BFAT$gam, shade = TRUE, shift = coef(VP.Window.60.BFAT$gam)[1],
     trans = exp, pages = 1, all.terms = TRUE, rug = FALSE)



# Thrush-like Wren 
VP.Window.60.TLWR <- VP.Window.60 %>% 
  filter(Species == "TLWR")

gghistogram(VP.Window.60.TLWR$VP, xlab = "VP")
ggqqplot(VP.Window.60.TLWR$VP, ylab = "VP")

VPgam.60.ar1.TLWR <- mgcv::gamm(cbind(VP, VA) ~ s(Minute) + s(Day, bs = "re") + s(Site, bs = "re"),
                                correlation = corAR1(form = ~ Minute | SiteDay),
                                family = binomial(link = "logit"),
                                method = "REML",
                                data = VP.Window.60.TLWR)

# checking for autocorrelation issues
performance::check_singularity(VP.Window.60.TLWR$gam)
acf(resid(VP.Window.60.TLWR$lme, type = "normalized")) # visually, we do not have autocorrelated residuals
pacf(resid(VP.Window.60.TLWR$lme, type = "normalized"))

# checking model diagnostics
par(mfrow=c(2,2))
gam.check(VP.Window.60.TLWR$gam) # k values are too small, but we can't change them, residuals look great
concurvity(VP.Window.60.TLWR$gam, full = TRUE) # no issues with concurvity
concurvity(VP.Window.60.TLWR$gam, full = FALSE) # no issues with concurvity

summary(VP.Window.60.TLWR$gam)
anova.gam(VP.Window.60.TLWR$gam)
# visualizing partial effects
plot(VP.Window.60.TLWR$gam, shade = TRUE, shift = coef(VP.Window.60.TLWR$gam)[1],
     trans = exp, pages = 1, all.terms = TRUE, rug = FALSE)


