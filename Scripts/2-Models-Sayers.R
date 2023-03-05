
# Chris Sayers
# updated March 4, 2023

# script designed for data visualizations and model fitting

library(ggpubr)
library(MuMIn)
library(glmmTMB)
library(performance)
library(car)
library(DHARMa)
library(lawstat)
library(emmeans)
library(multcomp)

# SPECIES RICHNESS --------------------------------------------------------

SR.Total <- read.csv("Outputs/SR.Total")
SR.Window.60 <- read.csv("Outputs/SR.Window.60")

# OUTLIERS & NORMALITY OF RESPONSE VARIABLES -----------------------------------
ggdensity(SR.Window.60$SR, xlab = "Species Richness")
gghistogram(SR.Window.60$SR, xlab = "Species Richness")
ggqqplot(SR.Window.60$SR, ylab = "Species Richness")
shapiro.test(SR.Window.60$SR) # W = 0.97624, p-value = 2.612e-12, not normal

library(ggplot2)
ggplot(data = SR.Window.60) +
  geom_point(mapping = aes(x = Time.Window, y = SR, color = Day)) +
  geom_smooth(mapping = aes(x = Time.Window, y = SR, color = Day)) 
  #facet_grid(Day ~ Hab2)

# SPECIES RICHNESS MODEL -------------------------------------------------------

SRmodel.60.1 <- glmmTMB(SR ~ poly(Time.Window, 2)*Day + Hab2 + (1 | Site),
                        data = SR.Window.60, family = "poisson", REML = F)

SRmodel.60.2 <- glmmTMB(SR ~ poly(Time.Window, 2)*Day*Site,
                      data = SR.Window.60, family = "poisson", REML = F)

AIC(SRmodel.60.1, SRmodel.60.2)

performance::r2(SRmodel.60)
car::Anova(SRmodel.60, type = 3)

summary(SRmodel.60)
as.data.frame(confint(SRmodel.60)) %>% 
  mutate(Estimate = exp(Estimate), `2.5 %` = exp(`2.5 %`), `97.5 %` = exp(`97.5 %`))
performance::r2(SRmodel.60)
car::Anova(SRmodel.60, type = 3)

# CHECKING MODEL ASSUMPTIONS -------------------------------------
# Checking for homogeneity of variance & normality of residuals
mean(residuals(SRmodel.60)) # VERY close to 0

simulateResiduals(SRmodel.60, plot = T, refit = F, use.u = T)
shapiro.test(residuals(SRmodel.60)) # W = 0.99486, p-value = 0.0009878, not normal
# residual plots look okay

## Checking for autocorrelation/independence
#acf(SR.Window.60$SR) # raw data is autocorrelated
#acf(residuals(SRmodel.60)) # random effects variable corrects for this
#runs.test(residuals(SRmodel.60)) # we do not have autocorrelated data
#
#library(DataCombine)
#SR.Window.lag <- data.frame(SR.Window.60, resid.mod = residuals(SRmodel.60)) %>% 
#  slide(Var = "resid.mod", NewVar = "lag1", slideBy = -1) %>% 
#  na.omit()
#
#SRmodel.lag <- glmmTMB(SR ~ poly(Time.Window, 2)*Day + lag1 + (1 + Day | Site),
#                   data = SR.Window.lag, family = "gaussian", REML = F)
#
#acf(SR.Window.lag$SR) # raw data is autocorrelated
#acf(residuals(SRmodel.lag)) # random effects variable corrects for this
#runs.test(residuals(SRmodel.lag)) # we do not have autocorrelated data


# MODEL SELECTION ---------------------------------------------------------
SRmodel.60 <- glmmTMB(SR ~ poly(Time.Window, 2)*Day + Hab2 + (1 | Site),
                      data = SR.Window.60, family = "poisson", REML = F)

options(na.action = "na.fail")
# computes marginal and conditional R^2
d.out <- MuMIn::dredge(SRmodel.60, extra = list("Rsq" = function(x){performance::r2(x)}))
View(d.out)
options(na.action = "na.omit")
write.csv(d.out, "Outputs/sr-model-selection.csv")

# 1st place model by a long-shot (R2 = 0.64, wi = 76%)
topSRmodel <- glmmTMB(SR ~ poly(Time.Window, 2)*Day + Hab2 + (1 | Site),
                      data = SR.Window.60, family = "poisson", REML = F)

summary(topSRmodel)
as.data.frame(confint(topSRmodel)) %>% 
  mutate(Estimate = exp(Estimate), `2.5 %` = exp(`2.5 %`), `97.5 %` = exp(`97.5 %`))
performance::r2(topSRmodel)
car::Anova(topSRmodel, type = 3)

# computing post-hoc comparisons to determine significant differences among the modeled means

emmeans(topSRmodel, "Day", type = "response") %>% 
  cld(Letter = "abcdefg")

emmeans(topSRmodel, "Hab2", type = "response") %>% 
  cld(Letter = "abcdefg")

# a cool visualization of the predicted means and variances
library(ggeffects)
plot(ggpredict(SRmodel.60, terms = c("Time.Window [all]", "Day", "Hab2"),
               type = "random", plot = T))






# TOTAL VOCAL PREVALENCE --------------------------------------------------------

TVP.Window.60 <- read.csv("Outputs/TVP.Window.60")

# OUTLIERS & NORMALITY OF RESPONSE VARIABLES -----------------------------------
# This distribution is to be expected when dealing with count data/ many zero counts
ggdensity(TVP.Window.60$TVP, xlab = "Total Vocal Prevalence")
gghistogram(TVP.Window.60$TVP, xlab = "Total Vocal Prevalence")
ggqqplot(TVP.Window.60$TVP, ylab = "Total Vocal Prevalence")
shapiro.test(TVP.Window.60$TVP) # W = 0.98492, p-value = 3.958e-09, not normal

library(ggplot2)
ggplot(data = TVP.Window.60) +
  geom_point(mapping = aes(x = Time.Window, y = TVP, color = Site)) +
  geom_smooth(mapping = aes(x = Time.Window, y = TVP)) + 
  facet_grid(~ Day)

# TOTAL VOCAL PREVALENCE MODEL -------------------------------------------------------

TVPmodel.60 <- glmmTMB(TVP ~ poly(Time.Window, 2)*Day + Hab2 + (1 | Site),
                      data = TVP.Window.60, family = "poisson", REML = F)

TVPmodel.60 <- glmmTMB(TVP ~ poly(Time.Window, 2)*Day*Site,
                       data = TVP.Window.60, family = "poisson", REML = F)

TVPmodel.60 <- glmmTMB(TVP ~ poly(Time.Window, 2)*Day + Site,
                       data = TVP.Window.60, family = "poisson", REML = F)

summary(TVPmodel.60)
as.data.frame(confint(TVPmodel.60)) %>% 
  mutate(Estimate = exp(Estimate), `2.5 %` = exp(`2.5 %`), `97.5 %` = exp(`97.5 %`))
performance::r2(TVPmodel.60)
car::Anova(TVPmodel.60, type = 3)


# CHECKING MODEL ASSUMPTIONS -------------------------------------
# Checking for homogeneity of variance & normality of residuals
mean(residuals(TVPmodel.60)) # VERY close to 0

simulateResiduals(TVPmodel.60, plot = T, refit = F, use.u = T)
runs.test(residuals(TVPmodel.60)) # W = 0.99464, p-value = 0.6516, normal!
# residual plots look okay

# Checking for autocorrelation/independence
acf(TVP.Window.60$TVP) # raw data is autocorrelated
acf(residuals(TVPmodel.60)) # random effects variable corrects for this
runs.test(residuals(TVPmodel.60)) # we have autocorrelated data

library(DataCombine)
TVP.Window.60.lag <- data.frame(TVP.Window.60, resid.mod = residuals(TVPmodel.60)) %>% 
  slide(Var = "resid.mod", NewVar = "lag1", slideBy = -1) %>% 
  na.omit()

TVPmodel.60.lag <- glmmTMB(TVP ~ poly(Time.Window, 2)*Day + Hab2 + lag1 + (1 | Site),
                   data = TVP.Window.60.lag, family = "poisson", REML = F)

acf(TVP.Window.60.lag$TVP) # raw data is autocorrelated
acf(residuals(TVPmodel.60.lag)) # random effects variable corrects for this
runs.test(residuals(TVPmodel.60.lag)) # we still have autocorrelated data

# Checking for homogeneity of variance & normality of residuals
mean(residuals(TVPmodel.60.lag)) # VERY close to 0

simulateResiduals(TVPmodel.60.lag, plot = T, refit = F, use.u = T)
runs.test(residuals(TVPmodel.60.lag)) # W = 0.99464, p-value = 0.6516, normal!
# residual plots look okay


# MODEL SELECTION ---------------------------------------------------------
TVPmodel.60 <- glmmTMB(TVP ~ poly(Time.Window, 2)*Day + Hab2 + (1 | Site),
                       data = TVP.Window.60, family = "poisson", REML = F)

options(na.action = "na.fail")
# computes marginal and conditional R^2
d.out <- MuMIn::dredge(TVPmodel.60, extra = list("Rsq" = function(x){performance::r2(x)}))
View(d.out)
options(na.action = "na.omit")
write.csv(d.out, "Outputs/tvp-model-selection.csv")

# 1st place model by a long-shot (R2 = 0.87, wi = 73%)
topTVPmodel.60 <- glmmTMB(TVP ~ poly(Time.Window, 2)*Day + (1 | Site),
                      data = TVP.Window.60, family = "poisson", REML = F)

summary(topTVPmodel.60)
as.data.frame(confint(topTVPmodel.60)) %>% 
  mutate(Estimate = exp(Estimate), `2.5 %` = exp(`2.5 %`), `97.5 %` = exp(`97.5 %`))
performance::r2(topTVPmodel.60)
car::Anova(topTVPmodel.60, type = 3)

# computing post-hoc comparisons to determine significant differences among the modeled means
# need to switch model to be as.factor(Time.Window) first before performing this comparison
emmeans(topTVPmodel.60, "Day", type = "response") %>% 
  cld(Letter = "abcdefg")

#emmeans(topTVPmodel.60, "Site", type = "response") %>% 
#  cld(Letter = "abcdefg")

# a cool visualization of the predicted means and variances
library(ggeffects)
plot(ggpredict(topTVPmodel.60, terms = c("Time.Window [all]", "Day"),
               type = "random", plot = T))






















# SPECIES-SPECIFIC VOCAL PREVALENCE --------------------------------------------------------

VP.Window.10 <- read.csv("Outputs/VP.Window.10")
VP.Window.60 <- read.csv("Outputs/VP.Window.60")
#TE.Window.10.spp <- read.csv("Outputs/TE.Window.10.spp")

# OUTLIERS & NORMALITY OF RESPONSE VARIABLES -----------------------------------
# This distribution is to be expected when dealing with count data/ many zero counts

ggdensity(VP.Window.10$VP, xlab = "Vocal Presence")
gghistogram(VP.Window.10$VP, xlab = "Vocal Presence")
ggqqplot(VP.Window.10$VP, ylab = "Vocal Presence")
shapiro.test(VP.Window.10$VP) # W = 0.97989, p-value = 4.43e-11, not normal

ggdensity(VP.Window.60$VP, xlab = "Vocal Prevalence")
gghistogram(VP.Window.60$VP, xlab = "Vocal Prevalence")
ggqqplot(VP.Window.60$VP, ylab = "Vocal Prevalence")
shapiro.test(VP.Window.60$VP) # W = 0.97989, p-value = 4.43e-11, not normal

library(ggplot2)
VP.Window.60 %>% 
  filter(VP != 0) %>% 
  ggplot() +
  geom_point(mapping = aes(x = Time.Window, y = VP, color = Site)) +
  geom_smooth(mapping = aes(x = Time.Window, y = VP, color = Site)) + 
  facet_grid(Day ~ Hab2)

# visualizing total encounters
ggplot(data = TE.Window.10.spp) +
  geom_boxplot(mapping = aes(x = Site:Day, y = TE)) +
  #geom_smooth(mapping = aes(x = Time.Window, y = VP, color = Site)) + 
  facet_grid(~ Hab2)

# SPECIES-SPECIFIC VOCAL PREVALENCE MODEL -------------------------------------------------------

#VPmodel.10 <- glmmTMB(VP ~ Time.Window + Day + Hab2 + Edge.Distance +
#                        (1 | Site) + (1 | Species),
#                      data = VP.Window.10, family = "binomial", REML = F)
#
#summary(VPmodel.10)
#as.data.frame(confint(VPmodel.10)) %>% 
#  mutate(Estimate = exp(Estimate), `2.5 %` = exp(`2.5 %`), `97.5 %` = exp(`97.5 %`))
#performance::r2(VPmodel.10)
#car::Anova(VPmodel.10, type = 3)

#VPmodel.60.1 <- glmmTMB(cbind(VP, VA) ~ Time.Window + Day + Site + (1 | Species),
#                        data = VP.Window.60, family = "binomial", REML = F)
#
#VPmodel.60.2 <- glmmTMB(cbind(VP, VA) ~ Time.Window*Day + Site + (1 | Species),
#                        data = VP.Window.60, family = "binomial", REML = F)
#
#VPmodel.60.3 <- glmmTMB(cbind(VP, VA) ~ poly(Time.Window, 2) + Day + Site + (1 | Species),
#                        data = VP.Window.60, family = "binomial", REML = F)
#
#VPmodel.60.4 <- glmmTMB(cbind(VP, VA) ~ poly(Time.Window, 2)*Day + Site + (1 | Species),
#                        data = VP.Window.60, family = "binomial", REML = F)
#
#VPmodel.60.5 <- glmmTMB(cbind(VP, VA) ~ Time.Window + Day + Hab2 + (1 | Site) + (1 | Species),
#                        data = VP.Window.60, family = "binomial", REML = F)
#
#VPmodel.60.6 <- glmmTMB(cbind(VP, VA) ~ poly(Time.Window, 2) + Day + Hab2 + (1 | Site) + (1 | Species),
#                        data = VP.Window.60, family = "binomial", REML = F)





VPmodel.60 <- glmmTMB(cbind(VP, VA) ~ poly(Time.Window, 2)*Day + Hab2 + (1 | Site) + (1 | Species),
                      data = VP.Window.60, family = "binomial", REML = F)

VPmodel.60.1 <- glmmTMB(cbind(VP, VA) ~ poly(Time.Window, 2)*Day + Site + (1 | Species),
                      data = VP.Window.60, family = "binomial", REML = F)

# won't converge
#VPmodel.60.2 <- glmmTMB(cbind(VP, VA) ~ poly(Time.Window, 2)*Day + Hab2 + (1 + Day | Site) + (1 + Time.Window | Species),
#                        data = VP.Window.60, family = "binomial", REML = F)

VPmodel.60.3 <- glmmTMB(cbind(VP, VA) ~ poly(Time.Window, 2)*Day + Hab2 + (1 | Site) + (1 + Time.Window | Species),
                      data = VP.Window.60, family = "binomial", REML = F)

VPmodel.60.4 <- glmmTMB(cbind(VP, VA) ~ poly(Time.Window, 2) + Day + Hab2 + (1 | Site) + (1 + Time.Window | Species),
                        data = VP.Window.60, family = "binomial", REML = F)

AIC(VPmodel.60.1, VPmodel.60.2, VPmodel.60.3, VPmodel.60)

summary(VPmodel.60)
as.data.frame(confint(VPmodel.60)) %>% 
  mutate(Estimate = exp(Estimate), `2.5 %` = exp(`2.5 %`), `97.5 %` = exp(`97.5 %`))
performance::r2(VPmodel.60)
car::Anova(VPmodel.60, type = 3)


#VPmodel.60 <- glmmTMB(cbind(VP, VA) ~ Time.Window + Day + Hab2 +
#                        (1 | Site) + (1 | Species),
#                        #(1 + Time.Window*Day | Species)
#                      data = VP.Window.60, family = "binomial", REML = F)
#
#VPmodel.60.2 <- glmmTMB(cbind(VP, VA) ~ Time.Window + Day + Hab2 +
#                        (1 | Site) + (1 + Time.Window | Species),
#                      #(1 + Time.Window*Day | Species)
#                      data = VP.Window.60, family = "binomial", REML = F)
#
#VPmodel.60.3 <- glmmTMB(cbind(VP, VA) ~ Time.Window + poly(Day, 2) + Hab2 +
#                          (1 | Site) + (1 + Time.Window | Species),
#                        #(1 + Day | Site) + (1 + Time.Window*Day | Species)
#                        data = VP.Window.60, family = "binomial", REML = F)

# CHECKING MODEL ASSUMPTIONS -------------------------------------
# Checking for homogeneity of variance & normality of residuals
mean(residuals(VPmodel.60)) # VERY close to 0

simulateResiduals(VPmodel.60, plot = T, refit = F, use.u = T)
runs.test(residuals(VPmodel.60)) # W = 0.99464, p-value = 0.6516, normal!
# residual plots look okay

## Checking for autocorrelation/independence
#acf(VP.Window.60$VP) # raw data is autocorrelated
#acf(residuals(VPmodel.60)) # random effects variable corrects for this
#runs.test(residuals(VPmodel.60)) # we do not have autocorrelated data
#
#library(DataCombine)
#VP.Window.60.lag <- data.frame(VP.Window.60, resid.mod = residuals(VPmodel.60)) %>% 
#  slide(Var = "resid.mod", NewVar = "lag1", slideBy = -1) %>% 
#  na.omit()
#
#VPmodel.60.lag <- glmmTMB(VP ~ Time.Window*Day*Site + lag1,
#                   data = VP.Window.60.lag, family = "binomial", REML = F)
#
#acf(VP.Window.60.lag$VP) # raw data is autocorrelated
#acf(residuals(VPmodel.60.lag)) # random effects variable corrects for this
#runs.test(residuals(VPmodel.60.lag)) # we do not have autocorrelated data


# MODEL SELECTION ---------------------------------------------------------
VPmodel.60 <- glmmTMB(cbind(VP, VA) ~ poly(Time.Window, 2)*Day + Hab2 + (1 | Site) + (1 | Species),
                      data = VP.Window.60, family = "binomial", REML = F)

options(na.action = "na.fail")
# computes marginal and conditional R^2
d.out <- MuMIn::dredge(VPmodel.60, extra = list("Rsq" = function(x){performance::r2(x)}))
View(d.out)
options(na.action = "na.omit")
write.csv(d.out, "Outputs/vp-model-selection.csv")

# 1st place model by a long-shot (R2 = 0.56, wi = 76%)
topVPmodel <- glmmTMB(cbind(VP, VA) ~ poly(Time.Window, 2)*Day + (1 | Site) + (1 | Species),
                      data = VP.Window.60, family = "binomial", REML = F)

summary(topVPmodel)
as.data.frame(confint(topVPmodel)) %>% 
  mutate(Estimate = exp(Estimate), `2.5 %` = exp(`2.5 %`), `97.5 %` = exp(`97.5 %`))
performance::r2(topVPmodel)
car::Anova(topVPmodel, type = 3)

# computing post-hoc comparisons to determine significant differences among the modeled means
# need to switch model to be as.factor(Time.Window) first before performing this comparison
emmeans(topVPmodel, "Day", type = "response") %>% 
  cld(Letter = "abcdefg")

#emmeans(topVPmodel, "Site", type = "response") %>% 
#  cld(Letter = "abcdefg")

