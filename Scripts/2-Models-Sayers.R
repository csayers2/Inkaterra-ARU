
# Chris Sayers
# updated March 4, 2023

# script designed for data visualizations and model fitting

library(ggplot2)
library(stringr)
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
SR.Window.60 <- read.csv("Outputs/SR.Window.60") %>% 
  mutate(Minute = Time.Window/60,
         fTime.Window = factor(Time.Window),
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
  facet_grid(~ Day)

# SPECIES RICHNESS MODEL -------------------------------------------------------

SRmodel.60 <- glmmTMB(SR ~ (Time.Window^2)*Day + (1 | Site),
                      data = SR.Window.60, family = "poisson", REML = F)

summary(SRmodel.60)
performance::check_singularity(SRmodel.60)

# CHECKING MODEL ASSUMPTIONS -------------------------------------
# Checking for homogeneity of variance & normality of residuals
mean(residuals(SRmodel.60)) # VERY close to 0

simulateResiduals(SRmodel.60, plot = T, refit = F, use.u = T)
shapiro.test(residuals(SRmodel.60)) # W = 0.99803, p-value = 0.2365, normal!
# residual plots look great

# Checking for autocorrelation/independence
acf(SR.Window.60$SR) # raw data is autocorrelated
acf(residuals(SRmodel.60)) # random effects variable does not correct for this
pacf(SR.Window.60$SR)
pacf(residuals(SRmodel.60))
runs.test(residuals(SRmodel.60)) # we have autocorrelated data
lmtest::dwtest(SRmodel.60) # we have autocorrelated data

library(DataCombine)
SR.Window.60.lag <- data.frame(SR.Window.60, resid.mod = residuals(SRmodel.60)) %>% 
  slide(Var = "resid.mod", GroupVar = "SiteDay", NewVar = "lag1", slideBy = -1) %>% 
  na.omit()

SRmodel.60.lag <- glmmTMB(SR ~ (Time.Window^2)*Day + lag1 + (1 | Site),
                          data = SR.Window.60.lag, family = "poisson", REML = F)

acf(SR.Window.60$SR) # raw data is autocorrelated
acf(residuals(SRmodel.60.lag)) # lagged residuals correct for this
pacf(SR.Window.60$SR)
pacf(residuals(SRmodel.60.lag)) # visually, we do not have autocorrelated data
runs.test(residuals(SRmodel.60.lag)) # we do not have autocorrelated data
lmtest::dwtest(SRmodel.60.lag) # we do not have autocorrelated data


performance::check_singularity(SRmodel.60.lag)
summary(SRmodel.60)
as.data.frame(confint(SRmodel.60)) %>% 
  mutate(Estimate = exp(Estimate), `2.5 %` = exp(`2.5 %`), `97.5 %` = exp(`97.5 %`))
performance::r2(SRmodel.60)
car::Anova(SRmodel.60, type = 3)


# MODEL SELECTION ---------------------------------------------------------

# checking the relative importance of site
SRmodel.60.lag <- glmmTMB(SR ~ (Time.Window^2)*Day + lag1 + (1 | Site),
                          data = SR.Window.60.lag, family = "poisson", REML = F)
SRmodel.60.lag.x <- glmmTMB(SR ~ (Time.Window^2)*Day + lag1,
                          data = SR.Window.60.lag, family = "poisson", REML = F)
# model with site as an RE performs a lot better
AIC(SRmodel.60.lag, SRmodel.60.lag.x) 

options(na.action = "na.fail")
# computes marginal and conditional R^2
d.out <- MuMIn::dredge(SRmodel.60.lag, extra = list("Rsq" = function(x){performance::r2(x)}))
View(d.out)
options(na.action = "na.omit")
write.csv(d.out, "Outputs/sr-model-selection.csv")
# model with lagged residuals 

# 1st place model by a long-shot (R2 = 0.64, wi = 76%)
topSRmodel <- glmmTMB(SR ~ poly(Time.Window, 2)*Day + Hab1 + (1 | Site),
                      data = SR.Window.60, family = "poisson", REML = F)

summary(topSRmodel)
as.data.frame(confint(topSRmodel)) %>% 
  mutate(Estimate = exp(Estimate), `2.5 %` = exp(`2.5 %`), `97.5 %` = exp(`97.5 %`))
performance::r2(topSRmodel)
car::Anova(topSRmodel, type = 3)

# computing post-hoc comparisons to determine significant differences among the modeled means

emmeans(SRmodel.60.lag, "Day", type = "response") %>% 
  cld(Letter = "abcdefg")

#emmeans(SRmodel.60.lag, "Hab1", type = "response") %>% 
#  cld(Letter = "abcdefg")

# a cool visualization of the predicted means and variances
library(ggeffects)
plot(ggpredict(SRmodel.60, terms = c("Time.Window [all]", "Day"),
               type = "random", plot = T))









# FIRST ORDER AUTOGRESSIVE STRUCTURE --------------------------------------


# creating a dummy variable for ar1 structure
SR.Window.60 <- SR.Window.60 %>% 
  mutate(group = 1, group = as.factor(group))

# adding a first order autoregressive covariance structure to account for
# temporal autocorrelation
SRmodel.60.ar1 <- glmmTMB(SR ~ (Minute^2)*Day + Site +
                            ar1(factor(Minute) - 1 | SiteDay),
                          data = SR.Window.60, family = "poisson", REML = F)

performance::check_singularity(SRmodel.60.ar1)

acf(SR.Window.60$SR) # raw data is autocorrelated
acf(residuals(SRmodel.60.ar1)) # visually, we do not have autocorrelated data
pacf(SR.Window.60$SR)
pacf(residuals(SRmodel.60.ar1)) # visually, we do not have autocorrelated data
runs.test(residuals(SRmodel.60.ar1)) # we do not have autocorrelated data
lmtest::dwtest(SRmodel.60.ar1) # we do not have autocorrelated data

summary(SRmodel.60.ar1)
performance::r2(SRmodel.60.ar1)
VarCorr(SRmodel.60.ar1, condVar = TRUE)
car::Anova(SRmodel.60.ar1, type = 3)













# GAM ---------------------------------------------------------------------

library(mgcv)
SRgam.60.ar1 <- mgcv::gamm(SR ~ s(Minute, by = Day) + Day + s(Site, bs = "re"),
           correlation = corAR1(form = ~ Minute | factor(SiteDay)),
           family = "poisson",
           method = "REML",
           data = SR.Window.60)
summary(SRgam.60.ar1$gam)

performance::check_singularity(SRgam.60.ar1)
acf(SR.Window.60$SR) # raw data is autocorrelated
acf(resid(SRgam.60.ar1$lme, type = "normalized")) # visually, we do not have autocorrelated data
pacf(resid(SRgam.60.ar1$lme, type = "normalized")) # visually, we do not have autocorrelated data
runs.test(resid(SRgam.60.ar1$lme, type = "normalized")) # we do not have autocorrelated data
lmtest::dwtest(SRgam.60.ar1$residuals) # we do not have autocorrelated data

plot(SRgam.60.ar1$gam, shade = TRUE, rug = FALSE, trans = exp, pages = 1,
     all.terms = TRUE)


par(mfrow=c(2,2))
gam.check(SRgam.60.ar1$gam)
concurvity(SRgam.60.ar1$gam, full = TRUE) # no issues with concurvity
concurvity(SRgam.60.ar1$gam, full = FALSE) # no issues with concurvity
summary(SRgam.60.ar1$gam)

options(na.action = "na.fail")
# computes marginal and conditional R^2
d.out <- MuMIn::dredge(SRmodel.60.ar1)
#, extra = list("Rsq" = function(x){performance::r2(x)}))
View(d.out)
options(na.action = "na.omit")
write.csv(d.out, "Outputs/sr-model-selection.csv")

plot(SRgam.60.ar1$gam, shade = TRUE, rug = FALSE, shift = coef(SRgam.60.ar1$gam)[1], trans = exp, pages = 1,
     all.terms = TRUE)




plot(SRgam.60.ar1$gam, shade = TRUE, rug = FALSE, residuals = TRUE, trans = exp,
     pch = 1, cex = 1, pages = 1, all.terms = FALSE)



plot(SRgam.60.ar1, seWithMean = TRUE, shift = coef(SRgam.60.ar1)[1], pages = 1)

acf(SR.Window.60$SR) # raw data is autocorrelated
acf(residuals(SRgam.60.ar1)) # visually, we do not have autocorrelated data
pacf(SRgam.60.ar1)
pacf(residuals(SRgam.60.ar1)) # visually, we do not have autocorrelated data


summary(SRgam.60.ar1)
plot(SRgam.60.ar1, resid = TRUE)
plot(SRgam.60.ar1, se = TRUE)

gam.check(SRgam.60.ar1)

plot(SRgam.60.ar1$)














# TOTAL VOCAL PREVALENCE --------------------------------------------------------

TVP.Window.60 <- read.csv("Outputs/TVP.Window.60") %>% 
  mutate(Minute = Time.Window/60,
         Site = as.factor(Site),
         Day = as.factor(Day),
         Hab1 = as.factor(Hab1),
         SiteDay = as.factor(str_c(Site, Day))) 

# OUTLIERS & NORMALITY OF RESPONSE VARIABLES -----------------------------------
# This distribution is to be expected when dealing with count data/ many zero counts
ggdensity(TVP.Window.60$TVP, xlab = "Total Vocal Prevalence")
gghistogram(TVP.Window.60$TVP, xlab = "Total Vocal Prevalence")
ggqqplot(TVP.Window.60$TVP, ylab = "Total Vocal Prevalence")
shapiro.test(TVP.Window.60$TVP) # W = 0.98492, p-value = 3.958e-09, not normal

library(ggplot2)
ggplot(data = TVP.Window.60) +
  geom_point(mapping = aes(x = Minute, y = TVP, color = Site)) +
  geom_smooth(mapping = aes(x = Minute, y = TVP)) + 
  facet_grid(~ Day)

ggplot(data = TVP.Window.60) +
  geom_point(mapping = aes(x = Minute, y = TVP, color = Day)) +
  geom_smooth(mapping = aes(x = Minute, y = TVP, color = Day)) + 
  facet_grid(~ Hab1)

# TOTAL VOCAL PREVALENCE MODEL -------------------------------------------------------

TVPmodel.60 <- glmmTMB(TVP ~ (Time.Window^2)*Day + (1 | Site),
                       data = TVP.Window.60, family = "poisson", REML = F)

# CHECKING MODEL ASSUMPTIONS -------------------------------------
# Checking for homogeneity of variance & normality of residuals
mean(residuals(TVPmodel.60)) # VERY close to 0

simulateResiduals(TVPmodel.60, plot = T, refit = F, use.u = T)
shapiro.test(residuals(TVPmodel.60)) # W = 0.99545, p-value = 0.002592, not normal
# residual plots look okay

# Checking for autocorrelation/independence
acf(TVP.Window.60$TVP) # raw data is autocorrelated
acf(residuals(TVPmodel.60)) # random effects variable does not correct for this
pacf(TVP.Window.60$TVP)
pacf(residuals(TVPmodel.60)) # visually, we do not have autocorrelated data
runs.test(residuals(TVPmodel.60)) # we have autocorrelated data
lmtest::dwtest(TVPmodel.60) # we have autocorrelated data

library(DataCombine)
TVP.Window.60.lag <- data.frame(TVP.Window.60, resid.mod = residuals(TVPmodel.60)) %>% 
  slide(Var = "resid.mod", GroupVar = "SiteDay", NewVar = "lag1", slideBy = -1) %>% 
  na.omit()

TVPmodel.60.lag <- glmmTMB(TVP ~ (Time.Window^2)*Day + lag1 + (1 | Site),
                          data = TVP.Window.60.lag, family = "poisson", REML = F)

performance::check_singularity(TVPmodel.60.lag)
acf(TVP.Window.60$TVP) # raw data is autocorrelated
acf(residuals(TVPmodel.60.lag)) # lagged residuals correct for this
pacf(TVP.Window.60$TVP)
pacf(residuals(TVPmodel.60.lag)) # visually, we do not have autocorrelated data
lmtest::dwtest(TVPmodel.60.lag) # we do NOT have autocorrelated data

# Checking for homogeneity of variance & normality of residuals
mean(residuals(TVPmodel.60.lag)) # VERY close to 0

simulateResiduals(TVPmodel.60, plot = T, refit = F, use.u = T)
shapiro.test(residuals(TVPmodel.60.lag)) # W = 0.99725, p-value = 0.06686, normal!
# residual plots look okay

summary(TVPmodel.60.lag)
as.data.frame(confint(TVPmodel.60.lag)) %>% 
  mutate(Estimate = exp(Estimate), `2.5 %` = exp(`2.5 %`), `97.5 %` = exp(`97.5 %`))
performance::r2(TVPmodel.60.lag)
car::Anova(TVPmodel.60.lag, type = 3)











# adding a first order autoregressive covariance structure to account for
# temporal autocorrelation
TVPmodel.60.ar1 <- glmmTMB(TVP ~ (Time.Window^2)*Day + (1 | Site) + 
                             ar1(factor(Time.Window) + 0 | Day:Site),
                           data = TVP.Window.60, family = "poisson")

acf(TVP.Window.60$TVP) # raw data is autocorrelated
acf(residuals(TVPmodel.60.ar1)) # visually, we do not have autocorrelated data
pacf(TVP.Window.60$TVP)
pacf(residuals(TVPmodel.60.ar1)) # visually, we do not have autocorrelated data
lmtest::dwtest(TVPmodel.60.ar1) # we have autocorrelated data

summary(TVPmodel.60.ar1)
performance::r2(TVPmodel.60.ar1)
VarCorr(TVPmodel.60.ar1, condVar = TRUE)
car::Anova(TVPmodel.60.ar1, type = 3)

performance::check_singularity(TVPmodel.60.ar1)









# MODEL SELECTION ---------------------------------------------------------
TVPmodel.60 <- glmmTMB(TVP ~ poly(Time.Window, 2)*Day + (1 | Site),
                       data = TVP.Window.60, family = "poisson", REML = F)

options(na.action = "na.fail")
# computes marginal and conditional R^2
d.out <- MuMIn::dredge(TVPmodel.60.lag, extra = list("Rsq" = function(x){performance::r2(x)}))
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
emmeans(TVPmodel.60.lag, "Day", type = "response") %>% 
  cld(Letter = "abcdefg")

emmeans(TVPmodel.60.lag, "Hab1", type = "response") %>% 
  cld(Letter = "abcdefg")

# a cool visualization of the predicted means and variances
library(ggeffects)
plot(ggpredict(TVPmodel.60.lag, terms = c("Time.Window [all]", "Day", "Hab1"),
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
  facet_grid(Day ~ Hab1)

# visualizing total encounters
ggplot(data = TE.Window.10.spp) +
  geom_boxplot(mapping = aes(x = Site:Day, y = TE)) +
  #geom_smooth(mapping = aes(x = Time.Window, y = VP, color = Site)) + 
  facet_grid(~ Hab1)

# SPECIES-SPECIFIC VOCAL PREVALENCE MODEL -------------------------------------------------------

#VPmodel.10 <- glmmTMB(VP ~ Time.Window + Day + Hab1 + Edge.Distance +
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
#VPmodel.60.5 <- glmmTMB(cbind(VP, VA) ~ Time.Window + Day + Hab1 + (1 | Site) + (1 | Species),
#                        data = VP.Window.60, family = "binomial", REML = F)
#
#VPmodel.60.6 <- glmmTMB(cbind(VP, VA) ~ poly(Time.Window, 2) + Day + Hab1 + (1 | Site) + (1 | Species),
#                        data = VP.Window.60, family = "binomial", REML = F)





VPmodel.60 <- glmmTMB(cbind(VP, VA) ~ poly(Time.Window, 2)*Day + Hab1 + (1 | Site) + (1 | Species),
                      data = VP.Window.60, family = "binomial", REML = F)

VPmodel.60.1 <- glmmTMB(cbind(VP, VA) ~ poly(Time.Window, 2)*Day + Site + (1 | Species),
                      data = VP.Window.60, family = "binomial", REML = F)

# won't converge
#VPmodel.60.2 <- glmmTMB(cbind(VP, VA) ~ poly(Time.Window, 2)*Day + Hab1 + (1 + Day | Site) + (1 + Time.Window | Species),
#                        data = VP.Window.60, family = "binomial", REML = F)

VPmodel.60.3 <- glmmTMB(cbind(VP, VA) ~ poly(Time.Window, 2)*Day + Hab1 + (1 | Site) + (1 + Time.Window | Species),
                      data = VP.Window.60, family = "binomial", REML = F)

VPmodel.60.4 <- glmmTMB(cbind(VP, VA) ~ poly(Time.Window, 2) + Day + Hab1 + (1 | Site) + (1 + Time.Window | Species),
                        data = VP.Window.60, family = "binomial", REML = F)

AIC(VPmodel.60.1, VPmodel.60.2, VPmodel.60.3, VPmodel.60)

summary(VPmodel.60)
as.data.frame(confint(VPmodel.60)) %>% 
  mutate(Estimate = exp(Estimate), `2.5 %` = exp(`2.5 %`), `97.5 %` = exp(`97.5 %`))
performance::r2(VPmodel.60)
car::Anova(VPmodel.60, type = 3)


#VPmodel.60 <- glmmTMB(cbind(VP, VA) ~ Time.Window + Day + Hab1 +
#                        (1 | Site) + (1 | Species),
#                        #(1 + Time.Window*Day | Species)
#                      data = VP.Window.60, family = "binomial", REML = F)
#
#VPmodel.60.2 <- glmmTMB(cbind(VP, VA) ~ Time.Window + Day + Hab1 +
#                        (1 | Site) + (1 + Time.Window | Species),
#                      #(1 + Time.Window*Day | Species)
#                      data = VP.Window.60, family = "binomial", REML = F)
#
#VPmodel.60.3 <- glmmTMB(cbind(VP, VA) ~ Time.Window + poly(Day, 2) + Hab1 +
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
VPmodel.60 <- glmmTMB(cbind(VP, VA) ~ poly(Time.Window, 2)*Day + Hab1 + (1 | Site) + (1 | Species),
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

