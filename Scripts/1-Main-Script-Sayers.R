
# Chris Sayers
# updated February 22, 2023

#---------------------- LOADING/MERGING THE DATA -------------------------------
library(tidyverse)
library(dplyr)

# resolve namespace conflicts and creating necessary functions
select <- dplyr::select
"%nin%" <- Negate("%in%")

# pulling in habitat information and community metrics calculated in Python
SiteData <- read.csv("Spreadsheets/Site-Data.csv") %>%
  select(Site, Day, Hab1, Hab2, Hab3, Edge.Distance) %>% 
  mutate(Day = as.factor(Day), Site = as.factor(Site))

# pulling in raw annotations from Raven Pro 1.5
Annotations <- read.csv("Spreadsheets/All-Annotations.csv") %>% 
  rename(Start = Begin.Time..s., End = End.Time..s., Species = species) %>% 
  mutate(Day = as.factor(Day), Site = as.factor(Site))
 
# SPECIES RICHNESS --------------------------------------------------------

# calculating species richness at each 60 s interval between 0600-0700
# creating a new data frame to populate
SR.Annotations <- Annotations

# defining the bounds of our windows and recording length
window.start <- 0
window.end <- 60
window.length <- 60
max.time <- 3600

# initializing loop to calculate species richness at each 60 s interval between 0600-0700
while(window.start < max.time) {
  
  # create a new column in df to populate with vocal presence/absence
  SR.Annotations[, ncol(SR.Annotations) + 1] <- NA
  names(SR.Annotations)[ncol(SR.Annotations)] <- paste0(window.end)
  
  for(i in 1:nrow(SR.Annotations)) {
    if(((SR.Annotations$Start[i] >= window.start) & (SR.Annotations$End[i] < window.end)) | # annotation lies completely within window
       ((SR.Annotations$Start[i] < window.start) & (SR.Annotations$End[i] >= window.start) & (SR.Annotations$End[i] < window.end)) | # annotation only intersects lower window bound
       ((SR.Annotations$Start[i] >= window.start) & (SR.Annotations$Start[i] < window.end) & (SR.Annotations$End[i] >= window.end)) | #annotation only intersects upper window bound
       ((SR.Annotations$Start[i] < window.start) & (SR.Annotations$End[i] >= window.end))) # annotation completely overlaps window
    {SR.Annotations[i, ncol(SR.Annotations)] <- 1}
    else {SR.Annotations[i, ncol(SR.Annotations)] <- 0}
  }
  # adjusting our counter
  window.start <- window.start + window.length
  window.end <- window.end + window.length
  
  # keep track of iteration
  cat(paste("done with iteration", window.start, "\n"))
}

write.csv(SR.Annotations, "Outputs/SR.Annotations")

# calculating species richness for each time window
SR.Window <- SR.Annotations %>% 
  pivot_longer(`60`:`3600`, names_to = "Time.Window", values_to = "Presence") %>%
  ## excluding individuals that were not identified with 100% confidence
  filter(exclusion.code < 3, background != 1) %>% 
  # calculating species richness per time window per day per site
  filter(Presence == 1) %>%
  group_by(Site, Day, Time.Window) %>% 
  summarize(SR = length(unique(Species)), exclusion.code = exclusion.code)

# zero-filling data frame to represent absences
SR.filler <- SR.Annotations %>% 
  pivot_longer(`60`:`3600`, names_to = "Time.Window", values_to = "Presence") %>%
  # excluding individuals that were not identified with 100% confidence
  filter(exclusion.code < 3, background != 1) %>% 
  select(Site, Day, Time.Window) %>% 
  distinct()

# putting everything together
SR.Window.60 <- left_join(SR.filler, SR.Window, by = c("Site", "Day", "Time.Window")) %>% 
  mutate(SR = replace_na(SR, 0)) %>% 
  left_join(SiteData, by = c("Site", "Day")) %>% 
  mutate(Time.Window = as.numeric(Time.Window))

write.csv(SR.Window.60, "Outputs/SR.Window.60")

# DATA VISUALIZATION ------------------------------------------------------
library(ggpubr)
library(MuMIn)
library(glmmTMB)
library(performance)
library(car)
library(DHARMa)
library(lawstat)
library(emmeans)
library(multcomp)

# much of the strategy below is from Zurr et al. (2010)  https://doi.org/10.1111/j.2041-210X.2009.00001.x

# OUTLIERS & NORMALITY OF RESPONSE VARIABLES -----------------------------------
ggdensity(SR.Window.60$SR, xlab = "Species Richness") # looks great
gghistogram(SR.Window.60$SR, xlab = "Species Richness")
ggqqplot(SR.Window.60$SR, ylab = "Species Richness") # tails stray from normal
shapiro.test(SR.Window.60$SR) # W = 0.97989, p-value = 4.43e-11, not normal

library(ggplot2)
ggplot(data = SR.Window.60) +
  geom_point(mapping = aes(x = Time.Window, y = SR, color = Site)) +
  geom_smooth(mapping = aes(x = Time.Window, y = SR, color = Site)) + 
  facet_grid(Day ~ Hab2)

# SPECIES RICHNESS MODEL -------------------------------------------------------

SRmodel.60 <- glmmTMB(SR ~ poly(Time.Window, 2)*Day + Hab2 + Edge.Distance + (1 | Site),
                        data = SR.Window.60, family = "poisson", REML = F)

SRmodel.60 <- glmmTMB(SR ~ Time.Window*Day + Hab2 + Edge.Distance + (1 | Site),
                      data = SR.Window.60, family = "poisson", REML = F)

performance::r2(SRmodel.60)
car::Anova(SRmodel.60, type = 3)

SRmodel.60 <- glmmTMB(SR ~ Time.Window*Day + Hab2 + Edge.Distance + (1 | Site),
                     data = SR.Window.60, family = "poisson", REML = F)

summary(SRmodel.60)
as.data.frame(confint(SRmodel.60)) %>% 
  mutate(Estimate = exp(Estimate), `2.5 %` = exp(`2.5 %`), `97.5 %` = exp(`97.5 %`))
performance::r2(SRmodel.60)
car::Anova(SRmodel.60, type = 3)

# CHECKING MODEL ASSUMPTIONS -------------------------------------
# Checking for homogeneity of variance & normality of residuals
mean(residuals(SRmodel.60)) # VERY close to 0

simulateResiduals(SRmodel.60, plot = T, refit = F, use.u = T)
shapiro.test(residuals(SRmodel.60)) # W = 0.99464, p-value = 0.6516, normal!
# residual plots look okay

## Checking for autocorrelation/independence
#acf(SpeciesRichnessWindow$Species.Richness.Window) # raw data is autocorrelated
#acf(residuals(SRmodel)) # random effects variable corrects for this
#runs.test(residuals(SRmodel)) # we do not have autocorrelated data
#
#library(DataCombine)
#SpeciesRichnessWindow.lag <- data.frame(SpeciesRichnessWindow, resid.mod = residuals(SRmodel)) %>% 
#  slide(Var = "resid.mod", NewVar = "lag1", slideBy = -1) %>% 
#  na.omit()
#
#SRmodel.lag <- glmmTMB(Species.Richness.Window ~ Time.Window*Day*Site + lag1,
#                   data = SpeciesRichnessWindow.lag, family = "gaussian", REML = F)
#
#acf(SpeciesRichnessWindow.lag$Species.Richness.Window) # raw data is autocorrelated
#acf(residuals(SRmodel.lag)) # random effects variable corrects for this
#runs.test(residuals(SRmodel.lag)) # we do not have autocorrelated data


# MODEL SELECTION ---------------------------------------------------------
SRmodel.60 <- glmmTMB(SR ~ Time.Window + Day + Hab2 + (1 | Site),
                      data = SR.Window.60, family = "poisson", REML = F)

options(na.action = "na.fail")
# computes marginal and conditional R^2
d.out <- MuMIn::dredge(SRmodel.60, extra = list("Rsq" = function(x){performance::r2(x)}))
View(d.out)
options(na.action = "na.omit")
write.csv(d.out, "Outputs/sr-model-selection.csv")

# 1st place model by a long-shot (R2 = 0.61, wi = 72%)
topSRmodel <- glmmTMB(SR ~ Time.Window + Day + Hab2 + (1 | Site),
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

library(ggeffects)
plot(ggpredict(SRmodel.60, terms = c("Time.Window [all]", "Day", "Hab2"),
          type = "random", plot = T))

# VOCAL PREVALENCE --------------------------------------------------------

# calculating vocal presence/absence for each species at 10 s time windows between 0600-0700
# there are a total of 5 ways that an annotation can interact with a time window:
# (1) the annotation lies completely within the window
# (2) the annotation only intersects the lower window bound
# (3) the annotation only intersects the upper window bound
# (4) the annotation completely overlaps the window
# (5) the annotation lies completely outside the window

# creating a new data frame to populate
VP.Annotations <- Annotations

# defining the bounds of our windows and recording length
window.start <- 0
window.end <- 10
window.length <- 10
max.time <- 3600

# initializing loop to calculate vocal presence/absence at each 10-second interval
# between 0600-0700
while(window.start < max.time) {
  
  # create a new column in df to populate with vocal presence/absence
  VP.Annotations[, ncol(VP.Annotations) + 1] <- NA
  names(VP.Annotations)[ncol(VP.Annotations)] <- paste0(window.end)
  
  for(i in 1:nrow(VP.Annotations)) {
    if(((VP.Annotations$Start[i] >= window.start) & (VP.Annotations$End[i] < window.end)) | # annotation lies completely within window
       ((VP.Annotations$Start[i] < window.start) & (VP.Annotations$End[i] >= window.start) & (VP.Annotations$End[i] < window.end)) | # annotation only intersects lower window bound
       ((VP.Annotations$Start[i] >= window.start) & (VP.Annotations$Start[i] < window.end) & (VP.Annotations$End[i] >= window.end)) | #annotation only intersects upper window bound
       ((VP.Annotations$Start[i] < window.start) & (VP.Annotations$End[i] >= window.end))) # annotation completely overlaps window
    {VP.Annotations[i, ncol(VP.Annotations)] <- 1}
    else {VP.Annotations[i, ncol(VP.Annotations)] <- 0}
  }
  # adjusting our counter
  window.start <- window.start + window.length
  window.end <- window.end + window.length
  
  # keep track of iteration
  cat(paste("done with iteration", window.start, "\n"))
}

write.csv(VP.Annotations, "Outputs/VP.Annotations")

VP.Window <- VP.Annotations %>% 
  pivot_longer(`10`:`3600`, names_to = "Time.Window", values_to = "VP") %>%
  mutate(Time.Window = as.numeric(Time.Window)) %>%
  # excluding individuals that were not identified with 100% confidence
  filter(exclusion.code < 9, background != 1) %>%
  group_by(Site, Day, Species, Time.Window) %>%
  summarize(VP = max(VP)) 
  ## filtering by present species only to make the df easier to loop over
  #filter(VP == 1)

# zero-filling each site day according to all species detected at Inkaterra
# dummy df of all the site-day-time combinations
time.combo <- VP.Annotations %>% 
  pivot_longer(`10`:`3600`, names_to = "Time.Window", values_to = "VP") %>%
  mutate(Time.Window = as.numeric(Time.Window)) %>% 
  select(Site, Day, Time.Window) %>% 
  distinct()


# initializing loop to convert 10 s windows to the appropriate 60 s window
# defining the bounds of our windows and recording length
window.start <- 0
window.end <- 60
window.length <- 60
max.time <- 3600

time.combo <- time.combo %>% 
  mutate(Time.Minute = NA)

while(window.start < max.time) {
  
  for(i in 1:nrow(time.combo)) {
    if((time.combo$Time.Window[i] > window.start) & (time.combo$Time.Window[i] <= window.end))
    {time.combo$Time.Minute[i] <- window.end}
    else {time.combo$Time.Minute[i] <- time.combo$Time.Minute[i]}
  }
  # adjusting our counter
  window.start <- window.start + window.length
  window.end <- window.end + window.length
  
  # keep track of iteration
  cat(paste("done with iteration", window.start, "\n"))
}


# dummy df that represents the entire community at Inkaterra
community <- VP.Annotations %>%
  pivot_longer(`10`:`3600`, names_to = "Time.Window", values_to = "VP") %>%
  mutate(Time.Window = as.numeric(Time.Window)) %>%
  # excluding individuals that were not identified with 100% confidence
  filter(exclusion.code < 9, background != 1) %>%
  select(Species, Time.Window) %>%
  distinct()

VP.filler <- left_join(time.combo, community, by = c("Time.Window"))

# putting everything together
VP.Window.10 <- left_join(VP.filler, VP.Window, by = c("Site", "Day", "Species", "Time.Window")) %>% 
  mutate(VP = replace_na(VP, 0)) %>% 
  left_join(SiteData, by = c("Site", "Day")) %>% 
  mutate(Time.Window = as.numeric(Time.Window))
 #%>% spread(Time.Window, VP)

write.csv(VP.Window.10, "Outputs/VP.Window.10")

# creating vocal prevalence df
VP.Window.60 <- VP.Window.10 %>% 
  group_by(Site, Day, Species, Time.Minute) %>% 
  summarize(VP = sum(VP)) %>%
  left_join(SiteData, by = c("Site", "Day")) %>% 
  rename(Time.Window = Time.Minute)

write.csv(VP.Window.60, "Outputs/VP.Window.60")

# DATA VISUALIZATION ------------------------------------------------------
library(ggpubr)
library(MuMIn)
library(glmmTMB)
library(performance)
library(car)
library(DHARMa)
library(lawstat)
library(emmeans)
library(multcomp)

# much of the strategy below is from Zurr et al. (2010)  https://doi.org/10.1111/j.2041-210X.2009.00001.x

# OUTLIERS & NORMALITY OF RESPONSE VARIABLES -----------------------------------
# This distribution is to be expected when dealing with count data/ many zero counts

ggdensity(VP.Window.10$VP, xlab = "Vocal Presence") # looks great
gghistogram(VP.Window.10$VP, xlab = "Vocal Presence")
ggqqplot(VP.Window.10$VP, ylab = "Vocal Presence") # tails stray from normal
shapiro.test(VP.Window.10$VP) # W = 0.97989, p-value = 4.43e-11, not normal

ggdensity(VP.Window.60$VP, xlab = "Vocal Prevalence") # looks great
gghistogram(VP.Window.60$VP, xlab = "Vocal Prevalence")
ggqqplot(VP.Window.60$VP, ylab = "Vocal Prevalence") # tails stray from normal
shapiro.test(VP.Window.60$VP) # W = 0.97989, p-value = 4.43e-11, not normal

library(ggplot2)
ggplot(data = VP.Window.60) +
  geom_point(mapping = aes(x = Time.Window, y = VP, color = Species)) +
  #geom_smooth(mapping = aes(x = Time.Window, y = VP, color = Site)) + 
  facet_grid(Day ~ Hab2)


TE.Window <- VP.Window.10 %>% 
  group_by(Site, Day) %>% 
  mutate(TE = sum(VP))

ggplot(data = TE.Window) +
  geom_point(mapping = aes(x = Day, y = TE)) +
  #geom_smooth(mapping = aes(x = Time.Window, y = VP, color = Site)) + 
  facet_grid(Site ~ Hab2)


# VOCAL PREVALENCE MODEL -------------------------------------------------------

VPmodel.10 <- glmmTMB(VP ~ Time.Window + Day + Hab2 + Edge.Distance +
                        (1 | Site) + (1 | Species),
                      data = VP.Window.10, family = "binomial", REML = F)

VPmodel.10 <- glmmTMB(VP ~ Time.Window*Day + Hab2 + Edge.Distance + (1 | Site),
                      data = SR.Window.60, family = "binomial", REML = F)

summary(VPmodel.10)
as.data.frame(confint(VPmodel.10)) %>% 
  mutate(Estimate = exp(Estimate), `2.5 %` = exp(`2.5 %`), `97.5 %` = exp(`97.5 %`))
performance::r2(VPmodel.10)
car::Anova(VPmodel.10, type = 3)


VPmodel.60 <- glmmTMB(cbind(VP, 6) ~ Time.Window*Day + Hab2 +
                      (1 + Day | Site) + (1 + Time.Window*Day | Species)
                      ,
                      data = VP.Window.60, family = "binomial", REML = F)

summary(VPmodel.60)
as.data.frame(confint(VPmodel.60)) %>% 
  mutate(Estimate = exp(Estimate), `2.5 %` = exp(`2.5 %`), `97.5 %` = exp(`97.5 %`))
performance::r2(VPmodel.60)
car::Anova(VPmodel.10, type = 3)


# CHECKING MODEL ASSUMPTIONS -------------------------------------
# Checking for homogeneity of variance & normality of residuals
mean(residuals(VPmodel.10)) # VERY close to 0

simulateResiduals(VPmodel.10, plot = T, refit = F, use.u = T)
runs.test(residuals(VPmodel.10)) # W = 0.99464, p-value = 0.6516, normal!
# residual plots look okay

## Checking for autocorrelation/independence
#acf(SpeciesRichnessWindow$Species.Richness.Window) # raw data is autocorrelated
#acf(residuals(SRmodel)) # random effects variable corrects for this
#runs.test(residuals(SRmodel)) # we do not have autocorrelated data
#
#library(DataCombine)
#SpeciesRichnessWindow.lag <- data.frame(SpeciesRichnessWindow, resid.mod = residuals(SRmodel)) %>% 
#  slide(Var = "resid.mod", NewVar = "lag1", slideBy = -1) %>% 
#  na.omit()
#
#SRmodel.lag <- glmmTMB(Species.Richness.Window ~ Time.Window*Day*Site + lag1,
#                   data = SpeciesRichnessWindow.lag, family = "gaussian", REML = F)
#
#acf(SpeciesRichnessWindow.lag$Species.Richness.Window) # raw data is autocorrelated
#acf(residuals(SRmodel.lag)) # random effects variable corrects for this
#runs.test(residuals(SRmodel.lag)) # we do not have autocorrelated data


# MODEL SELECTION ---------------------------------------------------------
SRmodel <- glmmTMB(Species.Richness.Window ~ Time.Window*Day*Site,
                   data = SpeciesRichnessWindow, family = "gaussian", REML = F)

options(na.action = "na.fail")
# computes marginal and conditional R^2
d.out <- MuMIn::dredge(SRmodel, extra = list("Rsq" = function(x){performance::r2(x)}))
View(d.out)
options(na.action = "na.omit")
write.csv(d.out, "Outputs/sr-model-selection.csv")

# 1st place model by a long-shot (R2 = 0.92, wi = 71%)
topSRmodel <- glmmTMB(Species.Richness.Window ~ Time.Window + Day + Site + 
                        Day*Site + Time.Window*Site,
                      data = WindowData, family = "gaussian", REML = F)

summary(topSRmodel)
as.data.frame(confint(topSRmodel)) %>% 
  mutate(Estimate = exp(Estimate), `2.5 %` = exp(`2.5 %`), `97.5 %` = exp(`97.5 %`))
performance::r2(topSRmodel)
car::Anova(topSRmodel, type = 3)

# computing post-hoc comparisons to determine significant differences among the modeled means
# need to switch model to be as.factor(Time.Window) first before performing this comparison
emmeans(topSRmodel, "Time.Window", type = "response") %>% 
  cld(Letter = "abcdefg")

emmeans(topSRmodel, "Day", type = "response") %>% 
  cld(Letter = "abcdefg")

emmeans(topSRmodel, "Site", type = "response") %>% 
  cld(Letter = "abcdefg")
