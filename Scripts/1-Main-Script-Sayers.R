
# Chris Sayers
# updated March 4, 2023

# script designed to produce relevant data frames for later modeling

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


# calculating species richness for each site-day
SR.Total <- SR.Annotations %>% 
  pivot_longer(`60`:`3600`, names_to = "Time.Window", values_to = "Presence") %>%
  ## excluding individuals that were not identified with 100% confidence
  filter(exclusion.code < 3, background != 1) %>% 
  group_by(Site, Day) %>% 
  summarize(SR = length(unique(Species)))

write.csv(SR.Total, "Outputs/SR.Total")

# calculating species richness for each time window
SR.Window <- SR.Annotations %>% 
  pivot_longer(`60`:`3600`, names_to = "Time.Window", values_to = "Presence") %>%
  ## excluding individuals that were not identified with 100% confidence
  filter(exclusion.code < 3, background != 1) %>% 
  # calculating species richness per time window per day per site
  filter(Presence == 1) %>%
  group_by(Site, Day, Time.Window) %>% 
  summarize(SR = length(unique(Species)))

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
  # calculating vocal absence for binomial model
  mutate(Time.Window = as.numeric(Time.Window))
 #%>% spread(Time.Window, VP)

write.csv(VP.Window.10, "Outputs/VP.Window.10")

# creating vocal prevalence df
VP.Window.60 <- VP.Window.10 %>% 
  group_by(Site, Day, Species, Time.Minute) %>% 
  summarize(VP = sum(VP)) %>%
  left_join(SiteData, by = c("Site", "Day")) %>% 
  rename(Time.Window = Time.Minute) %>%
  # calculating vocal absence for binomial model
  mutate(VA = 6 - VP)

write.csv(VP.Window.60, "Outputs/VP.Window.60")

# total encounters per species for each site-day = sum of VP across 10 s windows
#TE.Window.10.spp <- VP.Window.10 %>% 
#  group_by(Site, Day, Species) %>% 
#  summarize(TE = sum(VP)) %>% 
#  left_join(SiteData, by = c("Site", "Day")) %>%
#  mutate(TE.inv = 3600 - TE)
#
#write.csv(TE.Window.10.spp, "Outputs/TE.Window.10.spp")


# creating TOTAL Vocal Prevalence (TVP) data frame that sums 
# every 10-s detection across species for each 60-s window
TVP.Window.60 <- VP.Window.10 %>%
  group_by(Site, Day, Time.Minute) %>% 
  summarize(TVP = sum(VP)) %>% 
  rename(Time.Window = Time.Minute) %>%
  left_join(SiteData, by = c("Site", "Day"))

write.csv(TVP.Window.60, "Outputs/TVP.Window.60")
