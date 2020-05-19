## Code to analyze LEAP 2016 MS data
## Vincent Fugère 2019-2020
## This code tries out a number of SEMs

rm(list=ls())

library(tidyverse)
library(piecewiseSEM)

devtools::source_url('https://raw.githubusercontent.com/VFugere/Rfuncs/master/vif.R')
devtools::source_url('https://raw.githubusercontent.com/VFugere/Rfuncs/master/utils.R')

#### load data ####

load('~/Google Drive/Recherche/LEAP Postdoc/2016/MSdata.RData')

##### format treatment variables for models ####

#adding ordered factor for GAMs
merged.data$o.nut <- as.ordered(merged.data$nut.fac)
#rescaling pesticides gradients from 0 to 1 to compare effect with nutrient factor
merged.data$sc.gly <- scales::rescale(merged.data$gly, c(0,1))
merged.data$sc.imi <- scales::rescale(merged.data$imi, c(0,1))
#adding site factor
merged.data$site.f <- as.factor(merged.data$site)
#adding date factor
merged.data$date.f <- as.factor(merged.data$date)
#reordering
merged.data <- select(merged.data, date:pond.id,o.nut:date.f,everything())

#adding RUE
merged.data$RUE <- with(merged.data, total_zoo/total)
merged.data$RUE[55:56] <- 0


merged.data$nut.num <- merged.data$nut.num-1
merged.data$nut.num.std2 <- arm::rescale(merged.data$nut.num)
merged.data$sc.gly.std2 <- arm::rescale(merged.data$sc.gly)
merged.data$sc.imi.std2 <- arm::rescale(merged.data$sc.imi)

##### SEM for NEP, time point 5 ####

sub <- filter(merged.data, date == 15)
sub$phyto <- log10(sub$total)
sub$zoo <- log10p(sub$total_zoo)
sub$bacterio <- log10(sub$BA)
sub$delta.O2 <- log10p(sub$NEP)

sem1 <- psem(
  lm(zoo ~ sc.imi.std2 + sc.gly.std2, data = sub),
  lm(phyto ~ sc.gly.std2 + zoo, data = sub),
  lm(bacterio ~ sc.gly.std2, data = sub),
  lm(delta.O2 ~ sc.gly.std2 + phyto + bacterio + zoo, data = sub),
  data = sub
)

sem1  
basisSet(sem1)  
dSep(sem1)
fisherC(sem1)

summary(sem1)
