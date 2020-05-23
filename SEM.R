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

merged.data$phyto <- log10(merged.data$total)
merged.data$zoo <- log10p(merged.data$total_zoo)
merged.data$bacterio <- log10(merged.data$BA)
merged.data$delta.O2 <- log10p(merged.data$NEP)
merged.data$zp <- log10(merged.data$RUE)

results <- data.frame()

for(date.x in c(7,15,35,43)){
  
  sub <- filter(merged.data, date == date.x)
  
  sem <- psem(
    lm(zoo ~ sc.gly, data = sub),
    lm(phyto ~ sc.gly + nut.num + zoo, data = sub),
    lm(delta.O2 ~ sc.gly + phyto + nut.num, data = sub),
    data = sub
  ) 
  #adding bacterio to NEP never produces a significant path (although bacterio responds to treatments)
  #reversing the path from phyto to zoo does not change anything for other paths, and it is never a significant path
  
  cat(paste('#\n#\n================================ DAY',date.x,'================================\n#\n#'))
  print(summary(sem,.progressBar = F))
  
  subres <- cbind(date.x,coefs(sem))
  results <- rbind(results,subres)
  
}

# for plotting

results$arrow.size <- 0.5
results$arrow.size[results$P.Value < 0.05]  <- scales::rescale(results$Std.Estimate[results$P.Value < 0.05], c(1,3.5))
results

##### exploring other variables ####

#rue

for(date.x in c(7,15,30,35,43)){
  sub <- filter(merged.data, date == date.x)
  sub <- sub[is.finite(sub$zp),]
  sem <- psem(
    lm(zoo ~ sc.gly*sc.imi + phyto, data = sub),
    lm(phyto ~ sc.gly + nut.num, data = sub),
    lm(zp ~ phyto + zoo, data = sub),
    data = sub
  )
  cat(paste('#\n#\n================================ DAY',date.x,'================================\n#\n#'))
  print(summary(sem,.progressBar = F))
}

# zoo bottom up (phyto-zoo path only significant at one TP, )

for(date.x in c(7,15,30,35,43)){
  sub <- filter(merged.data, date == date.x)
  sub <- sub[is.finite(sub$zp),]
  sem <- psem(
    lm(zoo ~ sc.gly*sc.imi + phyto, data = sub),
    lm(phyto ~ sc.gly + nut.num, data = sub),
    data = sub
  ) 
  cat(paste('#\n#\n================================ DAY',date.x,'================================\n#\n#'))
  print(summary(sem,.progressBar = F))
}

# zoo top-down (worsens AIC)

for(date.x in c(7,15,30,35,43)){
  sub <- filter(merged.data, date == date.x)
  sub <- sub[is.finite(sub$zp),]
  sem <- psem(
    lm(zoo ~ sc.gly*sc.imi, data = sub),
    lm(phyto ~ sc.gly + nut.num + zoo, data = sub),
    data = sub
  ) 
  cat(paste('#\n#\n================================ DAY',date.x,'================================\n#\n#'))
  print(summary(sem,.progressBar = F))
}
