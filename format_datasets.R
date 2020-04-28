## Code to analyze LEAP 2016 MS data
## Vincent Fugère 2019-2020
## This script takes all datasets, formats them, and outputs an .RData file used by analysis scripts

rm(list=ls())

library(tidyverse)
library(magrittr)
library(readxl)
library(vegan)

devtools::source_url('https://raw.githubusercontent.com/VFugere/Rfuncs/master/utils.R')

#### Sampling schedule ####

Sampling.dates <-c('17/08/16','23/08/16','31/08/16','15/09/16','20/09/16','28/09/16')
Sampling.dates <- as.Date(Sampling.dates, format = '%d/%m/%y')
Sampling.dates <- format(Sampling.dates, '%j')
Sampling.dates <- as.numeric(Sampling.dates)
Sampling.dates <- Sampling.dates - 229 #day 1 of exp is Julian day 230

pulse.dates <- as.Date(c('22/08/2016','19/09/2016','30/09/2016'), format = '%d/%m/%Y')
pulse.dates <- format(pulse.dates, '%j')
pulse.dates <- as.numeric(pulse.dates) - 229

#### Treatments ####

treat <- read.csv('~/Google Drive/Recherche/LEAP Postdoc/2016/raw data/LEAP2016treatments.csv', stringsAsFactors = F) %>%
  rename('site' = pond,'nut' = nut.f, 'gly' = gly.lvl, 'imi' = imi.lvl, 'gly.target.ppb' = gly.conc, 'imi.target.ppb' = imi.conc) %>%
  select(-nut.ug.P,-nut.rel,-array,-nb)
treat$nut.fac <- factor(treat$nut, levels = c('low','high'))
treat$nut.num <- as.numeric(treat$nut.fac)
treat$pond.id <- 1:48

#### Nutrient and pesticide data (not included in merged dataset because < 6 samples/pond) ####

nut <- read_xlsx('~/Google Drive/Recherche/LEAP Postdoc/2016/raw data/2016nutrients.xlsx') %>%
  group_by(site,time.point,var) %>% summarize(day = mean(day), conc = mean(conc.ug.per.L, na.rm=T)) %>%
  ungroup %>% filter(var != 'SRP', time.point != 4) %>% select(-time.point) %>% spread(var,conc) %>% rename(date = day)

gly <- read_xlsx('~/Google Drive/Recherche/LEAP Postdoc/2016/raw data/Contaminants/gly_clean.xlsx') %>% select(-time.point, -gly.expected.ppb, -gly.max.ppb)
gly$date <- gly$date %>% as.Date(format = '%d.%m.%Y') %>% format('%j') %>% as.numeric
gly$date <- gly$date - 229

imi <- read_xlsx('~/Google Drive/Recherche/LEAP Postdoc/2016/raw data/Contaminants/imi_clean.xlsx') %>% select(-time.point, -notes)
imi$date <- imi$date %>% as.Date(format = '%d.%m.%Y') %>% format('%j') %>% as.numeric
imi$date <- imi$date - 229

#### HOBO data (also not merged) ####

hobodat <- read_xlsx('~/Google Drive/Recherche/LEAP Postdoc/2016/raw data/HOBO_temp_clean.xlsx') %>% 
  filter(time.num < 2065) #only keeping the first 43 days of the experiment, i.e. no phase 2

#### YSI & metab data ####

YSI <- read.csv('~/Google Drive/Recherche/LEAP Postdoc/2016/raw data/YSI_exp.csv',header=T,stringsAsFactors=F)
YSI$date <- sapply(strsplit(YSI$time, '[ ]'), function(x) x[1])
YSI$time <- sapply(strsplit(YSI$time, '[ ]'), function(x) x[2])
YSI[YSI$pond == 'F1','pond'] <- 'H3'
names(YSI)[7] <- 'site' #to be consistent with other data frames
YSI$date <- as.Date(YSI$date, format = '%d/%m/%Y')
YSI$date <- format(YSI$date, '%j')
YSI$date <- as.numeric(YSI$date)
YSI$date <- YSI$date - 229 #day 1 of exp is Julian day 230
YSI <- YSI[YSI$date < 45,]
# adding NEP
YSI$period <- as.factor(YSI$period)
YSI$TP <- as.factor(YSI$date)
YSI <- YSI %>% group_by(TP, site) %>%
  mutate(NEP = DO.mg.L[2] - DO.mg.L[1], pH.diff = pH[2] - pH[1], pH.mean = mean(pH[1],pH[2]), temp.mean = mean(temp.C[1],temp.C[2]), DO.mean = mean(DO.mg.L[1],DO.mg.L[2]), SPC.mean = mean(SPC.uS.cm[1],SPC.uS.cm[2])) %>%
  ungroup
YSI <- YSI[YSI$period != 'dusk',]
#correcting for offset in 5th time point due to calibration problems.
YSI$NEP[YSI$date == 35] <- YSI$NEP[YSI$date == 35] - 1
YSI <- YSI %>% select(date,site,NEP:SPC.mean)

#### fluoroprobe data ####

FP <- read.csv('~/Google Drive/Recherche/LEAP Postdoc/2016/raw data/fluoroprobe.csv',header=T,stringsAsFactors=F) %>% select(-X) %>% as.data.frame
FP$date <- as.Date(FP$date, format = '%d/%m/%Y')
FP$date <- format(FP$date, '%j')
FP$date <- as.numeric(FP$date)
FP$date <- FP$date - 229 #day 1 of exp is Julian day 230
FP[FP$site == 'F1','site'] <- 'H3'
FP <- filter(FP, site %!in% c('LAKE','WELL'))
FP[FP$date == 2,'date'] <- 1
FP[FP$date == 8,'date'] <- 7
FP <- FP[order(FP$date),]
FP <- select(FP, date,site,greens:cryptos,total) %>%
  filter(date %in% Sampling.dates)

#### flow cytometry data ####

FC <- read_xlsx('~/Google Drive/Recherche/LEAP Postdoc/2016/raw data/bacterial_counts.xlsx')
FC[FC$pond == 'F1','pond'] <- 'H3'
FC <- select(FC, day, pond, mean.cells.ul) %>% 
  filter(!is.na(mean.cells.ul)) %>%
  mutate('day' = as.numeric(day)) %>% 
  rename('date' = day, 'site' = pond, 'BA' = mean.cells.ul) %>%
  filter(date %in% Sampling.dates) %>% arrange(date,site)

#### ecoplates data ####

EP <- read_xlsx('~/Google Drive/Recherche/LEAP Postdoc/2016/raw data/Ecoplates/Ecoplates_clean.xlsx', sheet = 'by_substrate')
EP2 <- read_xlsx('~/Google Drive/Recherche/LEAP Postdoc/2016/raw data/Ecoplates/Ecoplates_clean.xlsx', sheet = 'by_guild') %>% select(date:Amines_amides)
EP <- left_join(EP, EP2, by = c('date','site')) %>%
  filter(date %in% Sampling.dates) %>% arrange(date,site)
rm(EP2)

#### zooplankton ####

ZOO <- read_xlsx('~/Google Drive/Recherche/LEAP Postdoc/2016/raw data/Zooplankton/LEAP2016-zoo-abundance_clean-perL_final.xlsx')
ZOO$day_adj <- ZOO$day_exp
ZOO$day_adj[ZOO$day_exp == 3] <- 1
ZOO$day_adj[ZOO$day_exp == 9] <- 7
ZOO$day_adj[ZOO$day_exp == 16] <- 15
ZOO$day_adj[ZOO$day_exp == 31] <- 30
ZOO$day_adj[ZOO$day_exp == 36] <- 35
ZOO$day_adj[ZOO$day_exp == 44] <- 43
ZOO <- select(ZOO, site:day_exp, day_adj, everything())
colnames(ZOO)[6:31] <- str_to_sentence(colnames(ZOO)[6:31])
colnames(ZOO) <- str_replace(colnames(ZOO), 'spp.', 'sp')
colnames(ZOO) <- str_replace(colnames(ZOO), 'spp', 'sp')
colnames(ZOO) <- str_replace(colnames(ZOO), 'adulte', 'adult')
colnames(ZOO) <- str_replace(colnames(ZOO), ' ', '_')

adult.species <- colnames(ZOO)[c(6:16,19:31)]
crustacean.species <- colnames(ZOO)[6:16]
rotifer.species <- colnames(ZOO)[19:31]

ZOO$richness <- specnumber(ZOO[,adult.species])
ZOO$crust.richness <- specnumber(ZOO[,crustacean.species])
ZOO$rot.richness <- specnumber(ZOO[,rotifer.species])

ZOO$alphadiv <- exp(diversity(ZOO[,adult.species]))
ZOO$crust.alphadiv <- exp(diversity(ZOO[,crustacean.species]))
ZOO$rot.alphadiv <- exp(diversity(ZOO[,rotifer.species]))

ZOO$evenness <- diversity(ZOO[,adult.species])/log(specnumber(ZOO[,adult.species]))
ZOO$crust.evenness <- diversity(ZOO[,crustacean.species])/log(specnumber(ZOO[,crustacean.species]))
ZOO$rot.evenness <- diversity(ZOO[,rotifer.species])/log(specnumber(ZOO[,rotifer.species]))

#### bind and clean ####

merged.data <- inner_join(FP,YSI, by = c('date','site')) %>%
  inner_join(FC, by = c('date','site')) %>%
  inner_join(EP, by = c('date','site')) %>%
  inner_join(select(ZOO, -date, -time_pt, -day_exp), by = c('site', 'date' = 'day_adj')) %>%
  left_join(treat, by = c('site')) %>% 
  select(-gly.target.ppb,-imi.target.ppb,-water) %>%
  mutate(nut = as.numeric(factor(nut, levels=c('low','high')))) %>%
  select(date, site, gly:pond.id, NEP:SPC.mean, BA, AWCD:Amines_amides, greens:total, total_zoo_adult:rot.evenness, Alona_sp:Monostyla_quadridentata, everything())

#### output data ####

save.image('~/Google Drive/Recherche/LEAP Postdoc/2016/MSdata.RData')
