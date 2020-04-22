## Code to analyze LEAP 2016 MS data
## Older code for UQAM LEAP talk in April 2019
## Plus some new raw code for data exploration
## Clean analyses are in script foodweb_paper.R

rm(list=ls())
library(tidyverse)
library(scales)
library(shape)
library(readxl)
library(RColorBrewer)
library(plotrix)
library(randomcoloR)
library(skimr)
library(magrittr)
library(corrgram)
library(readxl)
library(sjPlot)
library(sjlabelled)
library(sjmisc)

devtools::source_url('https://raw.githubusercontent.com/VFugere/Rfuncs/master/vif.R')
devtools::source_url('https://raw.githubusercontent.com/VFugere/Rfuncs/master/utils.R')

#### basic data

Sampling.dates <-c('17/08/16','23/08/16','31/08/16','15/09/16','20/09/16','28/09/16')
Sampling.dates <- as.Date(Sampling.dates, format = '%d/%m/%y')
Sampling.dates <- format(Sampling.dates, '%j')
Sampling.dates <- as.numeric(Sampling.dates)
Sampling.dates <- Sampling.dates - 229 #day 1 of exp is Julian day 230

pulse.dates <- as.Date(c('22/08/2016','19/09/2016','30/09/2016'), format = '%d/%m/%Y')
pulse.dates <- format(pulse.dates, '%j')
pulse.dates <- as.numeric(pulse.dates) - 229

treat <- read.csv('~/Google Drive/Recherche/LEAP Postdoc/2016/raw data/LEAP2016treatments.csv', stringsAsFactors = F) %>%
  rename('site' = pond,'nut' = nut.f, 'gly' = gly.lvl, 'imi' = imi.lvl, 'gly.target.ppb' = gly.conc, 'imi.target.ppb' = imi.conc) %>%
  select(-nut.ug.P,-nut.rel,-array,-nb)
treat$nut.fac <- factor(treat$nut, levels = c('low','high'))
treat$nut.num <- as.numeric(treat$nut.fac)
treat$pond.id <- 1:48

#three.cols <- c('#058A51','#130CA6','#EB0C00')
#library(wesanderson)
#three.cols <- wes_palette('Zissou1',5)[c(3,1,5)]
three.cols <- brewer.pal(3,'Dark2')[c(2,3,1)]
glycolfunc <- colorRampPalette(c("gray80", three.cols[1]))
gly.cols <- glycolfunc(8)
imicolfunc <- colorRampPalette(c("gray80", three.cols[2]))
imi.cols <- imicolfunc(8)
bothcolfunc <- colorRampPalette(c("gray80", three.cols[3]))
both.cols <- bothcolfunc(8)
allcols <- c(gly.cols,gly.cols,imi.cols,imi.cols,both.cols,both.cols)
pchs <- c(16,15)

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

#### bind and clean ####

data <- inner_join(FP,YSI, by = c('date','site')) %>%
  inner_join(FC, by = c('date','site')) %>%
  inner_join(EP, by = c('date','site')) %>%
  left_join(treat, by = c('site')) %>% 
  select(-gly.target.ppb,-imi.target.ppb,-water) %>%
  mutate(nut = as.numeric(factor(nut, levels=c('low','high')))) %>%
  select(date, site, gly:pond.id, NEP:SPC.mean, greens:total, BA, AWCD:Amines_amides, everything())

#adding ordered factor for GAMs
data$o.nut <- as.ordered(data$nut.fac)

#rescaling pesticides gradients from 0 to 1 to compare effect with nutrient factor
data$sc.gly <- rescale(data$gly, c(0,1))
data$sc.imi <- rescale(data$imi, c(0,1))

#adding site factor
data$site.f <- as.factor(data$site)

#adding date factor
data$date.f <- as.factor(data$date)

#reordering
data <- select(data, date:pond.id,o.nut:date.f,everything())

#### data exploration: correlations in the dataset ####

corrgram(data, order=F, lower.panel=panel.shade,
         upper.panel=panel.pie, text.panel=panel.txt,
         main=NULL)

#hard to see because of non-linear effects. Splitting by sampling date
pdf('~/Desktop/corrgrams.pdf',height = 12,width = 12,onefile = T,pointsize = 8)
for(d in Sampling.dates){
  data.sub <- filter(data, date == d) %>% select(-date,-site)
  data.sub <- data.sub[,apply(data.sub, 2, sd) != 0]
  corrgram(data.sub, order=F, lower.panel=panel.shade,
           upper.panel=panel.pie, text.panel=panel.txt,
           main=paste0('Day ',d))
}
dev.off()


#### finding the optimal gamm ####

#all smooths + interactions
ba.model <- gam(log10(BA) ~ o.nut + s(date,k=4) + s(sc.gly, k = 4) + s(sc.imi, k = 4) + ti(date,sc.gly, k=4) + ti(date,sc.gly, by = o.nut, k=4) + ti(date,sc.imi, k=4) + ti(date,sc.imi, by = o.nut, k=4) + ti(date,sc.gly,sc.imi, k=4) + ti(date,sc.gly,sc.imi, by = o.nut, k=4) + s(date, site.f, bs='fs',k=3, m=2), data=data, method = 'REML')
fit <- fitted(ba.model)
summary(ba.model)
gam.check(ba.model)

#no 3 way interaction
ba.model2 <- gam(log10(BA) ~ o.nut + s(date,k=4) + s(sc.gly, k = 4) + s(sc.imi, k = 4) + ti(date,sc.gly, k=4) + ti(date,sc.gly, by = o.nut, k=4) + ti(date,sc.imi, k=4) + ti(date,sc.imi, by = o.nut, k=4) + ti(date,sc.gly,sc.imi, k=4) + s(date, site.f, bs='fs',k=3, m=2), data=data, method = 'REML')
fit2 <- fitted(ba.model2)
summary(ba.model2)
gam.check(ba.model2)
plot(fit2~fit)

#not eactly sure what the te() in this one does...
ba.model3 <- gam(log10(BA) ~ nut.fac + te(date,sc.gly,sc.imi, by = nut.fac, k=4) + s(date, site.f, bs='fs',k=3, m=2), data=data, method = 'REML')
summary(ba.model3)
gam.check(ba.model3)
fit3 <- fitted(ba.model3)
plot(fit3~fit)

#no smooth for overall effect of pesticides, only interactions with time. Not sure this is appropriate
ba.model4 <- gam(log10(BA) ~ o.nut + s(date,k=4) + ti(date,sc.gly, k=4) + ti(date,sc.gly, by = o.nut, k=4) + ti(date,sc.imi, k=4) + ti(date,sc.imi, by = o.nut, k=4) + ti(date,sc.gly,sc.imi, k=4) + ti(date,sc.gly,sc.imi, by = o.nut, k=4) + s(date, site.f, bs='fs',k=3, m=2), data=data, method = 'REML')
fit4 <- fitted(ba.model4)
summary(ba.model4)
gam.check(ba.model4)
plot(fit4~fit)

#te approach instead
ba.model5 <- gam(log10(BA) ~ nut.fac + te(date,sc.gly, by = nut.fac, k=4) + te(date,sc.imi, by = nut.fac, k=4) + ti(date,sc.gly,sc.imi, by = nut.fac, k=4) + s(date, site.f, bs='fs',k=3, m=2), data=data, method = 'REML')
fit5 <- fitted(ba.model5)
summary(ba.model5)
gam.check(ba.model5)
plot(fit5~fit)

#changing s to ti
ba.model6 <- gam(log10(BA) ~ o.nut + ti(date,k=4) + ti(sc.gly, k = 4) + ti(sc.imi, k = 4) + ti(date,sc.gly, k=4) + ti(date,sc.gly, by = o.nut, k=4) + ti(date,sc.imi, k=4) + ti(date,sc.imi, by = o.nut, k=4) + ti(date,sc.gly,sc.imi, k=4) + ti(date,sc.gly,sc.imi, by = o.nut, k=4) + s(date, site.f, bs='fs',k=3, m=2), data=data, method = 'REML')
summary(ba.model6)
gam.check(ba.model6)
plot(predict(ba.model6)~predict(ba.model))

AIC(ba.model,ba.model2,ba.model3,ba.model4,ba.model5,ba.model6)

#for talk
vars <- c('total','BA','NEP','use','Glycogen')
var.names <- c(expression(chlorophyll~italic(a)~(mu*g~L^-1)),
               expression(bacterial~abundance~(cells~mu*l^-1)),
               expression(NEP~(Delta*mg~DO~L^-1)),
               expression(C~substrates~used),
               expression(glycogen~respiration~(OD~units~above~blanks)))
var.log <- c(T,T,F,F,F)

# #ecoplates exploratory analysis
# vars <- names(data[13:52])
# var.names <- names(data[13:52])
# var.log <- rep(F,40)

lwds <- seq(from = 1, to = 2.5, length.out = 8)
ln.alpha <- seq(from = 0.6, to = 0.9, length.out = 8)

#pdf('~/Desktop/MSplots.pdf',width=8,height = 6,pointsize = 12,onefile = T)
par(mfrow=c(3,2),mar=c(2,2,1,1),oma=c(2.5,2.5,1,1),cex=1)

for(v in 1:length(vars)){
  
  tmp <- data[,c(colnames(data)[1:5],vars[v])]
  
  for(letter in c('C','D','E','H','J','K')){
    if(var.log[v] == T){
      plotfunctions::emptyPlot(xlim=c(1,43),ylim=range(tmp[,6]),yaxt='n',xaxt='n',ann=F, bty='l', log = 'y')
    }else{
      plotfunctions::emptyPlot(xlim=c(1,43),ylim=range(tmp[,6]),yaxt='n',xaxt='n',ann=F, bty='l')
    }
    axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
    axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
    tmp.sub <- filter(tmp, str_detect(site, letter))
    if(letter %in% c('C','D')){
      cols <- gly.cols
    }else if(letter %in% c('E','H')){
      cols <- imi.cols
    }else{
      cols <- both.cols
    }
    for(i in 1:8){
      pond <- unique(tmp.sub$site)[i]
      sub.sub <- filter(tmp.sub, site == pond)
      points(x=sub.sub$date,y=sub.sub[,6],type='l',lwd=lwds[i],col=alpha(cols[i],ln.alpha[i]))
    }
    abline(v=pulse.dates[1:2],lty=3)
  }
  mtext('day of experiment',side=1,outer=T,line=1)
  mtext(var.names[v],side=2,outer=T,line=1)
  
}

#dev.off()

#for Jesse

ba <- select(data, date:imi, BA) %>% mutate_at(vars('site':'imi'), as.factor)
hist(ba$BA) #gamma
plot(BA~date,ba)
library(lme4)
f1 <- formula(BA ~ date+gly+imi+nut)
f1 <- update(f1, . ~ . + .:date)
f1 <- update(f1, . ~ . + (1+date|site))
# m1 <- glmer(f1, ba, family=Gamma(link='log'))
# plot(m1)
# anova(m1)
# f2 <- update(f1, log(BA) ~ .)
# m2 <- lmer(f1, ba)
# summary(m2)
# plot(m2) #bad

library(brms)



m3 <- brm(f1, ba, family=Gamma(link='log'), cores=4, iter = 10000)
summary(m3)
plot(m3)
pp_check(m3)
marginal_effects(m3)
m3 -> m.fact

ba <- ba %>% mutate_at(vars('gly','imi'), as.numeric)
m3 <- brm(f1, ba, family=Gamma(link='log'), cores=8, iter = 10000)
summary(m3)
plot(m3)
pp_check(m3)
marginal_effects(m3)

#### Ecoplates ordination (in progress) ####

#com2 <- select(data, date:imi,Polymers:Amines_amides)
com2 <- select(data, date:imi,`Pyruvic acid methyl ester`:Putrescine)
com2$fill.alpha <- 0.1
com2$fill.alpha[com2$date == 43] <- 0.8
com2$pt.cex <- 1
com2$pt.cex[com2$date > 1 & com2$date < 43] <- 0.3
com2$col <- rep(c(gly.cols,gly.cols,imi.cols,imi.cols,both.cols,both.cols),6)
com2$pch <- rep(rep(c(21,22,21,22,22,21),each=8),6)
com2$site <- as.factor(com2$site)

#spe <- com2 %>% select(Polymers:Amines_amides) %>% as.matrix
spe <- com2 %>% select(`Pyruvic acid methyl ester`:Putrescine) %>% as.matrix
spe <- spe-min(spe)+0.01
spe <- log1p(spe)
row.names(spe) <- com2$site

library(vegan)
ordi <- metaMDS(spe, distance = 'bray', k = 2, autotransform = FALSE, trymax = 1000)

pdf('~/Desktop/ordi.pdf',height = 5,width = 6,pointsize = 12)
par(mfrow=c(1,1),mar=c(4,4,1,1),oma=c(0,0,0,0))
g<-ordi$points[,1:2]
plot(g[,2] ~ g[,1], type = "n",yaxt='n',xaxt='n',ann=F,xlim=c(-0.05,0.05),ylim=c(-0.05,0.05))
title(xlab='NMDS dimension 1',cex.lab=1,line = 2.5)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
title(ylab='NMDS dimension 2',cex.lab=1,line = 2.5)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
for(i in 1:nlevels(com2$site)){
  subdat <- g[row.names(g) == levels(com2$site)[i],]
  subdat <- cbind(subdat, com2[com2$site == levels(com2$site)[i],c('pch','col')])
  lncol <- com2$col[com2$site == levels(com2$site)[i]][1]
  #arrows(x0 = subdat[1,1], x1 = subdat[5,1], y0 = subdat[1,2], y1 = subdat[5,2], lty=1, lwd=1, col = alpha(lncol,1))
  lines(x = subdat[,1], y = subdat[,2], lty=1, lwd=1, col = alpha(lncol,0.5))
  #segments(x0=subdat[1,1],x1=subdat[5,1],y0=subdat[1,2],y1=subdat[5,2], col = alpha(cols[subdat$col],0.8))
  #subdat <- subdat[c(1,5),]
  #subdat$pch[1] <- subdat$pch[1] - 15
  points(subdat[6,2]~subdat[6,1], pch = subdat$pch, col=1, bg = subdat$col, cex = 1.3)
}
stress.val <- round(ordi$stress,2)
legend('topright',bty='n',legend='Stress = 0.11')
dev.off()