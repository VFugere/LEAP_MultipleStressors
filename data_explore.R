## Code to analyze LEAP 2016 MS data
## Vincent Fugère 2019-2020
## This is a raw/rough data exploration script:
## Older code for UQAM LEAP talk in April 2019
## Plus some new raw code for data exploration
## Clean analyses are in script foodweb_paper.R

rm(list=ls())

library(tidyverse)
library(skimr)
library(magrittr)

library(scales)
library(shape)
library(RColorBrewer)
library(viridis)
library(plotrix)

library(corrgram)
library(lme4)
library(performance)
library(sjPlot)
library(sjlabelled)
library(sjmisc)
library(mgcv)
library(itsadug)
library(mgcViz)

devtools::source_url('https://raw.githubusercontent.com/VFugere/Rfuncs/master/vif.R')
devtools::source_url('https://raw.githubusercontent.com/VFugere/Rfuncs/master/utils.R')

#### load data ####

load('~/Google Drive/Recherche/LEAP Postdoc/2016/MSdata.RData')

#### colour palette & plotting parameters ####

# mauve <- brewer.pal(9,'BuPu')[8]
# jaune <- brewer.pal(11,'RdYlBu')[6]
# vert <- brewer.pal(11,'RdYlGn')[11]
# rouge <- brewer.pal(11,'RdYlBu')[1]
# four.cols <- c(jaune,rouge,mauve,vert)

#three.cols <- c('#058A51','#130CA6','#EB0C00')
#library(wesanderson)
#three.cols <- wes_palette('Zissou1',5)[c(3,1,5)]
#three.cols <- brewer.pal(3,'Dark2')[c(2,3,1)]
#four.cols <- c('gray80',three.cols)

# glycolfunc <- colorRampPalette(c(four.cols[1], four.cols[2]))
# gly.cols <- glycolfunc(8)
# imicolfunc <- colorRampPalette(c(four.cols[1], four.cols[3]))
# imi.cols <- imicolfunc(8)
# bothcolfunc <- colorRampPalette(c(four.cols[1], four.cols[4]))
# both.cols <- bothcolfunc(8)

pchs <- c(1,0)

gly.cols <- c('gray90',brewer.pal(9, 'Reds')[2:8])
imi.cols <- c('gray90',brewer.pal(9, 'Blues')[2:8])
both.cols <- c('gray90',brewer.pal(9, 'Greens')[2:8])

glycolfunc <- colorRampPalette(gly.cols)
imicolfunc <- colorRampPalette(imi.cols)
bothcolfunc <- colorRampPalette(both.cols)

allcols <- c(gly.cols,gly.cols,imi.cols,imi.cols,both.cols,both.cols)

##### format treatment variables for models ####

#adding ordered factor for GAMs
merged.data$o.nut <- as.ordered(merged.data$nut.fac)
#rescaling pesticides gradients from 0 to 1 to compare effect with nutrient factor
merged.data$sc.gly <- rescale(merged.data$gly, c(0,1))
merged.data$sc.imi <- rescale(merged.data$imi, c(0,1))
#adding site factor
merged.data$site.f <- as.factor(merged.data$site)
#adding date factor
merged.data$date.f <- as.factor(merged.data$date)
#reordering
merged.data <- select(merged.data, date:pond.id,o.nut:date.f,everything())

#### data exploration: correlations in the dataset ####

corrgram(merged.data, order=F, lower.panel=panel.shade,
         upper.panel=panel.pie, text.panel=panel.txt,
         main=NULL)

#hard to see because of non-linear effects. Splitting by sampling date
pdf('~/Desktop/corrgrams.pdf',height = 12,width = 12,onefile = T,pointsize = 8)
for(d in Sampling.dates){
  data.sub <- filter(merged.data, date == d) %>% select(-(date:o.nut), -site.f, -date.f)
  data.sub <- data.sub[,apply(data.sub, 2, sd) != 0]
  corrgram(data.sub, order=F, lower.panel=panel.shade,
           upper.panel=panel.pie, text.panel=panel.txt,
           main=paste0('Day ',d))
}
dev.off()

##### Time series plot ####

# vars <- c('total','BA','NEP','use','Glycogen')
# var.names <- c(expression(chlorophyll~italic(a)~(mu*g~L^-1)),
#                expression(bacterial~abundance~(cells~mu*L^-1)),
#                expression(net~change~'in'~DO~(Delta*mg~DO~L^-1)),
#                expression(C~substrates~used),
#                expression(glycogen~respiration~(OD~units~above~blanks)))
# var.log <- c(T,T,F,F,F)

vars <- c('BA','greens','NEP','use')
var.names <- c(expression(bacterial~abundance~(cells~mu*L^-1)),
               expression(chlorophyll~italic(a)~(mu*g~L^-1)),
               expression(daytime~Delta*DO~(mg~L^-1)),
               expression(C~substrates~used))
var.log <- c(T,T,F,F)

#pdf('~/Desktop/MSplots.pdf',width=9,height = 4.5,pointsize = 12,onefile = T)
pdf('~/Desktop/MSplots.pdf',width=14,height = 10,pointsize = 12,onefile = T)
par(mfrow=c(4,3),mar=c(2,2,1,1),oma=c(2.5,2.5,1,1),cex=1)
#par(mfrow=c(2,3),mar=c(2,2,1,1),oma=c(2.5,2.5,1,1),cex=1)
pchs.2 <- c(16,15)
#pchs.2 <- c(16,15)

for(v in 1:length(vars)){
  
  tmp <- merged.data[,c(colnames(merged.data)[1:13],vars[v])]
  #tmp$date.idx <- tmp$date
  tmp$date.idx <- tmp$date - 0.6
  tmp$date.idx[tmp$nut.num == 2] <- tmp$date[tmp$nut.num == 2] + 0.6
  for(letter in c('C|D','E|H','J|K')){
  #for(letter in c('C','E','K','D','H','J')){
    if(var.log[v] == T){
      plotfunctions::emptyPlot(xlim=c(1,43),ylim=range(tmp[,14]),yaxt='n',xaxt='n',ann=F, bty='l', log = 'y')
    }else{
      plotfunctions::emptyPlot(xlim=c(1,43),ylim=range(tmp[,14]),yaxt='n',xaxt='n',ann=F, bty='l')
    }
    axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
    axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
    tmp.sub <- filter(tmp, str_detect(site, letter))
    for(i in 1:16){
    #for(i in 1:8){
      pond <- unique(tmp.sub$site)[i]
      sub.sub <- filter(tmp.sub, site == pond)
      points(x=sub.sub$date.idx,y=sub.sub[,14],type='o',lwd=1.2, pch=pchs.2[sub.sub$nut.num], col=alpha(allcols[sub.sub$pond.id],1))
    }
    abline(v=pulse.dates[1:2],lty=3)
    if(letter == 'C|D'){mtext(var.names[v],side=2,line=3)}
  }
  #mtext(var.names[v],side=2,outer=T,line=1)
  #mtext('day of experiment',side=1,outer=T,line=1)
}

mtext('day of experiment',side=1,outer=T,line=1)

dev.off()


