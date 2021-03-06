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

#adding RUE
merged.data$RUE <- with(merged.data, total_zoo/total)
merged.data$RUE[55:56] <- 0

#### Plot data ####

plot.data <- merged.data
plot.data <- plot.data %>% mutate_at(vars(BA,greens,total,RUE), log10) %>%
  mutate_at(vars(NEP,cyanos:cryptos,total_zoo:Monostyla_quadridentata), log10p)
plot.data <- select(plot.data, date:date.f, pH.diff:SPC.mean, BA:use, greens:total, NEP, total_zoo:rot.evenness, RUE, everything())
vars <- colnames(plot.data)[14:ncol(plot.data)]

##### Time series (6 panels) ####

pdf('~/Desktop/MSplots_6pan.pdf',width=8,height = 6,pointsize = 12,onefile = T)
par(mfrow=c(3,2),mar=c(2,2,1,1),oma=c(2.5,2.5,1,1),cex=1)

for(v in 1:length(vars)){
  
  tmp <- plot.data[,c(colnames(plot.data)[1:13],vars[v])]
  tmp <- drop_na(tmp)
  tmp <- filter(tmp, is.finite(tmp[,14]))
  
  for(letter in c('C','D','E','H','J','K')){
    
    plotfunctions::emptyPlot(xlim=c(1,43),ylim=range(tmp[,14]),yaxt='n',xaxt='n',ann=F, bty='l')
    axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
    axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
    tmp.sub <- filter(tmp, str_detect(site, letter))
    for(i in 1:n_distinct(tmp.sub$site)){
      pond <- unique(tmp.sub$site)[i]
      sub.sub <- filter(tmp.sub, site == pond) %>% as.data.frame
      points(x=sub.sub$date,y=sub.sub[,14],type='l',lwd=1,col=alpha(allcols[sub.sub$pond.id],1))
    }
    abline(v=pulse.dates[1:2],lty=3)
  }
  mtext('day of experiment',side=1,outer=T,line=1)
  mtext(vars[v],side=2,outer=T,line=1)
}

dev.off()

##### Time series (3 panels) ####

# vars <- c('total','BA','NEP','use','Glycogen')
# var.names <- c(expression(chlorophyll~italic(a)~(mu*g~L^-1)),
#                expression(bacterial~abundance~(cells~mu*L^-1)),
#                expression(net~change~'in'~DO~(Delta*mg~DO~L^-1)),
#                expression(C~substrates~used),
#                expression(glycogen~respiration~(OD~units~above~blanks)))
# var.log <- c(T,T,F,F,F)

pdf('~/Desktop/MSplots_3pan.pdf',width=14,height = 10,pointsize = 12,onefile = T)
par(mfrow=c(4,3),mar=c(2,2,1,1),oma=c(2.5,2.5,1,1),cex=1)
pchs.2 <- c(16,15)

for(v in 1:length(vars)){
  
  tmp <- plot.data[,c(colnames(plot.data)[1:13],vars[v])]
  tmp <- drop_na(tmp)
  tmp <- filter(tmp, is.finite(tmp[,14]))
  
  tmp$date.idx <- tmp$date - 0.6
  tmp$date.idx[tmp$nut.num == 2] <- tmp$date[tmp$nut.num == 2] + 0.6
  for(letter in c('C|D','E|H','J|K')){
  #for(letter in c('C','E','K','D','H','J')){
    plotfunctions::emptyPlot(xlim=c(1,43),ylim=range(tmp[,14]),yaxt='n',xaxt='n',ann=F, bty='l')
    axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
    axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
    tmp.sub <- filter(tmp, str_detect(site, letter))
    for(i in 1:n_distinct(tmp.sub$site)){
    #for(i in 1:8){
      pond <- unique(tmp.sub$site)[i]
      sub.sub <- filter(tmp.sub, site == pond)
      points(x=sub.sub$date.idx,y=sub.sub[,14],type='o',lwd=1.2, pch=pchs.2[sub.sub$nut.num], col=alpha(allcols[sub.sub$pond.id],1))
    }
    abline(v=pulse.dates[1:2],lty=3)
    if(letter == 'C|D'){mtext(vars[v],side=2,line=3)}
  }
  #mtext(var.names[v],side=2,outer=T,line=1)
  #mtext('day of experiment',side=1,outer=T,line=1)
}

mtext('day of experiment',side=1,outer=T,line=1)

dev.off()

#### scatterplot arrays ####

pdf('~/Desktop/scatterplots.pdf',width=5.5,height = 3,pointsize = 8,onefile = T)
par(mfrow=c(3,6),mar=c(0.1,0.1,0.1,0.1),oma=c(4,4,2,0.5),cex=1,xpd=T)
pchs.2 <- c(16,15)

for(v in 1:length(vars)){
  tmp <- plot.data[,c(colnames(plot.data)[1:13],vars[v])]
  tmp <- drop_na(tmp)
  tmp <- filter(tmp, is.finite(tmp[,14]))
  ylims <- range(tmp[,14])
  for(letter in c('C|D','E|H','J|K')){
    for(date.x in Sampling.dates){
      sub <- tmp %>% filter(date == date.x, str_detect(site, letter))
      sub$pesticide <- as.numeric(str_remove(sub$site, letter))
      if(date.x == 1 & letter == 'J|K'){
        plot(y=sub[,14],x=sub[,15],xlab=NULL,ylab=NULL,bty='o',type='p',ylim=ylims,pch=pchs.2[sub$nut.num],col=allcols[sub$pond.id])
      }else if(date.x == 1 & letter != 'J|K'){
        plot(y=sub[,14],x=sub[,15],xlab=NULL,ylab=NULL,bty='o',type='p',ylim=ylims,pch=pchs.2[sub$nut.num],col=allcols[sub$pond.id],xaxt='n')
      }else if(date.x != 1 & letter == 'J|K'){
        plot(y=sub[,14],x=sub[,15],xlab=NULL,ylab=NULL,bty='o',type='p',ylim=ylims,pch=pchs.2[sub$nut.num],col=allcols[sub$pond.id],yaxt='n')
      }else{
        plot(y=sub[,14],x=sub[,15],xlab=NULL,ylab=NULL,bty='o',type='p',ylim=ylims,pch=pchs.2[sub$nut.num],col=allcols[sub$pond.id],yaxt='n',xaxt='n')
      }
      #   sub.low <- filter(sub, nut==1)
    #   lines(predict(loess(sub.low[,14]~sub.low[,15])), col=alpha(1,0.5), lwd=1.5,lty=1)
    #   sub.high <- filter(sub, nut==2)
    #   lines(predict(loess(sub.high[,14]~sub.high[,15])), col=alpha(1,0.5), lwd=1.5,lty=2)
    }
  }
  mtext(vars[v],side=2,outer=T,line=2.5,cex=1.2)
  mtext('pesticide concentration',side=1,outer=T,line=2.5,cex=1.2)
  mtext(paste('day',Sampling.dates,' '),side=3,outer=T,line=0.1,at=seq(0.1,0.93,length.out = 6),adj=0.5)
}

dev.off()

#### Scatterplots all in one panel ####

scatterplot.array <- function(var, varname, xticks=F){
  tmp <- plot.data[,c(colnames(plot.data)[1:13],var)]
  tmp <- drop_na(tmp)
  tmp <- filter(tmp, is.finite(tmp[,14]))
  ylims <- range(tmp[,14])
  for(date.x in Sampling.dates){
    sub <- tmp %>% filter(date == date.x)
    sub$pesticide <- rescale(as.numeric(str_remove(sub$site, 'C|D|E|H|J|K')),c(0,1))
    if(date.x == 1 & xticks==T){
      plot(y=sub[,14],x=sub[,15],bty='o',xlim=c(-0.05,1.05),type='p',ylim=ylims,pch=pchs.2[sub$nut.num],col=1,bg=allcols[sub$pond.id],xlab=NULL,ylab=NULL, xaxt='n')
      mtext(varname,side=2,outer=F,line=2.5,cex=1.1)
      axis(1,at=rescale(1:8,c(0,1)),labels=1:8,lwd=0,lwd.ticks = 1)
    }else if(date.x == 1 & xticks==F){
      plot(y=sub[,14],x=sub[,15],bty='o',xlim=c(-0.05,1.05),type='p',ylim=ylims,pch=pchs.2[sub$nut.num],col=1,bg=allcols[sub$pond.id],xlab=NULL,ylab=NULL, xaxt='n')
      mtext(varname,side=2,outer=F,line=2.5,cex=1.1)
    }else if(date.x != 1 & xticks==T){
      plot(y=sub[,14],x=sub[,15],bty='o',xlim=c(-0.05,1.05),type='p',ylim=ylims,pch=pchs.2[sub$nut.num],col=1,bg=allcols[sub$pond.id],xlab=NULL,ylab=NULL,yaxt='n', xaxt='n')
      axis(1,at=rescale(1:8,c(0,1)),labels=1:8,lwd=0,lwd.ticks = 1)
    }else{
      plot(y=sub[,14],x=sub[,15],bty='o',xlim=c(-0.05,1.05),type='p',ylim=ylims,pch=pchs.2[sub$nut.num],col=1,bg=allcols[sub$pond.id],xlab=NULL,ylab=NULL,yaxt='n',xaxt='n')
    }
  }
}

pdf('~/Desktop/scatterplots.pdf',width=8.5,height = 5.5,pointsize = 10)
par(mfrow=c(3,6),mar=c(0.1,0.1,0.1,0.1),oma=c(4,4,2,0.5),cex=1,xpd=T)
scatterplot.array(var='BA',varname=expression(log[10]~bact.~(cells/mu*L)))
scatterplot.array(var='total',varname=expression(log[10]~chl.~italic(a)~(mu*g/L)))
scatterplot.array(var='total_zoo',varname=expression(log[10](1+zoo)~(mu*g/L)),xticks=T)
mtext(paste('day',Sampling.dates,' '),side=3,outer=T,line=0.1,at=seq(0.1,0.93,length.out = 6),adj=0.5,cex=1.1)
mtext('pesticide nominal concentration (dose 1 to 8)',side=1,outer=T,line=2.5,cex=1.1)
dev.off()

#### scatterplot arrays with gamms ####

load('~/Google Drive/Recherche/LEAP Postdoc/2016/GAMMs.RData')

scattergam <- function(var, varname, model.name, xticks=F){
  tmp <- plot.data[,c(colnames(plot.data)[1:13],var)]
  tmp <- drop_na(tmp)
  tmp <- filter(tmp, is.finite(tmp[,14]))
  ylims <- range(tmp[,14])
  for(date.x in Sampling.dates){
    sub <- tmp %>% filter(date == date.x)
    sub$pesticide <- rescale(as.numeric(str_remove(sub$site, 'C|D|E|H|J|K')),c(0,1))
    if(date.x == 1 & xticks==T){
      plot(y=sub[,14],x=sub[,15],bty='o',xlim=c(-0.05,1.05),type='p',ylim=ylims,pch=pchs[sub$nut.num],col=allcols[sub$pond.id],xlab=NULL,ylab=NULL, xaxt='n')
      mtext(varname,side=2,outer=F,line=2.5,cex=1.1)
      axis(1,at=rescale(1:8,c(0,1)),labels=1:8,lwd=0,lwd.ticks = 1)
    }else if(date.x == 1 & xticks==F){
      plot(y=sub[,14],x=sub[,15],bty='o',xlim=c(-0.05,1.05),type='p',ylim=ylims,pch=pchs[sub$nut.num],col=allcols[sub$pond.id],xlab=NULL,ylab=NULL, xaxt='n')
      mtext(varname,side=2,outer=F,line=2.5,cex=1.1)
    }else if(date.x != 1 & xticks==T){
      plot(y=sub[,14],x=sub[,15],bty='o',xlim=c(-0.05,1.05),type='p',ylim=ylims,pch=pchs[sub$nut.num],col=allcols[sub$pond.id],xlab=NULL,ylab=NULL,yaxt='n', xaxt='n')
      axis(1,at=rescale(1:8,c(0,1)),labels=1:8,lwd=0,lwd.ticks = 1)
    }else{
      plot(y=sub[,14],x=sub[,15],bty='o',xlim=c(-0.05,1.05),type='p',ylim=ylims,pch=pchs[sub$nut.num],col=allcols[sub$pond.id],xlab=NULL,ylab=NULL,yaxt='n',xaxt='n')
    }
    plot_smooth(model.name, view="sc.gly", cond=list('date'=date.x,'sc.imi'=0,'o.nut'='low'), rm.ranef=TRUE, rug=F,col=alpha(gly.cols[8],0.7),lty=1,add=T,print.summary = F,se=0)
    plot_smooth(model.name, view="sc.gly", cond=list('date'=date.x,'sc.imi'=0,'o.nut'='high'), rm.ranef=TRUE, rug=F,col=alpha(gly.cols[8],0.7),lty=2,lwd=1.5,add=T,print.summary = F,se=0)
    plot_smooth(model.name, view="sc.imi", cond=list('date'=date.x,'sc.gly'=0,'o.nut'='low'), rm.ranef=TRUE, rug=F,col=alpha(imi.cols[8],0.7),lty=1,add=T,print.summary = F,se=0)
    plot_smooth(model.name, view="sc.imi", cond=list('date'=date.x,'sc.gly'=0,'o.nut'='high'), rm.ranef=TRUE, rug=F,col=alpha(imi.cols[8],0.7),lty=2,lwd=1.5,add=T,print.summary = F,se=0)
    fitvals <- predict.gam(model.name, newdata = list('date' = rep(date.x,200), 'sc.gly' = rep(seq(0,1,length.out=100),2), 'sc.imi' = rep(seq(0,1,length.out=100),2), 'o.nut' = c(rep('low',100),rep('high',100)), 'site.f' = rep('C1',200)), exclude = s(date,site.f))
    points(fitvals[1:100]~seq(0,1,length.out=100),type='l',col=alpha(both.cols[8],0.7),lty=1)  
    points(fitvals[101:200]~seq(0,1,length.out=100),type='l',col=alpha(both.cols[8],0.7),lty=2,lwd=1.5)  
  }
}

pdf('~/Desktop/scattergams.pdf',width=8.5,height = 5.5,pointsize = 10)
par(mfrow=c(3,6),mar=c(0.1,0.1,0.1,0.1),oma=c(4,4,2,0.5),cex=1,xpd=T)
scattergam(var='BA',varname=expression(log[10]~BA~(cells/mu*L)),model.name=ba.m)
scattergam(var='total',model.name=chla.m, varname=expression(log[10]~chl.~italic(a)~(mu*g/L)))
scattergam(var='total_zoo',model.name=zoo.m, varname=expression(log[10](1+zoo)~(mu*g/L)),xticks=T)
mtext(paste('day',Sampling.dates,' '),side=3,outer=T,line=0.1,at=seq(0.1,0.93,length.out = 6),adj=0.5,cex=1.1)
mtext('pesticide nominal concentration (dose 1 to 8)',side=1,outer=T,line=2.5,cex=1.1)
dev.off()

pdf('~/Desktop/scattergams_EF.pdf',width=8.5,height = 5.5,pointsize = 10)
par(mfrow=c(3,6),mar=c(0.1,0.1,0.1,0.1),oma=c(4,4,2,0.5),cex=1,xpd=T)
scattergam(var='use',varname=expression(C~sources~used),model.name=use.m)
scattergam(var='NEP',model.name=nep.m, varname=expression(log[10]~(1+Delta*DO)~(mu*g/L)))
scattergam(var='RUE',model.name=rue.m, varname=expression(log[10]~RUE~(mu*g/mu*g)),xticks=T)
mtext(paste('day',Sampling.dates,' '),side=3,outer=T,line=0.1,at=seq(0.1,0.93,length.out = 6),adj=0.5,cex=1.1)
mtext('pesticide nominal concentration (dose 1 to 8)',side=1,outer=T,line=2.5,cex=1.1)
dev.off()

#### data exploration: correlations in the dataset ####

less.data <- plot.data %>% select(date:rot.richness)

# corrgram(merged.data, order=F, lower.panel=panel.shade,
#          upper.panel=panel.pie, text.panel=panel.txt,
#          main=NULL)

#hard to see because of non-linear effects. Splitting by sampling date
pdf('~/Desktop/corrgrams.pdf',height = 12,width = 12,onefile = T,pointsize = 8)
for(d in Sampling.dates){
  data.sub <- filter(less.data, date == d) %>% select(-(date:o.nut), -site.f, -date.f)
  data.sub <- data.sub[,apply(data.sub, 2, sd) != 0]
  corrgram(data.sub, order=F, lower.panel=panel.shade,
           upper.panel=panel.pie, text.panel=panel.txt,
           main=paste0('Day ',d))
}
dev.off()

#### data exploration: regression trees ####

library(party)

tree.data <- select(plot.data, -site, -nut, -nut.num, -pond.id, -(o.nut:date.f), -(evenness:rot.evenness))

# pdf('~/Desktop/reg_trees.pdf',width=12,height = 8,pointsize = 8,onefile = T)
# par(mfrow=c(1,1),mar=c(4,4,1,1),oma=c(0,0,0,0),cex=1)
# 
# for(v in 5:ncol(tree.data)){
#   var <- colnames(tree.data)[v]
#   sub.dat <- tree.data %>% select(var, everything())
#   colnames(sub.dat)[1] <- 'response'
#   fit <- ctree(response ~ ., data = sub.dat, controls = ctree_control(testtype = 'MonteCarlo', maxdepth = 3))
#   plot(fit, inner_panel=node_inner(fit,pval = T), terminal_panel=node_boxplot(fit, width=0.4,fill='white',ylines=3,id=F),main=var)
# }
# 
# dev.off()

sem.data <- select(tree.data, date:NEP, BA, use, greens, diatoms, total_crustacean_adult, Nauplii, Copepodite, total_rotifer, RUE)

pdf('~/Desktop/reg_trees_sem.pdf',width=12,height = 8,pointsize = 8,onefile = T)
par(mfrow=c(1,1),mar=c(4,4,1,1),oma=c(0,0,0,0),cex=1)

for(v in 5:ncol(sem.data)){
  var <- colnames(sem.data)[v]
  sub.dat <- sem.data %>% select(var, everything())
  colnames(sub.dat)[1] <- 'response'
  fit <- ctree(response ~ ., data = sub.dat, controls = ctree_control(testtype = 'MonteCarlo', maxdepth = 3))
  plot(fit, inner_panel=node_inner(fit,pval = T), terminal_panel=node_boxplot(fit, width=0.4,fill='white',ylines=3,id=F),main=var)
}

dev.off()

##### do i need to include a diatoms vs. greens throughout? ####

par(mfrow=c(3,6),mar=c(0.1,0.1,0.1,0.1),oma=c(4,4,2,0.5),cex=1,xpd=T)
pchs.2 <- c(16,15)

tmp <- merged.data[,c(colnames(merged.data)[1:13])]
tmp$prop.diatoms <- with(merged.data, diatoms/total)
tmp <- drop_na(tmp)
ylims <- range(tmp[,14])
for(letter in c('C|D','E|H','J|K')){
  for(date.x in Sampling.dates){
    sub <- tmp %>% filter(date == date.x, str_detect(site, letter))
    sub$pesticide <- as.numeric(str_remove(sub$site, letter))
    if(date.x == 1 & letter == 'J|K'){
      plot(y=sub[,14],x=sub[,15],xlab=NULL,ylab=NULL,bty='o',type='p',ylim=ylims,pch=pchs.2[sub$nut.num],col=allcols[sub$pond.id])
    }else if(date.x == 1 & letter != 'J|K'){
      plot(y=sub[,14],x=sub[,15],xlab=NULL,ylab=NULL,bty='o',type='p',ylim=ylims,pch=pchs.2[sub$nut.num],col=allcols[sub$pond.id],xaxt='n')
    }else if(date.x != 1 & letter == 'J|K'){
      plot(y=sub[,14],x=sub[,15],xlab=NULL,ylab=NULL,bty='o',type='p',ylim=ylims,pch=pchs.2[sub$nut.num],col=allcols[sub$pond.id],yaxt='n')
    }else{
      plot(y=sub[,14],x=sub[,15],xlab=NULL,ylab=NULL,bty='o',type='p',ylim=ylims,pch=pchs.2[sub$nut.num],col=allcols[sub$pond.id],yaxt='n',xaxt='n')
    }
    #   sub.low <- filter(sub, nut==1)
    #   lines(predict(loess(sub.low[,14]~sub.low[,15])), col=alpha(1,0.5), lwd=1.5,lty=1)
    #   sub.high <- filter(sub, nut==2)
    #   lines(predict(loess(sub.high[,14]~sub.high[,15])), col=alpha(1,0.5), lwd=1.5,lty=2)
  }
}
mtext('proportion diatoms',side=2,outer=T,line=2.5,cex=1.2)
mtext('pesticide concentration',side=1,outer=T,line=2.5,cex=1.2)
mtext(paste('day',Sampling.dates,' '),side=3,outer=T,line=0.1,at=seq(0.1,0.93,length.out = 6),adj=0.5)
