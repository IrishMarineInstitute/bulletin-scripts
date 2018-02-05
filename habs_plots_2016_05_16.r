# edited CC 13/07/2015
# Used Yvones new query so that it is not constrained by publication date, rather by date sample analysed in the laboratory.
# Old Query was called "HABS_AZIMUTH_GET_ALL_RESULTS_BY_DATE_PIVOT"
# New query is called "HABS_AZIMUTH_GET_ALL_RESULTS_BY_DATE_PIVOT_YVTEST"

# edited Hans [troubleshooting date range on bubble maps..]
# There were a couple of differences between the habs_plots.R and habs_plots.txt files. 
# The line:
# ndays <- round(difftime(as.POSIXct(end_date),as.POSIXct(t0),units='days'))
# was changed to: 
# ndays <- difftime(as.POSIXct(end_date),as.POSIXct(t0),units='days')
# I changed that back.
# Also in the loop for(y in 1:0){…
# ndays <- 21 causes problems because ndays is already defined above. So I changed it back to:
# numdays <- 21
# Q <- paste("",end_date1,"', ",numdays,sep='')
# That should fix it. I think only the labels, not the data, were affected.Hans

# edited CC 11/05/2015
# changed pathlength 

# edited CC 11/03/2015
# changed scientific notation for cells/L 

# edited HG 21/11/14
# new stored procedure and new spp.
# also added donut plots

# edited HG 26/01/15
# historic plots

# edited CC 15/09/2014
# New-Servervaulttwo

# Plot HABs maps and trends
#
# Based on stored procedure: HABS_AZIMUTH_GET_ALL_RESULTS_BY_DATE_PIVOT_KML_TEST
# The person running this code needs to have execute permission on that stored procedure
# (See Yvonne McFadden about this)
# 
# Hans Gerritsen April 2013



# Fill in end date here (2016-07-09 format)
end_date <- '2017-07-22'

# Set the directory for the plots, note: double backslash
# old setwd('')
setwd('')


# Thats all, now just run the rest of the code:


# Check that the required packages are installed, if not install them (only the first time you run this code)
# if R cant access the internet, navigate to C:\Program Files\R\R-X.X.X\etc
# and add the following lines to the file Rprofile.site:
#   library(utils)
#   setInternet2(TRUE)
required.packages <- c('RODBC','mapplots')
i <- match(required.packages,installed.packages()[,1])
if(is.na(sum(i))) install.packages(required.packages[is.na(i)])
library(RODBC)
library(mapplots)
data(coast)

# determine current year
year <- substring(end_date,1,4)

# jan 1st of the current year
t0 <- paste(year,'-01-01',sep='')

# number of days to extract (end-date - jan 1st)
ndays <- round(difftime(as.POSIXct(end_date),as.POSIXct(t0),units='days'))


# sql code to run
Q <- paste("",end_date,"', ",ndays,sep='')

# connect to database and extract the data
channel <- odbcDriverConnect("")
  habs <- sqlQuery(channel,Q)

  # do the same but now for the current week in the last 10 years
  habs10y <- habs
  for(y in 1:9){
    cat('Historic data for', as.numeric(year)-y,'\n')
    # because the week numbers are messy, just take the data from a 21 day window around the end_date for each year
    end_date1 <- paste(as.numeric(year)-y,substring(end_date,6,10),sep='-')
    end_date1 <-as.POSIXct(end_date1) + 7 * 60 * 60 * 24
    numdays <- 21
    Q <- paste("",end_date1,"', ",numdays,sep='')
    habs10y <- rbind(habs10y,sqlQuery(channel,Q))
  }

close(channel)


# create new column with all dinophysis species
# not needed anymore
# habs$DinophysisAll <- rowSums(habs[,c('Dinophysis acuta','Dinophysis acuminata','Dinophysis spp.')],na.rm=T)

# create new column with week number (starts on monday by default; UK convention)
# so if you add a day the week number will shift the other direction... bit tricky
# anyway adding a day results in a week starting on sunday (morning)
# finally add one week at the end because by default the week before the first monday (or sunday in this case)
# will be week 0 and you want it to be week 1, this will go wrong when the year
# starts on a sunday i think. but hey.
habs$week <- as.numeric(strftime(as.POSIXct(habs$sample_date)+1*60*60*24,format='%W')) + 1
habs10y$week <- as.numeric(strftime(as.POSIXct(habs10y$sample_date)+1*60*60*24,format='%W')) + 1

# or if you want Jan1 is start of week1
#habs$week <- as.numeric(ceiling((difftime(as.POSIXct(habs$sample_date),t0,units='days')+1)/7))


# set up a dataframe that has all possible dates and week numbers
# this is used for the label: week x (dd/mm/yy - dd/mm/yy)
weekdate <- data.frame(date=seq(as.POSIXct(t0,tz='GMT'),as.POSIXct(end_date,tz='GMT'),length=ndays+1))
#weekdate$week <- as.numeric(ceiling((difftime(weekdate$date,t0,units='days')+0)/7))
weekdate$week <- as.numeric(strftime(as.POSIXct(weekdate$date)+1*60*60*24,format='%W')) + 1


# remove any data from the previous year (if present)
habs <- subset(habs,substring(sample_date,1,4)==year)

# plot trends in toxins over time
# first set up the png file
png('plot1.png',7,10,units='in',res=300)
  # set up the number of plots and margins
  par(mfcol=c(5,3),mar=c(4,5,3,1))
  # loop through the species
  for(sp in c('Mytilus edulis','Crassostrea gigas','Ostrea edulis')){
    # loop through the toxins
    for(tox in c('AZP','DSP','PTX','ASP','PSP')){
      # column name is combination of toxin and species
      colname <- paste(tox,'in',sp)
      # define the regulatroy limit for each toxin
      if(tox=='AZP') lim <- 0.16
      if(tox=='DSP') lim <- 0.16
      if(tox=='PTX') lim <- 0.16    
      if(tox=='PSP') lim <- 800
      if(tox=='ASP') lim <- 20
      #set up the plot
      xlim=c(0,52)
      i <- which(substring(names(habs),1,3)%in%tox)
      ylim=c(0,max(lim,max(habs[,i],na.rm=T)))
      plot(NA,xlim=xlim,ylim=ylim,main=colname,xlab='Week',ylab='Concentration')
      points(habs$week,habs[,colname],cex=0.5)
      abline(h=lim,col=2,lty=2)
    # end loop through toxins
    }
  # end loop through species
  }
#close the png file
dev.off()

# plot trends in algae over time
png('plot2.png',7,7,units='in',res=300)
  par(mfcol=c(3,3),mar=c(4,5,3,1))
  for(sp in c('Dinophysis acuminata','Dinophysis acuta','All Dinophysis spp.','Karenia mikimotoi','Alexandrium spp.','All Pseudo-nitzschia spp.','Pseudo-nitzschia delicatissima complex','Pseudo-nitzschia seriata complex','Azadinium-like cells')){
    xlim=c(0,52)
    main <- ifelse(sp=='Pseudo-nitzschia delicatissima complex','P. delicatissima complex',
                   ifelse(sp=='Pseudo-nitzschia seriata complex','P. seriata complex',sp))
    plot(habs$week,habs[,sp],cex=0.5,xlim=xlim,main=main,xlab='Week',ylab='Concentration')
  }
dev.off()


# maps of each toxin
# loop through the toxins
for(tox in c('AZP','DSP','PTX','ASP','PSP')){
  # set up png file
  png(paste(tox,'.png',sep=''),5,6.5,units='in',res=300)
  # set up number of plots and margins
  par(mfrow=c(2,2),mar=c(0,0,0,0))
  # identify the column numbers for the toxin
  i <- which(substring(names(habs),1,3)%in%tox)
  # temporary function: find the maximum value, or return NA if all values are NA
  fun <- function(x) if(sum(is.na(x))==length(x)) return(NA) else  return(max(x,na.rm=T))
  # find the maximum level of the toxin over the 3 species
  level <- apply(habs[,i],1,fun)
  # sometimes there are negative values, replace these with NA
  level <- ifelse(level<0,NA,level)
  
  # set up the breakpoints for the bubble plot
  if(tox=='AZP') lim <- breaks <- c(0,0.01,0.08,0.16,0.5,1,2,5)
  if(tox=='DSP') lim <- breaks <- c(0,0.01,0.08,0.16,0.5,1,2,5)
  if(tox=='PTX') lim <- breaks <- c(0,0.01,0.08,0.16,0.5,1,2,5)  
  if(tox=='PSP') lim <- breaks <- c(0,1,400,800,1500,3000,5000)
  if(tox=='ASP') lim <- breaks <- c(0,.01,10,20,50,100,500,1000)
  
  # find the breakpoint interval for each observation
  i <- findInterval(level,breaks)
  # set up the plotting symbol (pch), colour (col) and size (cex)
  n <- length(breaks)
  pch <- c(4,rep(16,times=n))
  col <- c('black','black','orange',rep('red',times=n-2))
  cex <- c(1,seq(1.0,5,length=n))
  habs$pch <- pch[i]
  habs$col <- col[i]
  habs$cex <- cex[i]
  # set up a blank map
  basemap(c(-11,-5),c(51.25,56),axes=F,xlab='',ylab='',bg='lightgrey')  
  # work out the legend text
  legend <- c('=0',paste(paste('>=',breaks[-1],sep=''),c(paste('<',breaks[-(1:2)],sep=''),'')))
  # plot the legend
  legend('center',legend=legend,col=col,pch=pch,pt.cex=cex,bg=NULL,bty='n',y.intersp=1.8)
  # legend header
  mtext(tox,3,-1.5,cex=0.8,font=2)
  mtext(expression(paste(,'\n(',mu,'g/g)')),3,-2.5,cex=0.8)
  
  # find the most recent week
  maxweek <- max(weekdate$week)
  # loop through the most recent 3 weeks
  for(w in (maxweek-2):maxweek){
    # start and end date for the week
    date1 <- format(min(subset(weekdate,week==w)$date),'%d/%m/%Y')
    date2 <- format(max(subset(weekdate,week==w)$date),'%d/%m/%Y')
    # blank map
    basemap(c(-11,-5),c(51.25,56),axes=F,xlab='',ylab='',bg='lightgrey')   
    # add the coastline
    draw.shape(coast,border=NA,col='darkgrey')
    # plot the bubbles for the relevant week
    with(subset(habs,as.numeric(as.character(week))==w),
      points(longitude,latitude,pch=pch,col=col,cex=cex)
    )
    # plot header
    mtext(paste(' week',w),3,-1.5,font=2,adj=0)
    mtext(paste('',date1,'-',date2),3,-2.5,font=1,adj=0,cex=0.8)
  }
 # close png file
  dev.off()
# end of loop
}


# same thing for the algae
for(sp in c('Dinophysis acuminata','Dinophysis acuta','All Dinophysis spp.'
            ,'Karenia mikimotoi','Alexandrium spp.','All Pseudo-nitzschia spp.'
            ,'Pseudo-nitzschia delicatissima complex','Pseudo-nitzschia seriata complex'
            ,'Azadinium-like cells')){

  cat(sp,'\n')
  png(paste(sp,'.png',sep=''),5,6.5,units='in',res=300)
  par(mfrow=c(2,2),mar=c(0,0,0,0))
  level <- habs[,sp]
  options(scipen=10)
  if(sp=='Dinophysis acuminata') lim <- breaks <- c(0,0.04,.12,.36,1,5,10,50)*1000
  if(sp=='Dinophysis acuta') lim <- breaks <- c(0,0.04,.12,.36,1,5,10,50)*1000
  if(sp=='All Dinophysis spp.') lim <- breaks <- c(0,0.04,.12,.36,1,5,10,50)*1000
  if(sp=='Karenia mikimotoi') lim <- breaks <- c(0,0.04,5,10,50,100,500,1000)*1000
  if(sp=='Alexandrium spp.') lim <- breaks <- c(0,0.04,.12,1,5,10,50,100)*1000
  if(sp=='All Pseudo-nitzschia spp.') lim <- breaks <- c(0,0.04,5,20,50,100,500,1000)*1000
  if(sp=='Pseudo-nitzschia delicatissima complex') lim <- breaks <- c(0,0.04,5,20,50,100,500,1000)*1000
  if(sp=='Pseudo-nitzschia seriata complex') lim <- breaks <- c(0,0.04,5,20,50,100,500,1000)*1000
  if(sp=='Azadinium-like cells') lim <- breaks <- 
c(0,0.04,.12,1,10,50,100,1000)*1000

  main <- ifelse(sp=='Pseudo-nitzschia delicatissima complex','P. delicatissima complex',
                 ifelse(sp=='Pseudo-nitzschia seriata complex','P. seriata complex',sp))
  
  
  i <- findInterval(level,breaks)
  n <- length(breaks)
  pch <- c(4,rep(16,times=n))
  col <- c('black',rep('red',times=n))
  cex <- c(1,seq(0.5,5,length=n))
  habs$pch <- pch[i]
  habs$col <- col[i]
  habs$cex <- cex[i]
  basemap(c(-11,-5),c(51.25,56),axes=F,xlab='',ylab='',bg='lightgrey')   
  legend <- c('=0',paste(paste('>=',breaks[-1],sep=''),c(paste('<',breaks[-(1:2)],sep=''),'')))
  legend('center',legend=legend,col=col,pch=pch,pt.cex=cex,bg=NULL,bty='n',y.intersp=1.8)
  mtext(main,3,-1.5,cex=0.8,font=2)
  mtext('(cells per litre)',3,-2.5,cex=0.8)  
  
  maxweek <- max(weekdate$week)
  for(w in (maxweek-2):maxweek){
    date1 <- format(min(subset(weekdate,week==w)$date),'%d/%m/%Y')
    date2 <- format(max(subset(weekdate,week==w)$date),'%d/%m/%Y')
    basemap(c(-11,-5),c(51.25,56),axes=F,xlab='',ylab='',bg='lightgrey')   
    draw.shape(coast,border=NA,col='darkgrey')
    with(subset(habs,as.numeric(as.character(week))==w),
      points(longitude,latitude,pch=pch,col=col,cex=cex)
    )
  mtext(paste(' week',w),3,-1.5,font=2,adj=0)
  mtext(paste('',date1,'-',date2),3,-2.5,font=1,adj=0,cex=0.8)
  }
  dev.off()
}
  
#### donut plots

# read the looukup site table (you can edit this in excel if needed)
lut <- read.csv('SiteLookupTable.csv')

# merge the two data frames and use data for the current week only
habs1 <- merge(lut,subset(habs,week==maxweek),by.x='Site',by.y='location',all.y=T)
# this is the number of sites that could not be linked, i.e. the site could not be found in the lookup table
if(sum(is.na(habs1$Site))>0)
   warning('At least one site could not be linked to the lookup table!') else
   message('All sites could be linked to the lookup table')

# turn it into an ordered factor (so it plots it in the right order)
habs1$Geographic.Location_ASIMUTH <- ordered(habs1$Geographic.Location_ASIMUTH,levels=c('north','east','south','southwest','west'))

# loop through the toxins
for(tox in c('AZP','DSP','PTX','ASP','PSP')){

  # define the regulatroy limit for each toxin
  if(tox=='AZP') lim <- 0.16
  if(tox=='DSP') lim <- 0.16
  if(tox=='PTX') lim <- 0.16    
  if(tox=='PSP') lim <- 800
  if(tox=='ASP') lim <- 20

  # make a table with number of samples below and above threshold
  # first find columns with toxin
  i <- grep(tox,names(habs1)  )
  # then stack those columns
  habs2 <- data.frame(area=habs1$'Geographic.Location_ASIMUTH',Site=habs1$Site,stack(habs1[,i]))
  
  # this is to be able to count the sites
  habs2$site1 <- ifelse(!is.na(habs2$values),habs2$Site,NA)
  habs2$site2 <- ifelse(habs2$values>=lim,habs2$Site,NA)
  fun <- function(x) length(unique(na.omit(x)))
  table1 <- with(habs2,rbind(
    tapply(site1,list(area),fun)-tapply(site2,list(area),fun)
    ,tapply(site2,list(area),fun)
    ))
  table1
  # background colour for plot
  col.bg <- 'dodgerblue'
  # text colour
  col.txt <- 'white'
  # colours of the slices
  col.pie <- c('lightgrey','red')
  
  png(paste(tox,'_donut.png',sep=''),3,2.7,'in',8,res=300)
    # plotting area with zero margins and background colour
    par(mar=c(0,0,0,0),bg=col.bg)
    # blank plot
    plot(NA,xlim=0:1,ylim=c(0.1,1),asp=1,ann=F,axes=F)
    # check if there is any data (all areas combined)
    if(sum(table1>0)) {
      # draw pie plot
      add.pie(rowSums(table1),0.25,0.5,NA,0.25,col=col.pie,border=NA)
      # add big circle in the middle to make it into a donut
      points(0.25,0.5,cex=18,pch=16,col=col.bg)
    } else {
      # draw empty pie plot
      add.pie(1,0.25,0.5,NA,0.25,col=col.bg,border='white')
      points(0.25,0.5,cex=18,pch=21,bg=col.bg,col='white')
    }
    # add the text for number of samples etc
    text(0.25,0.5,paste(sum(table1),'sites',sep='\n'),cex=2,col=col.txt,adj=c(0.5,0.25))
    text(0.25,0.5+.3,'All areas',cex=1.66,col=col.txt)
  
    # separate donuts for each area
    for(i in 1:5){
      # position of the pie plots
      x <- c(0.775,0.9,0.9,0.65,0.65)[i]
      y <- c(0.75,0.5,0.25,0.25,0.5)[i]
      # check if there is data etc.
      if(sum(table1[,i])>0) {
        add.pie(table1[,i],x,y,NA,0.09,col=col.pie,border=NA) 
        points(x,y,cex=6,pch=16,col=col.bg)
        } else {
        add.pie(1,x,y,NA,0.09,col=col.bg,border='white')
        points(x,y,cex=6,pch=21,bg=col.bg,col='white')
        }
        # ass the text
        text(x,y,sum(table1[,i]),cex=1.5,col=col.txt)
        text(x,y+0.125,colnames(table1)[i],cex=1.25,col=col.txt)          
    # end of loop trough the 5 areas
    }    
    # add main title, week number and legend
    text(0,1,tox,cex=4,col=col.txt,adj=c(0,1))
    text(1,1,paste('Week',habs1$week[1]),cex=1.25,col=col.txt,adj=c(1,1))
    labels <- paste(c('Not toxic (','Toxic ('),round(100*rowSums(table1)/sum(table1)),c('%)','%)'),sep='')
    legend(0,0.25,labels,fill=col.pie,bty='n',text.col=col.txt,cex=1.25,border=NA)
  dev.off()
}


#### historic plots

# subset data for current weeknumber only
maxweek <- max(weekdate$week)
habs10y1 <- merge(lut,subset(habs10y,week==maxweek),by.x='Site',by.y='location',all.y=T)

# find column names that have the toxins
i <- grep(' in ',names(habs10y1)  )
# find the column names that we want to keep 
j <- which(names(habs10y1) %in% c('Site','Geographic.Location_ASIMUTH','sample_date','week'))
# stack the data
habs10y2 <- cbind(habs10y1[,j],stack(habs10y1[,i]))
# work out which toxin and the limits
habs10y2$tox <- substring(habs10y2$ind,1,3)
habs10y2$lim <- NA
habs10y2$lim <- ifelse(habs10y2$tox=='AZP',0.16,habs10y2$lim)
habs10y2$lim <- ifelse(habs10y2$tox=='DSP',0.16,habs10y2$lim)
habs10y2$lim <- ifelse(habs10y2$tox=='PTX',0.16,habs10y2$lim)
habs10y2$lim <- ifelse(habs10y2$tox=='PSP',800,habs10y2$lim)
habs10y2$lim <- ifelse(habs10y2$tox=='ASP',20,habs10y2$lim)
# get year from sample date
habs10y2$year <- as.numeric(substring(habs10y2$sample_date,1,4))

# trick to work out sites that were teseted tested and sites that tested positive
habs10y2$site1 <- ifelse(!is.na(habs10y2$values),habs10y2$Site,NA)
habs10y2$site2 <- ifelse(habs10y2$values>=habs10y2$lim,habs10y2$Site,NA)

# function to plot the history plots, loc is location, e.g. 'north', 'east' etc
thisweek <- function(loc) {
  # function to count the sites
  fun <- function(x) length(unique(na.omit(x)))
  # table with all tested sites
  t1a <- with(subset(habs10y2,Geographic.Location_ASIMUTH==loc),
              tapply(site1,list(year,tox),fun))
  # table with all positive sites
  t1b <- with(subset(habs10y2,Geographic.Location_ASIMUTH==loc),
              tapply(site2,list(year,tox),fun))
  # proportions (for colour)
  t1 <- ifelse(t1a==0,-1,t1b/t1a)
  # break points for colours
  breaks <- c(-1,-1e-6,seq(0,1,by=0.1))
  breaks[3] <- 1e-6
  col <- c('grey','lightgreen',rev(heat.colors(10)))
  # rotate axis labels
  par(las=2)
  # 'heatmap'
  image(x=1:nrow(t1),y=1:ncol(t1),z=t1,breaks=breaks,col=col,axes=F,ann=F)
  axis(1,at=1:nrow(t1),labels=rownames(t1))
  axis(2,at=1:ncol(t1),labels=colnames(t1))
  # axis labels back to normal
  par(las=0)
  # headers
  mtext('sites toxic / sites measured',3,1)
  mtext(loc,3,2.5,font=2,cex=1.5)
  # grid lines
  abline(v=(0:nrow(t1))+0.5,col='white')
  abline(h=(0:ncol(t1))+0.5,col='white')
  # box around plot
  box()
  # data frame for text in plot
  d1 <- data.frame(expand.grid(x=1:nrow(t1),y=1:ncol(t1)),total=c(t1a),toxic=c(t1b))
  d1$lab <- with(d1,ifelse(total==0,NA,paste0(toxic,'/',total)))
  with(d1,text(x,y,lab,cex=0.9))
}

# save the plot as png
png('HistoricPlots.png',6,7,'in',10,res=600)
  par(mar=rep(0.1,4))
  # background map
  basemap(c(-11,-5),c(51.25,56),axes=F,xlab='',ylab='',bg='lightgrey')  
  draw.shape(coast,border=NA,col='darkgrey')
  # titles
  mtext(paste('Week',maxweek),3,-3,cex=2.5,font=2,adj=0.93)
  mtext(paste('in the last 10 years'),3,-5,cex=2,font=2,adj=0.9)
  # colour legend
  legend('topright',c('No samples','Not toxic','Some toxic sites','Many toxic sites','All sites toxic'),fill=c('grey','lightgreen',rev(heat.colors(3))),inset=c(0.05,0.15),bg='lightgrey')
  # new=T allows you to plot on top of existing plot; layout positions plot
  par(mar=c(5,4,4,2),new=T)
  layout(matrix(c(0,0,0,0,1,0,0,0,0),ncol=3),c(0,50,50),c(5,30,65))
  thisweek('north')
  par(new=T)
  layout(matrix(c(0,0,0,0,1,0,0,0,0),ncol=3),c(50,50,0),c(35,30,35))
  thisweek('east')
  par(new=T)
  layout(matrix(c(0,0,0,0,1,0,0,0,0),ncol=3),c(0,50,50),c(35,30,35))
  thisweek('west')
  par(new=T)
  layout(matrix(c(0,0,0,0,1,0,0,0,0),ncol=3),c(0,50,50),c(65,30,5))
  thisweek('southwest')
  par(new=T)
  layout(matrix(c(0,0,0,0,1,0,0,0,0),ncol=3),c(50,50,0),c(65,30,5))
  thisweek('south')
dev.off()
