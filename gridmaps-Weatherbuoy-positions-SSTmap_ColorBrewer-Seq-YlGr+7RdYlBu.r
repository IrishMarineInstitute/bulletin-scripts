# HG edits 21/11/14
# - week number fixed

# install the required libraries; you only need to run this line once:
# install.packages(c('sp','raster','rgdal','ncdf','mapplots','shapefiles','RColorBrewer'))
library(sp)
library(raster)
library(rgdal)
library(ncdf)
library(mapplots)
library(shapefiles)
library(RColorBrewer)

# read a shapefile with the coastline
coast <- read.shapefile("")

boys <- data.frame(boy=c("M5","M4","M3","M2","M6")
                   ,lon=c(-6.701,-9.992,-10.551,-5.431,-15.881)
                   ,lat=c(51.689,54.998,51.217,53.484,53.075)
)



# date for which you want the plots (YYYY-MM-DD)
#d <- as.POSIXct("2015-11-18")
#d <- as.POSIXct("2015-11-21")


# this bit is to loop through all days
dates <- seq(as.POSIXct("2016-03-03"),as.POSIXct("2016-03-07"),by=60*60*24)
for(d in format(dates)){ d <- as.POSIXct(d)

# the directories with the grid files
dir1 <- ""
dir2 <- ""

# set the working directory (this is where the output file will be written to)
setwd('')


# get the filenames
fchl <- list.files(dir1,"chl.*.grd")
fca <- list.files(dir1,"CA.*.grd")
fsst <- list.files(dir2,"sst.*.grd")

# get the date
dchl <- as.POSIXct(substring(fchl,5,12),format="%Y%m%d")
dca <- as.POSIXct(substring(fca,4,11),format="%Y%m%d")
dsst <- as.POSIXct(substring(fsst,5,12),format="%Y%m%d")

# find the latest file before the date specified
i <- which.min(sqrt(c(d-dchl))) # use sqrt to get rid of negative numbers
fchl1 <- fchl[i]
dchl1 <- format(dchl[i],"%d %B %Y") # this is for the header in the plot
i <- which.min(sqrt(c(d-dca)))
fca1 <- fca[i]
dca1 <- paste("Week",as.numeric(strftime(dca[i]+1*60*60*24,format='%W')) + 1)
i <- which.min(sqrt(c(d-dsst)))
fsst1 <- fsst[i]
dsst1 <- format(dsst[i],"%d %B %Y")


# read the grid files
chl <- raster(file.path(dir1,fchl1))
ca <- raster(file.path(dir1,fca1))
sst <- raster(file.path(dir2,fsst1))

# tidy the data 
values(chl) <- ifelse(values(chl)==-999,NA,values(chl))
values(ca) <- ifelse(values(ca)==-999,NA,values(ca))
values(sst) <- ifelse(values(sst)==-999,NA,values(sst))

# create a function to plot the maps with legends
mapfun <- function(r,zlim,pal,zlog=TRUE,main='',contour=F){
  # r is the raster object to be plotted
  # zlim are the limits of the colour scale
  # pal is the palette of the colour scale
  # main is the title
  
  # map limits
  xlim <- c(-12,-4)
  ylim <- c(48,58)
  # margins of the legend and map
  mar1 <- c(1,2,6,0.5)
  mar2 <- c(2,2,0,0.5)

  # fix any oulying values
  values(r) <- ifelse(values(r)<=zlim[1],zlim[1]*1.0000001,values(r))
  values(r) <- ifelse(values(r)>zlim[2],zlim[2],values(r))
  
  # range of the colour scale
   if(zlog) breaks <- exp(seq(log(zlim[1]),log(zlim[2]),length=100)) else 
     breaks <- seq(zlim[1],zlim[2],length=100)
  # colours
  col <- pal(length(breaks)-1)
  # set margins for colour scale  
  
  par(mar=mar1)
  # plot scale
  image(x=breaks,z=matrix(breaks),breaks=breaks,xlim=range(breaks),col=col,xaxs='i',axes=F,ann=F,log=ifelse(zlog,'x',''))
  axis(3)
  box()
  # add title
  mtext(main,3,2.5)
  # set margins for map
  par(mar=mar2)    
  # set up plotting area
  plot(NULL,xlim=xlim,ylim=ylim,xaxs='i',yaxs='i',ann=F,axes=F)
  axis(1,-seq(4,12,by=2),paste0(seq(4,12,by=2),'°W'))
  axis(2,48:58,paste0(48:58,'°N'))
  # plot the grid file
  image(r,breaks=breaks,col=col,add=T)
  # add contours if required
  if(contour) contour(sst,add=T)
  box()
}


filename <- paste0(format(d,format="%Y%m%d"),".png")
# open graphics device
png(filename=filename,width=12,height=9,units='in',pointsize=16,res=300)
  
  # set outer margins and character size
  par(oma=c(2,2,0,0.5),cex=0.3)
  # layout of the 6 plots (3 legends + 3 maps)
  layout(matrix(1:6,nrow=2),heights=rep(c(1,5),3))
  
  # colour palette
  pal <- colorRampPalette (brewer.pal(9,'YlGn'))

# make the map
  # title
  main <- paste("Chlorophyll a (mg/m\u00b3)\n",dchl1)
  mapfun(chl,c(1,20),pal,TRUE,main,F)
  # add the coastline
  draw.shape(coast,col='white',border='black')

  # title of next plot (same colour palette)
  main <- paste("Chlorophyll a - Anomaly (mg/m\u00b3)\n",dca1)
  mapfun(ca,c(1,20),pal,TRUE,main,F)
  draw.shape(coast,col='white',border='black')

  # title of next plot
  main <- paste("Sea Surface Temperature (°C)\n",dsst1)
  pal <- colorRampPalette(rev(brewer.pal(11,'RdYlBu')))
  mapfun(sst,c(4,20),pal,FALSE,main,T)
  draw.shape(coast,col='black',border='black')
  with(boys, points(lon,lat,pch=16,cex=4,col='white'))
  with(boys, text(lon,lat,boy,font=2))

  # axis labels
  mtext('Longitude',1,0.5,outer=T)
  mtext('Latitude',2,0.5,outer=T)

  # clear memory
  gc(reset=T)

# close graphics device
dev.off()

}

