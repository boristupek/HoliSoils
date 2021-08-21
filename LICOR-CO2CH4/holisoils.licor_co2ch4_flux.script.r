#######################################################
## File for processing LICOR flux data ###############
#
# boris.tupek@luke.fi 
# June 2021

#this is an example file for processing raw data from LICOR portable gas analyser LI-COR LI-7810 (CH4, CO2, H2O)
#1) read and preprocess based on ambient concentrations from LICOR file
#2) read times of field measurements from field notes file
#3) match times from field notes and from licor to select the time series of measured concetration data
#4) calculate the fluxes for CO2 and CH4 based on the molecular weight, time, concetrations, chamber dimensions, and temperature
#5) plot some preliminary figure for comparison between treatments (e.g. management, fertilization) and
#   control plots (undisturbed forest floor gas exchange) and trenched plots (where tree roots were killed by trenching but were left to decompose inside the soil)

#   notes) 
#   i) field notes were saved in excel and converted to comma delimited csv format for easy reading by R
#   ii) there is 10 plots for the group of measurements (1:4 are plots from the trenched area, and 5:10 from the control area). 
#       If needed for exampe in Karstula FI (in plots 1:2 also aboveground understory vegetation was removed and understory roots were cut 
#       and new ingrowth was prevented by the fabric installed to 10 cm depth along the measurement pint circumference) the )
#   i) running this code would analyse the example data, produce figures to FIGS folder and save fluxes in table in OUT folder

#


rm(list=ls())


## DATA ###############

#define working path (modify this paths to folder in your pc where you extracted the code)
w.path <- "D:/LUKE/Holisoils/METHODS/"

#input data
path.in <- paste0(w.path,"LICOR-CO2CH4/IN/")
#output for figures
path.figs <- paste0(w.path,"LICOR-CO2CH4/FIGS/")
#output for results
path.results <- paste0(w.path,"LICOR-CO2CH4/OUT/")

#list licor files in the folder
lic.f <- list.files(path = path.in,
                        pattern = "data", 
                        all.files = T, #F #If FALSE, only the names of visible files are returned
                        recursive = F, full.names = T,
                        include.dirs = FALSE) #search for results file in all subdirectories
lic.f
#[1] "D:/LUKE/Holisoils/METHODS/LICOR-CO2CH4/IN/SK-DOBROC TG10-01100-2021-05-06T080000.data"

###read dobroc data ###########

#read dobroc file to make a header (1 file as in list above) 
lic.vars.units <- data.frame(read.delim(lic.f[1], header = T, skip = 5, sep="\t"))[1,]
#View(lic.vars.units)

#read licor file
d <- data.frame(read.delim(lic.f[1], header = F, skip = 7, sep="\t")) 
names(d)<-names(lic.vars.units )
#View(d)
dim(d)


library(Hmisc) 

#print basic statistics
describe(d$CO2)
describe(d$CH4)

#figure 
#plot raw versus preprocessed data
par(mfrow=c(2,2))
plot(d$CO2, type="l",main="raw CO2")
plot(d$CH4, type="l",main="raw CH4")


#date time
time.l <- as.POSIXct(paste(d$DATE,d$TIME), format="%Y-%m-%d %H:%M:%S")
d$date.time <- as.POSIXct(paste(d$DATE,d$TIME), format="%Y-%m-%d %H:%M:%S")


#remove warmaup/instability sequences
d$CO2[which(d$CO2 < 400)] <- NA
d$CO2[which(d$CO2 > 1100)] <- NA #reasonable fluxes should be lower, around 400 - 800 for these small samples
d$CH4[which(d$CH4 < 0)] <- NA
d$CH4[which(d$CH4 > 15000)] <- NA #reasonable fluxes should be lower, around 2 - 10 ppm (2 000 - 10 000 ppb) for these small samples

#add to figure for quick check of data, 
#some signal from Licor warming instability has been removed in lower panes and plotted agains time
plot(d$date.time,d$CO2, type="l",main="raw CO2")
plot(d$date.time,d$CH4, type="l",main="raw CH4")


### find changing points in concentration data ################

#notes metadata times
#load full measurements metadata details
notes.f <- list.files(path = path.in,
                    pattern = "csv", 
                    all.files = T, #F #If FALSE, only the names of visible files are returned
                    recursive = F, full.names = T,
                    include.dirs = FALSE) #search for results file in all subdirectories
notes.f

plot.times <- read.csv(notes.f[1], #paste(d.path,"Dobroc GHG field campaigns May2021.csv",sep=""), 
                  header=T, sep = ",",
                  stringsAsFactors = FALSE, na.strings=c("NA",""))

names(plot.times)
#[1] "Date"                                        
#[2] "Site..DS...dobroc.spruce..DM...dobroc.mixed."
#[3] "Trench.group"                                
#[4] "PlotID"                                      
#[5] "PlotN............trench.1.4..control.5.10."  
#[6] "LICOR.time.start"                            
#[7] "LICOR.time.end"                              
#[8] "Tchamber"                                    
#[9] "Tsoil5cm"      
#View(plot.times)
#simplify names
names(plot.times)[c(2,5)]<- c("Site","Plot")
#remove empty rows
plot.times<-plot.times[complete.cases(plot.times),]
pd <-plot.times

#extract time of measurements (recorded in licor times)
library(lubridate)

#time series
time.ms <-  as.POSIXct(paste(pd$Date,pd$LICOR.time.start), format="%d/%m/%Y %H:%M:%S")+ lubridate::seconds(15)
time.me <- as.POSIXct(paste(pd$Date,pd$LICOR.time.end), format="%d/%m/%Y %H:%M:%S")

names(pd)

times <- data.frame(id = do.call(paste, c(pd[,c("Site","Trench.group","Plot")], sep="-")),
                    DATE = pd$Date,  temperature = pd$Tchamber, time.ms,time.me)
View(times)


#https://rstudio.github.io/dygraphs/
#install.packages("dygraphs")
library(dygraphs)
library(xts)          # To make the convertion data-frame / xts format
#library(tidyverse)
#library(lubridate)library(scales)
plotRange1 <- function(xMin, xMax, yMin, yMax,rangeColor){
  polygon(x=c(xMax,xMin, xMin,xMax),y=c(yMin,yMin,yMax,yMax),
          col=alpha(rangeColor,0.7),border=NA)
}

#skip this if you accept the times as recorded #################################
manual.selection <-function(time,conc){
  print(paste("identify n number of fluxes:",length(dis.ts)))
  #manual identification
  nt<- identify(time,conc)
  if (length(nt) == length(dis.ts)){ 
    return(nt)
  } else{
    print("different number of fluxes (measured /selected): check and repeat the selection")
    if(askYesNo("check and repeat the selection (Yes = repeat, No = Exit)?")==TRUE){
      Recall(time,conc)
    }else{
      break
    }
    
  }
}

##process dates and series separately ###########

#determine only starts and end the periods after 3 mins## 
mantime.ms <- as.POSIXct(rep(NA,length(time.ms)), format="%Y-%m-%d %H:%M:%S") 
mantime.me <- as.POSIXct(rep(NA,length(time.ms)), format="%Y-%m-%d %H:%M:%S")


uniq.measured.dates <- unique(as.Date(pd$Date,format="%d/%m/%Y")) 
uniq.measured.dates
#[1] "2021-05-06" #data was test measured once 

#View(d)
#select date# note:starte from day x if previous days done)
#the looop will be useful if data consists of more days
for(D in 1:length(uniq.measured.dates)){
  #D = 1
di<-subset(d,d$DATE==uniq.measured.dates[D])
#par(mfrow=c(2,1), mar=c(2,4,2,1))
#plot(di$date.time,di$CO2, type="l",main=paste("flux separation - raw CO2:",unique(d$DATE)[D]))
#plot(di$date.time,di$CH4, type="l",main=paste("flux separation - raw CH4:",unique(d$DATE)[D]))
#View(di)

##zoom in to data
di.co2 <- cbind(di$date.time,di$CO2)
# Check type of variable
## date-time format!
#ymd_hms(di$date.time)
# Then you can create the xts necessary to use dygraph
dix.co2 <- xts(x = di$CO2, order.by = di$date.time)

#https://www.r-graph-gallery.com/318-custom-dygraphs-time-series-example.html
#review the data (support for closer zoome for identify)
#note: this opens data in "Viewer" panel in RStudio, 
#move over to see the date and time of measurements
dygraph(dix.co2) %>% dyRangeSelector()

#time series 1:3
par(mfrow=c(1,1), mar=c(2,4,2,1))
#for(s in 1:3){
s=1
times.ms <- get(c("time.ms","time.ms","time.ms")[s])
times.me <- get(c("time.me","time.me","time.me")[s])
nts.meas <- which(as.Date(times.ms)==uniq.measured.dates[D])
nte.meas <- which(as.Date(times.me)==uniq.measured.dates[D])
dis.ts <- times.ms[nts.meas]
dis.te <- times.me[nte.meas]

#add 10 min before and after
xi.min <- match(min(dis.ts, na.rm=T),di$date.time)
xi.max <- match(max(dis.te, na.rm=T),di$date.time)
if(is.na(xi.min)){
  xi.min <- 1
}
if(is.na(xi.max)){
  xi.max <- length(di$date.time)
}   
   
yi.lim <- range(di$CO2[c(xi.min:xi.max)], na.rm=T)
#draw the concentratios in time
plot(di$date.time,di$CO2, type="l",main=paste("flux separation - raw CO2 series",s,":",uniq.measured.dates[D]),
     xlim = c(min(dis.ts, na.rm=T)-3*60, max(dis.te, na.rm=T)+3*60 ), # +3*60 substract/add 3 minute
     ylim = yi.lim, 
     xaxt="n")
axis.POSIXct(side=1, at = seq(min(dis.ts, na.rm=T), max(dis.te, na.rm=T), length.out = 7), format="%H:%M:%S")
abline(v=dis.ts, col = 2, lwd = 2)
abline(v=dis.te, col = 4, lwd = 2, lty = 3)
print(paste0(uniq.measured.dates[D],", series = ",s))

#choose to either accept the starting points or to modifythem if needed
#if(askYesNo("accept defualt selection (Yes = Proceed to next, No = select manually)?")==TRUE){
#    ndis.ts <- dis.ts
#    return(ndis.ts)
#  } else {
#click on starts manually
#ndis.ts <- manual.selection(di$date.time,di$CO2) #starts
#}

ndis.ts<-dis.ts
#ndis.ts[1]<-dis.ts[1]-20 #modify start if needed
abline(v=ndis.ts, col = "orange", lwd = 2)

ndis.te <- dis.te
#ndis.te[1] <- dis.te[1]+20 #ndis.ts+2*60 #modify start if needed

#highlight the new selected periods for flux calculations
for(i in 1:length(ndis.ts)){
plotRange1(ndis.ts[i], ndis.te[i], yi.lim[1], yi.lim[2],"yellow")
}

#}
#save the plot for evaluation of selected concetrations for time series for flux calculations
dev.print(pdf, file=paste0(path.figs ,"check.selected.concentrations_",uniq.measured.dates[D],".pdf"),
          width=11, height=7, pointsize=12)

if(D ==1){
  ddis.ts = ndis.ts
  ddis.te = ndis.te
}else{
  ddis.ts <- c(ddis.ts,ndis.ts)
  ddis.te <- c(ddis.te,ndis.te)
}

}

selected.timeperiods <- data.frame(start1=ddis.ts,end1=ddis.te)


#View(selected.timeperiods)
#save 
#write.table(selected.timeperiods , file = paste(d.path, "selected.timeperiods_03.06.21.csv", sep=""),
#            sep = ",", col.names = NA,qmethod = "double")


#index meta data times with licor times
ix.s <- match(ddis.ts, time.l) # start from 30 sec
ix.e <- match(ddis.te, time.l) #2 min fluxes
#rows of  seconds of measurments
for(i in 1:length(ix.s)){
  print(i)
  #i = 1
  if( !is.na(ix.s[i])&!is.na(ix.e[i]) ){
  ix.i <- ix.s[i]:ix.e[i]
  if(i == 1){
    ix <-ix.i 
  }else{
  ix <- c(ix, ix.i)
    }
  }
}
length(ix)

#select only  the licor CHAMBER-FLUX data - MEASURED (not warm ups ...)
#raw data
dim(d)
#[1] 12084    22
#processed data
dm<-d[ix,]
dim(dm)
#[1] 4400   22

#figure
#plot only selected measured concentrations for time series calculations
par(mfrow=c(2,1))
plot(dm$CO2, type= "l",main="CO2 processed chamber flux concentrations ")
plot(dm$CH4, type= "l",main="CH4 processed chamber flux concentrations ")

#View(dm)


## 2) CO2.CH4 FLUXES (SITE.PLOT FIGURES) #######################################################################################
#read, calculate fluxes, output to tables, draw figures

# Molar mass [g/mol]####
# ch4 
M.ch4 <- 16.5
# co2 
M.co2 <- 44.01

#chamber dimensions####
chr=0.15 #radius,m 
chh = 0.315 # height,m
#chamber volume 
chvol <- pi*chr^2*chh # m3
chvol 
#chamber area
cha <- pi*chr^2 #m2
cha

############################
## flux calculation loop ###
#

#air temperature
names(pd)
pd<-pd[complete.cases(pd),]
ta <- pd$Tchamber

n.d <- dim(pd)[1]
names(selected.timeperiods)

############################
## flux calculation loop ###
#
# note:
# in code below exponential fit is used only if (1) nls converged,
# and (2) if kappa of exponential slope is significant, otherwise linear fit is used!!
#
#predefine empty table for fluxes
lic.flux <- data.frame(matrix(NA, nrow = n.d, ncol = 3))
#collumn names
names(lic.flux) <- c("date.time", "co2_gm2h","ch4_mgm2h") #ch4_ugm2h


dummy_fig <- paste(path.figs ,"check.fitted.fluxes_CO2CH4.pdf", sep = "")
pdf(dummy_fig, width=11, height=7) #open figure connection
par(mfrow=c(1,2), mar=c(5,5,2,1))

for(i in 1:n.d){ #change example of n to n.m = length of all
    #i = 1
    
    
      ix.s <- match(selected.timeperiods$start1[i], dm$date.time)
      ix.e <- match(selected.timeperiods$end1[i], dm$date.time)
      
      ix <- seq(ix.s,ix.e)
      dm.i <- dm[ix,c("date.time","CO2","CH4")]
      #convert CH4 to ppm
      dm.i$CH4 <- dm.i$CH4/1e3 
    
    
    #small data fix 1
    #1. removes duplicates
    if(length(which(duplicated(dm.i)))>0){
      dm.i <- dm.i[-which(duplicated(dm.i)),]}
    
    #dim(dm.i)
    ##View(dm.i)
    
    
    dt.1 <- dm.i$date.time[1]
    dt.e <- tail(dm.i$date.time,1)[1]
    
    #detect time difference (seconds) 
    dt.i <- as.numeric(difftime(dt.e,dt.1), units="secs")
    dtime.i <- seq(0,dt.i,1) #1 s records of concentrations over 5 min period
    #dt.i
    
    #chamber temp 
    chtemp.i <- ta[i]
    
    par(mfrow=c(1,2))
    #select 2 ghg
    for(f in 1:2){ #f = gas (co2 or ch4)
      #f = 1
      #select ghg
      dconc.if <- dm.i[,c("CO2","CH4")[f]]
      dtc.i.all <- data.frame(dt=dtime.i, dc=dconc.if)
      
      
      #exclude first and last 11 sec (start remove initial stabilization, end reduce recorded time (seconds) uncertainty)
      ex.s <- 0:10 
      n0.dc <- length(dtc.i.all$dc)
      ex.e <- (n0.dc-11):n0.dc #and optionally last 11 sec
      dtc.i.ex <- dtc.i.all[-c(ex.s, ex.e),] #,ex.e5
      dtc.i <- dtc.i.ex[complete.cases(dtc.i.ex),]
      #plot(dtc.i.all$dt,dtc.i.all$dc)
      #lines(dtc.i$dt,dtc.i$dc, col = 2, lwd=3)
      
      # fix 2
      #2. removes data after an early flux end 
      # condition, accepts last half concentrations only if they are higher than the mean of first 100
      n.dc <- length(dtc.i$dc)
      dc.ex <- tail(dtc.i$dc,n.dc/2) >=  mean(dtc.i$dc[1:100])
      dtc.i$dc[((n.dc/2):n.dc)] [which(dc.ex == FALSE)]<- NA
      dtc.i <- dtc.i[!is.na(dtc.i$dc),]
      
      M.mass <- c(M.co2,M.ch4)[f] #[g/mol]
      
      #gas LINEAR slope [ppm/s] as dt is in s
      slope.lin.if = as.numeric(lsfit(dtc.i$dt, dtc.i$dc, intercept = T)$coef)[2] 
      
      #linear model
      lm.fit.i<- lm(dtc.i$dc ~ dtc.i$dt)
      lin.ps <- as.numeric(coef(lm.fit.i))
      
      #gas EXPONENTIAL slope [ppm/s] as dt is in s
      nls.cx <- nls(dc ~ c8 + (c0-c8)*exp(-kx*dt), 
                    data= dtc.i,
                    start = c(c8 = max(dtc.i$dc, na.rm=T),
                              c0 = min(dtc.i$dc, na.rm=T),
                              kx = 0.000127),
                    control = list(maxiter = 50000, minFactor=1/2000, warnOnly=T),
                    algorithm="port", 
                    lower=c(min(dtc.i$dc),
                            min(dtc.i$dc)/2,
                            0.000001),
                    upper=c(max(dtc.i$dc)*2,
                            max(dtc.i$dc),
                            0.1))
      #check if nls model converged
      #see nls code #https://github.com/SurajGupta/r-source/blob/master/src/library/stats/R/nls.R
      nls.isconv <- nls.cx$convInfo$isConv
      
      summary(nls.cx)
      
      exp.ps <- summary(nls.cx)$coeff[,1]
      exp.kx.p <- summary(nls.cx)$coeff[3,3] #p-value of kappa
      #slope from exponential flux is Sexp = (c8-c0)K 
      x1.c8 <- exp.ps[1]
      x1.c0 <- exp.ps[2]
      x1.k  <- exp.ps[3]
      slope.exp.if <- as.numeric((x1.c0-x1.c8)*(-x1.k))
      
      exp.fun <- function(exp.ps,dt){
        c8=exp.ps[1]
        c0=exp.ps[2]
        kx=exp.ps[3]
        c8+(c0-c8)*exp(-kx*dt)
      }
      
      
      dc=dtc.i$dc
      dt=dtc.i$dt
      
      plot(dt, dc, ylim =c(min(dc),max(dc)), 
           ylab=c(paste(c("CO2","CH4","N20")[f], "concentration [ppm]", sep=" ")),
           xlab = "time [s]", main = dt.1)
      lines(dt, lin.ps[1]+lin.ps[2]*dt, col = 3, lwd=3)
      #exponential fit only if nls converged!!!
      if(nls.isconv == T){
        lines(dt, exp.fun(exp.ps,dt), col = 4, lwd= 3)
      } 
      
      
      #use exponential fit if nls model converged and kappa is significant
      if(nls.isconv == F){
        slope.if <- slope.lin.if
        selected.fit = "linear"
      } else if(exp.kx.p < 0.05){ #if kappa of exponential slope is significant use exponential, otherwise use linear fit
        slope.if <- slope.exp.if
        selected.fit = "exponential"
      } else{
        slope.if <- slope.lin.if
        selected.fit = "linear"
      }
      
      legend("topleft", paste( "fit =", selected.fit, sep = ""), bty = "n" )#show in plot what fit type is used
      
      
      #flux
      #covert units [ug/m2/s] -> *60*60 [ug/m2/h] -> /1e6 [g/m2/h]
      #note:co2 convert to [g/m2/h] and and ch4 in [mg/m2/h]
      unit.convers.f <- c(1e6,1e3)[f] #if 1 units in  [ug/m2/h]
      
      flux.if <- (chvol*M.mass*101325)/(8.314*(chtemp.i+273.15)*cha)*slope.if*60*60/unit.convers.f #*24 #multiply by 24 for results per day!!!
      legend("bottomright", 
             legend= paste(round(flux.if,3),c("[g/m2/h]","[mg/m2/h]")[f], sep =" "), 
             bty="n")
      
      #update temp fluxes
      assign(c("fco2.i","fch4.i")[f],flux.if)
      
      
    }
    
    
    #update flux table ############
    #if(i == 1){
      lic.flux[i,1] <- as.character(dt.1)
      lic.flux[i,2:3] <- c(fco2.i,fch4.i)
    #}else{
    #  lic.flux[(i+n.m),1] <- as.character(dt.1)
    #  lic.flux[(i+n.m),2:3] <- c(fco2.i,fch4.i)
    #}
    
  } #end i

dev.off() #save site.plot.date fit figure into the fig.raw folder
#note: ignore nls convergence failor messages (if not converged it was not used)


#View(lic.flux)

#merge fluxes with plots and times data 
pdflux <-cbind(pd,lic.flux)
#View(pdflux)
lic.flux <-pdflux

names(lic.flux)
View(lic.flux)

#rename site to treatment (as site represents managment)
lic.flux$Treatment <- lic.flux$Site
lic.flux$Treatment <- sub("DS","monoculture",lic.flux$Treatment)
lic.flux$Treatment <- sub("DM","mixed",lic.flux$Treatment)

ti <- which(lic.flux$Plot< 5)
ci <- which(lic.flux$Plot> 4)

lic.flux$tc <- NA
lic.flux$tc[ti]<-"trench"
lic.flux$tc[ci]<-"control"

#save fluxes table
write.table(lic.flux, 
            file = paste(path.results,
                         "dobroc.fluxes_co2ch4_may2020.csv", sep=""),
            row.names=F, na="",col.names=T, sep=",")

## plot fluxes ###############################

lf.co2<-lic.flux
lf.ch4<-lic.flux
#https://www.r-graph-gallery.com/9-ordered-boxplot.html#grouped
# Reorder varieties (group) (mixing low and high treatments for the calculations)
#new_order <- with(data, reorder(variety , note, mean , na.rm=T))
new_co2 <- with(lf.co2, reorder(Treatment, co2_gm2h, mean , na.rm=T))
new_ch4 <- with(lf.ch4, reorder(Treatment, ch4_mgm2h, mean , na.rm=T))
# Then make the boxplot, asking to use the 2 factors : Treatment (in the good order) AND gasconcentration :

## boxplots all data ##
par(mfrow=c(2,1), mar=c(3,5,1,1), oma=c(1,1,1,1))
myplot<-boxplot(co2_gm2h ~ tc*new_co2,  data = lf.co2,
        boxwex=0.5 , #ylab="sickness",
        main="SK-Dobroc (Forest managements)" , 
        col=c("slateblue1" , "tomato") ,  
        xaxt="n",ylab = "Forest flloor CO2  [g m-2 h-1]", xlab = "")
my_names <- sapply(strsplit(myplot$names , '\\.') , function(x) x[[2]] )
my_names <- my_names[seq(1 , length(my_names) , 2)]
axis(1, 
     at = seq(1.5 , 4 , 2), 
     labels = my_names , 
     tick=FALSE , cex=0.3)
# Add the grey vertical lines
for(i in seq(0.5 , 10 , 2)){ 
  abline(v=i,lty=1, col="grey")
}
# Add a legend
legend("topleft", legend = c("control", "trench"), 
       col=c("slateblue1" , "tomato"),
       pch = 15, bty = "n", pt.cex = 3, cex = 1.2,  horiz = F, inset = c(0.1, 0.1))

myplot<-boxplot(ch4_mgm2h ~ tc*new_ch4,  data = lf.ch4,
                boxwex=0.5 , #ylab="sickness",
                col=c("slateblue1" , "tomato") ,  
                xaxt="n",ylab = "Forest flloor CH4  [mg m-2 h-1]", xlab = "")
my_names <- sapply(strsplit(myplot$names , '\\.') , function(x) x[[2]] )
my_names <- my_names[seq(1 , length(my_names) , 2)]
axis(1, 
     at = seq(1.5 , 4 , 2), 
     labels = my_names , 
     tick=FALSE , cex=0.3)
# Add the grey vertical lines
for(i in seq(0.5 , 10 , 2)){ 
  abline(v=i,lty=1, col="grey")
}
# Add a legend
legend("bottomright", legend = c("control", "trench"), 
       col=c("slateblue1" , "tomato"),
       pch = 15, bty = "n", pt.cex = 3, cex = 1.2,  horiz = F, inset = c(0.1, 0.1))


dev.print(pdf, file=paste(path.figs,
                          "figure_boxplot_co2ch4_fluxes.pdf", sep=""),
          width=9, height=7, pointsize=12)

