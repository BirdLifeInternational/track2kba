#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

## TRYING TO FIND THE BEST H VALUE ###

#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

## Steffen oppel, 10 March 2019

library(ggplot2)
library(tidyverse)


## description of how a kernel is calculated: https://cran.r-project.org/web/packages/kdensity/vignettes/tutorial.html

kernel.dens  = function(y, x, h) dnorm((y-x)/h)

### SHORT GRAPHIC EVALUATION ###

x<-0
y<- seq(0,500,0.5)
h<-seq(10,60,10)
heval<-data.frame()
for (i in h){
  dens<-kernel.dens(y,x,i)
  out<-data.frame(dist=y,use=dens, h=i)
  out$min<-min(out$dist[out$use<0.0001])
  heval<-rbind(heval,out)
}

heval %>% arrange(h,dist) %>% filter(dist<301) %>%
ggplot(aes(x=dist,y=use)) + geom_line(size=1.5) +
  facet_wrap("h") +
  geom_vline(aes(xintercept=min), col='red')



#### EXHAUSTIVE EVALUATION OF H VALUE IN RELATION TO CELL SIZE ###

x<-0
y<- seq(0,1000,0.5)
h<-seq(0.5,200,0.5)
heval<-data.frame()
for (i in h){
  dens<-kernel.dens(y,x,i)
  out<-data.frame(dist=y,use=dens)
  out<-data.frame(h=i,dist0UD=min(out$dist[out$use<0.0001]))
  heval<-rbind(heval,out)
}

### PLOT THE RELATIONSHIP BETWEEN H AND MIN CELL SIZE

heval %>% mutate(min.cellsize=2*dist0UD) %>%
  ggplot(aes(y=h,x=min.cellsize)) + geom_line(size=1.2)

### HOW TO FIND MIN H VALUE
heval<-heval %>% mutate(min.cellsize=2*dist0UD) 
summary(lm(h~min.cellsize, data=heval))


#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
####        CONCLUSION            ####
#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

# The min H parameter should be 0.1228 times the cell.size of the SpatialPixels over which kernel density is estimated.
# If the H parameter is smaller, then 99.99% of the kernel density will be within a single cell (if the location is at the centre of the cell)
