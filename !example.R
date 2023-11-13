library(lubridate)
library(GIGrvg)
library(lemon)
source("!qrdhs.R")

set.seed(1)
# Number of cores used in parallel MCMC process
cpu = 10

#load data
library(readxl)
data_raw <- read_excel("CPI_20231030.xlsx", 
                       col_types = c("date", "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric"))

#data processing
trnct <- F
if(trnct){
  start = which(as.Date(data_raw$Dates)=="2014-10-01")
  end = which(as.Date(data_raw$Dates)=="2023-09-01")
  # bandwidth = 150
  data <- data_raw[start:end,]
}else{
  data <- data_raw[]
}
rm(data_raw)
# time plot of data
for (i in 2:5) {
  plot(ts(data[,i], start = c(2000,3,1), frequency = 12),main=names(data)[i])
}

#predict and validation setting
simu <- T
# The interested quantiles
grid.p <- c(seq(0.05,0.95,by=0.05))

colNum <- 3
if(simu){
  fhorz <- 10
  Xout <- as.matrix(cbind(intercept=rep(1,fhorz), 
                          data[(nrow(data)-fhorz):(nrow(data)-1),colNum]))
  X <- as.matrix(cbind(intercept=rep(1,nrow(data)-1-nrow(Xout)), scale(data[1:(nrow(data)-fhorz-1),colNum])))
  Y <- matrix(400*diff(log(unlist(data$CPI)))[1:(nrow(data)-fhorz-1)],ncol=1)
  Yreal <- matrix(400*diff(log(unlist(data$CPI)))[(nrow(data)-fhorz):(nrow(data)-1)],ncol=1)
  T <- nrow(Y) + nrow(Yreal)
}else{
  Y <- matrix(400*diff(log(unlist(data$CPI))),ncol=1)
  X <- as.matrix(cbind(intercept=rep(1,nrow(data)-1), data[-nrow(data),colNum]))
  T <- nrow(Y)
}


# grid.p <- c(0.05,0.1,0.25,0.5,0.9,0.95)

# MCMC settings
nburn <- 3000
nsave <- 9000
thinfac <- 3

prior <- "dhs"
sl.cn <- "PRC"
quant_ls <- list()


# ---------------------------------------------------------------------------------------------------------
# insample
# Y <- matrix(400*diff(log(unlist(data$CPI)))[-(nrow(data)-1)],ncol=1)
# plot(ts(Y, start = c(2000,3,1), frequency = 12))

# causal discovery
library(causalXtreme)
# cd <- ease(data_raw[,-2:-1]);cd
# X <- as.matrix(cbind(intercept=rep(1,nrow(data_raw)-1), 400*scale(unlist(data_raw[-1,cd[1:3]+2]))))
# 
# cd <- ease(data_raw[,-1]);cd
# X <- as.matrix(cbind(intercept=rep(1,nrow(data_raw)-1), data_raw[-1,cd[1]+1]))



# Xout <- as.matrix(cbind(intercept=rep(1,1), data[nrow(data)-1,2]))
# X <- as.matrix(cbind(intercept=rep(1,nrow(data)-1-nrow(Xout)), data[c(-1,-(nrow(data)-1)),2]))



# X <- as.matrix(cbind(intercept=rep(1,nrow(data)-2),
#                      apply(data[,3], 2, function(x){diff(diff(diff(x)))})))


# cd <- ease(data[,1:5]);names(data)[cd+1]
# X <- as.matrix(cbind(intercept=rep(1,nrow(data)-1), data[-1,cd[1]+1]))
# X <- as.matrix(cbind(intercept=rep(1,nrow(data)-1), apply(data[,cd[1:2]+1], 2, function(x){diff(log(x))})))
# X <- as.matrix(cbind(intercept=rep(1,nrow(data_raw)-1), 400*diff(unlist(data_raw[,cd[1]+2]))))
# insample
# Y <- matrix(unlist(data_raw[,2]),ncol=1)
# X <- as.matrix(cbind(intercept=rep(1,nrow(data_raw)), data_raw[,-2:-1]))

# X <- cbind(intercept=rep(1,nrow(data_raw)-1), 400*apply(data_raw[,-2:-1], 2, function(x){log(x)[-1]}))
# X <- cbind(intercept=rep(1,nrow(data_raw)-1), apply(data_raw[,-2:-1], 2, function(x){diff(log(x))}))
# X <- cbind(intercept=rep(1,nrow(data_raw)),apply(data_raw[,-2:-1],2,function(x){x-mean(x)}))
# X <- cbind(intercept=rep(1,nrow(data_raw)), data_raw[,-2:-1])
# T <- nrow(Y)
dates <- as.character(data$Dates)[-1]
# dates <- as.character(format(date_decimal(as.numeric(time(data_raw[[sl.cn]]))),"%Y-%m-%d")[-1])

# UC-QR
quant_store <- array(NA,dim=c(2,T,length(grid.p)))
dimnames(quant_store) <- list(c("UC-QR","UC-QR-SV"),dates,paste0("p",grid.p*100))

for(sv in c(FALSE)){
  message("Estimating: ",sl.cn,", sv=",sv,".")
  Yt <- Y
  Xt <- X
  est <- tvpqr.grid(Y=Yt,X=Xt,Xout=Xout,p=grid.p,cpu=cpu,tvp=prior,sv=sv,fhorz=nrow(Xout),
                    nburn=nburn,nsave=nsave,thinfac=thinfac,out="mcmc",pred.dsp=T,pred.sv=T)
  bt_store <- est$bt
  sig2_store <- est$sig2
  fcst <- est$fcst
  fcstsig2 <- est$fcstsig2
  # if(sv){
  #   quant_store["UC-QR-SV",,] <- t(apply(apply(bt_store,c(2,3),remove_outliers),c(2,3),mean,na.rm=TRUE))
  # }else{
  #   quant_store["UC-QR",,] <- t(apply(apply(bt_store,c(2,3),remove_outliers),c(2,3),mean,na.rm=TRUE))
  # }
  if(sv){
    # quant_store["UC-QR-SV",,] <- t(apply(apply(bt_store,c(2,3),remove_outliers),c(2,3),mean,na.rm=TRUE))
    quant_store["UC-QR-SV",,] <- apply(apply(apply(bt_store,c(2,3,4),remove_outliers),
                                             c(2,3,4),mean,na.rm=TRUE),1,function(x){rowSums(x*X)})
  }else{
    quant_store["UC-QR",,] <- apply(apply(apply(bt_store,c(2,3,4),remove_outliers),
                                          c(2,3,4),mean,na.rm=TRUE),1,function(x){rowSums(x*X)})
  }
  rm(est, bt_store, sig2_store)
}
quant_ls[[sl.cn]] <- quant_store

# ------------------------------------
# output
library(reshape2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

shade.col <- "#1974D2"
r80col <- "#2c7bb6"
r90col <- "#d7191c"
shade.col.ser <- c('#f7fbff','#deebf7','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#08519c','#08306b')
alpha <- 0.9
sl.mods <- c("UCQR-TIS","UCQR-TVS")

# Y <- matrix(data_raw[[sl.cn]],ncol=1)[-1,]
quant_tmp <- melt(quant_ls[[sl.cn]]) %>% rename("Model"="Var1","Time"="Var2","Quantile"="Var3") %>% pivot_wider(names_from=Quantile,values_from=value)
levels(quant_tmp$Model) <- c("UCQR-TIS","UCQR-TVS")

quant_long <- melt(quant_tmp) %>% rename("Quantile"="variable")
quant_tmp <- cbind(quant_tmp,Y)

pp <- quant_tmp %>%
  subset(Model %in% sl.mods) %>%
  ggplot(aes(x=as.Date(Time))) +
  
  geom_ribbon(aes(ymin=p5,ymax=p95),alpha=alpha,fill=shade.col.ser[1]) + 
  geom_ribbon(aes(ymin=p10,ymax=p90),alpha=alpha,fill=shade.col.ser[2]) +
  geom_ribbon(aes(ymin=p15,ymax=p85),alpha=alpha,fill=shade.col.ser[3]) +
  geom_ribbon(aes(ymin=p20,ymax=p80),alpha=alpha,fill=shade.col.ser[4]) +
  geom_ribbon(aes(ymin=p25,ymax=p75),alpha=alpha,fill=shade.col.ser[5]) +
  geom_ribbon(aes(ymin=p30,ymax=p70),alpha=alpha,fill=shade.col.ser[6]) +
  geom_ribbon(aes(ymin=p35,ymax=p65),alpha=alpha,fill=shade.col.ser[7]) +
  geom_ribbon(aes(ymin=p40,ymax=p60),alpha=alpha,fill=shade.col.ser[8]) +
  geom_ribbon(aes(ymin=p45,ymax=p55),alpha=alpha,fill=shade.col.ser[9]) +
  geom_line(aes(y=value,group=Quantile),color="black",size=0.7,data=subset(subset(quant_long,Model %in% sl.mods),Quantile %in% c("p5","p95"))) +
  
  geom_hline(yintercept=0,size=0.5,color="red") +
  geom_line(aes(y=p50),color="white",size=0.7) +
  facet_rep_wrap(.~Model,scales="fixed",repeat.tick.labels = 'all') + xlab("") + ylab(expression(pi[t])) +
  coord_cartesian(expand=FALSE) +
  theme_cowplot() + 
  theme(axis.line.x = element_blank(), 
        panel.grid.major=element_line(color="grey80",linetype="solid",size=0.5), panel.grid.minor=element_line(color="grey80",linetype="dashed",size=0.5),
        strip.background = element_blank(), strip.text = element_text(face="bold",size=17),
        panel.spacing.x = unit(1,"cm"),
        axis.text = element_text(size=16), axis.title = element_text(size=17),
        plot.margin = unit(c(0,1,0,0), "cm"))

print(pp)