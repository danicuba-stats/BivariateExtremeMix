################################################################################
##### Bivariate mixed coregionalization model - Case study #####################
##### Date: 30th of November, 2022
##### Author: Daniela Cuba (m.cuba.1@research.gla.ac.uk)
################################################################################


### Load libraries ----
library("INLA")
library("parallel")
library("dplyr")
library("evd")
library("gtools")

### Load functions ----
source("./bivmix_functions.R")

### Set up parameters ----
n.cores<-1
elem.name<-c("Cr-Pb")
elem_pairs<-strsplit(elem.name,"-")[[1]]
p.grid<-c(0.98,0.95)

# Load in data ----
data<-read.csv("./my_data.csv") # Data can be requested from the British Geological Survey
pred_data<-read.csv("./my_predictions.csv") # Data can be requested from the British Geological Survey
# Both of these datasets include covariate information at every observation location. 

# Define functions ---- 
fit.model<-function(ps,d=data,e.name=elem.name,pred) {
  p.grid<-ps
  
  # Obtain initial values
  w_init<-list()
  init_idx<-list()
  pars_init<-list()
  
  for(m in 1:2){
    all.norm<-vector("numeric",
                     nrow(d$data$var1))
    all.gpd<-vector("numeric",
                    nrow(d$data$var2))
    idx.groups<-vector("numeric",
                       nrow(d$data$var1))
    # Normal distribution - body
    norm.mean<-d$data[[m]]$y.data %>% mean
    norm.sd<-d$data[[m]]$y.data %>% sd
    # Extreme distribution - tail
    gpd.pars<-evd::fpot(d$data[[m]]$y.data,
                        threshold=quantile(d$data[[m]]$y.data,
                                           p.grid[m]))#data$var1$y.data
    
    # Categorise
    for(j in 1:nrow(d$data[[m]])){
      pr.norm<-dnorm(d$data[[m]]$y.data[j],
                     mean=norm.mean,
                     sd=norm.sd)
      pr.gpd<-dgpd(as.numeric(d$data[[m]]$y.data[j]),#-quantile(d$data[[m]]$y.data,p.grid[i,m])),
                   loc=quantile(d$data[[m]]$y.data,p.grid[m]),
                   scale=gpd.pars$param[1],
                   shape = gpd.pars$param[2])
      prop<-min(1,pr.norm/pr.gpd)
      idx.groups[j]<-ifelse(runif(1)>prop,2,1)
    }
    
    w_init[[m]]<-table(idx.groups)/length(idx.groups)
    pars_init[[m]]<-list(norm=c(norm.mean,
                                norm.sd),
                         gpd=c(quantile(d$data[[m]]$y.data,p.grid[m]),
                               gpd.pars$param[1],
                               gpd.pars$param[2]))
    init_idx[[m]]<-idx.groups
  }

  init<-list(w = w_init,# mu = means_init, taus=tau_init,
             pars=pars_init,
             init.idx = init_idx)

  # Fit model
  init_pred<-fit.inlaQPRED(data=d$data,
                           data_pred=pred,
                           idx.groups=init$init.idx,
                           covs=d$covariate,
                           p=ps,
                           final=T,
                           pars=init$pars)
  
  inla.q.cs<-list(pred=init_pred,
                  init.vals=init) # Compile information
  
  return(inla.q.cs)
}



## Case Study - Cr-Pb
elems<-c("Cr","Pb")
d<-data[,c("X_COORD","Y_COORD",
           paste0(elems,"_XRF"),
           "elevation","slope","aspect", "plan.curve", "profile.curve",
           "twi","mrvbf","mrrtf",
           "Distance.To.Road","Distance.To.MAB")]
d$var1<-d[,paste0(elems[1],"_XRF")] %>% as.numeric %>% log
d$var2<-d[,paste0(elems[2],"_XRF")] %>% as.numeric %>% log
d<-d[complete.cases(d),]
log1<-d[,paste0(elems[1],"_XRF")] %>% as.numeric %>% log
log2<-d[,paste0(elems[2],"_XRF")] %>% as.numeric %>% log

# Standardise Covariates
standard.covs<-apply(d[,c("elevation","slope","aspect", "plan.curve", 
                          "profile.curve", "twi","mrvbf","mrrtf","Distance.To.Road",
                          "Distance.To.MAB")],2,stand.covs) %>% as.data.frame 
# "Distance.To.Road","Distance.To.MAB")],2,stand.covs) %>% as.data.frame 
dats<-list(data=list(var1=data.frame(X=d$X_COORD,
                                     Y=d$Y_COORD,
                                     y.data=log1),
                     var2=data.frame(X=d$X_COORD,
                                     Y=d$Y_COORD,
                                     y.data=log2)),
           covariate=standard.covs)
# Standarise Covaraites for predction locations
pred_data<-pred_data[,c("X","Y","elevation","slope","aspect", "plan.curve", 
                        "profile.curve", "twi","mrvbf","mrrtf","Distance.To.Road",
                        "Distance.To.MAB")]
pred_data <- pred_data[complete.cases(pred_data),]
pred.covs<-apply(pred_data[,c("elevation","slope","aspect", "plan.curve", 
                              "profile.curve", "twi","mrvbf","mrrtf","Distance.To.Road",
                              "Distance.To.MAB")],2,stand.covs) %>% as.data.frame 
pred_data<-cbind(pred_data[,c("X","Y")],
                 pred.covs)

# Fit model
final<- fit.model(ps=p.grid,
                  d=dats,
                  pred=pred_data)

