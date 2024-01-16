
# # # Fit inla function ( no predictions )
# # data: data[loc_1,loc_2,data]
# # idx.groups: list length 2 of indeces
# # cov: list length 2 of covariates
# # p: probabilities (p1,p2)
# # pars: hyperparameters


# define fit.inlaQ function
fit.inlaQ<-function(data,idx.groups,covs,final=F,p,pars){
  
  # Membership (body or tail)
  idx.groups1<-idx.groups[[1]]
  idx.groups2<-idx.groups[[2]]
  
  ##  Var 1
  # y1_body
  y_b1<-data[[1]]
  y_b1[idx.groups1==2,"y.data"]<-NA
  cov_b1<-data.frame(covs)
  cov_b1[idx.groups1==2,]<-NA
  
  # y1_tail
  y_t1<-data[[1]]
  qs.t1<-pgpd(y_t1$y.data[which(idx.groups1==2)],
              loc=pars[[1]]$gpd[1],
              scale=pars[[1]]$gpd[2],
              shape=pars[[1]]$gpd[3])
  qs.t1<-qnorm(qs.t1)
  y_t1$y.data<-NA#y_t1$y.data-quantile(data[[1]]$y.data,p[1])
  y_t1[idx.groups1==2,"y.data"]<-qs.t1
  cov_t1<-data.frame(covs)
  cov_t1[idx.groups1==1,]<-NA
  
  # Var 2
  # y2_body
  y_b2<-data[[2]]
  y_b2[idx.groups2==2,"y.data"]<-NA
  cov_b2<-data.frame(covs)
  cov_b2[idx.groups2==2,]<-NA
  
  # y2_tail
  y_t2<-data[[2]]
  qs.t2<-pgpd(y_t2$y.data[which(idx.groups2==2)],
              loc=pars[[2]]$gpd[1],
              scale=pars[[2]]$gpd[2],
              shape=pars[[2]]$gpd[3])
  qs.t2<-qnorm(qs.t2)
  y_t2$y.data<-NA
  y_t2[idx.groups2==2,"y.data"]<-qs.t2
  cov_t2<-data.frame(covs)
  cov_t2[idx.groups2==1,]<-NA

  ## Mesh
  mesh_all<-inla.mesh.2d(loc=data[[1]][,c("X","Y")],
                         loc.domain = matrix(data=c(min(data[[1]][,c("X")]),
                                                    min(data[[1]][,c("Y")]),
                                                    min(data[[1]][,c("X")]),
                                                    max(data[[1]][,c("Y")]),
                                                    max(data[[1]][,c("X")]),
                                                    min(data[[1]][,c("Y")]),
                                                    max(data[[1]][,c("X")]),
                                                    max(data[[1]][,c("Y")])),
                                             ncol=2,byrow=T),
                         max.edge=c(100000,100000),
                         offset=1)
  
  ## SPDE
  spde_b1<-inla.spde2.pcmatern(mesh=mesh_all,
                               alpha=2,
                               prior.range=c(20000,0.05),
                               prior.sigma=c(1,0.01))
  spde_t1<-inla.spde2.pcmatern(mesh=mesh_all,
                               alpha=2,
                               prior.range=c(30000,0.05),
                               prior.sigma=c(1,0.95))
  
  spde_b2<-inla.spde2.pcmatern(mesh=mesh_all,
                               alpha=2,
                               prior.range=c(50000,0.5),
                               prior.sigma=c(1,0.95))
  spde_t2<-inla.spde2.pcmatern(mesh=mesh_all,
                               alpha=2,
                               prior.range=c(30000,0.05),
                               prior.sigma=c(1,0.95))
  spde_all<-inla.spde2.pcmatern(mesh=mesh_all,
                                alpha=2,
                                prior.range=c(50000,0.05),
                                prior.sigma=c(1,0.95))
  
  
  # Projection Matrices
  # Var 1
  A_b1<- inla.spde.make.A(mesh=mesh_all,
                          loc=as.matrix(y_b1[,c("X","Y")]))
  A_t1<- inla.spde.make.A(mesh=mesh_all,
                          loc=as.matrix(y_t1[,c("X","Y")])) # same mesh using all coordinates and a single mesh
  # Var 2
  A_b2<- inla.spde.make.A(mesh=mesh_all,
                          loc=as.matrix(y_b2[,c("X","Y")]))
  A_t2<- inla.spde.make.A(mesh=mesh_all,
                          loc=as.matrix(y_t2[,c("X","Y")]))
  A_all<- inla.spde.make.A(mesh=mesh_all,
                           loc=as.matrix(data[[1]][,c("X","Y")]))
  
  ## Stacks
  # Var 1
  # Stack 1 - bulk
  stack.covb1<-lapply(1:ncol(cov_b1),function(x) cov_b1[,x])
  names(stack.covb1)<-paste0(names(cov_b1),".b1")
  stack_b1<-inla.stack(data=list(y=cbind(as.vector(y_b1$y.data),NA,NA,NA)),
                       A=list(A_b1,1),
                       # A=list(1),
                       effects = list(list(z1=1:spde_b1$n.spde),
                                      c(intercept1=1,stack.covb1)), tag="y_b1.data")
  
  
  # Stack 2 - tail
  stack.covt1<-lapply(1:ncol(cov_t1),function(x) cov_t1[,x])
  names(stack.covt1)<-paste0(names(cov_t1),".t1")
  stack_t1<-inla.stack(data=list(y=cbind(NA,as.vector(y_t1$y.data),NA,NA)),
                       A=list(A_t1,1),
                       effects = list(list(z2=1:spde_t1$n.spde),
                                      c(stack.covt1)), tag="y_t1.data")
  # Var 2
  # Stack 3 - bulk
  stack.covb2<-lapply(1:ncol(cov_b2),function(x) cov_b2[,x])
  names(stack.covb2)<-paste0(names(cov_b2),".b2")
  stack_b2<-inla.stack(data=list(y=cbind(NA,NA,as.vector(y_b2$y.data),NA)),
                       A=list(A_b2,1),
                       effects = list(list(z3=1:spde_b2$n.spde),
                                      c(intercept3=1,stack.covb2)), tag="y_b2.data")
  
  
  # Stack 4 - tail
  stack.covt2<-lapply(1:ncol(cov_t2),function(x) cov_t2[,x])
  names(stack.covt2)<-paste0(names(cov_t2),".t2")
  stack_t2<-inla.stack(data=list(y=cbind(NA,NA,NA,as.vector(y_t2$y.data))),
                       A=list(A_t2,1),
                       effects = list(list(z24=1:spde_t1$n.spde,
                                           z4=1:spde_t2$n.spde),
                                      c(stack.covt2)), tag="y_t2.data")
  
  # Stack all - all data
  all.covs<-lapply(1:ncol(covs),function(x) covs[,x])
  names(all.covs)<-names(covs)
  stack_all<-inla.stack(data=list(y=cbind(rep(NA,nrow(data[[1]][,c("X","Y")])),
                                          NA,NA,NA)),
                        A=list(A_all,1),
                        effects=list(list(z1=1:spde_all$n.spde,
                                          z2=1:spde_all$n.spde,
                                          z3=1:spde_all$n.spde,
                                          z24=1:spde_all$n.spde,
                                          z4=1:spde_all$n.spde),
                                     c(intercept1=1,
                                       intercept3=1,
                                       all.covs)), tag="y_all.data") # do I need this?
  
  # Merge stacks
  monster.stack<-inla.stack(stack_b1,
                            stack_t1,
                            stack_b2,
                            stack_t2,
                            stack_all)
  ## Formula
  form<-y ~ -1 + intercept1 + intercept3 +
    f(z1,model=spde_b1) + f(z3,model=spde_b2) +
    f(z2,model=spde_t1)+ f(z4,model=spde_t2) + f(z24,copy='z2', fixed=F) +
    elevation.b1 + slope.b1 + aspect.b1 + plan.curve.b1 + profile.curve.b1 +
    twi.b1 + mrvbf.b1 + mrrtf.b1 + Distance.To.Road.b1 + Distance.To.MAB.b1 +
    elevation.t1 + slope.t1 + aspect.t1 + plan.curve.t1 + profile.curve.t1 +
    twi.t1 + mrvbf.t1 + mrrtf.t1 + Distance.To.Road.t1 + Distance.To.MAB.t1 +
    elevation.b2 + slope.b2 + aspect.b2 + plan.curve.b2 + profile.curve.b2 +
    twi.b2 + mrvbf.b2 + mrrtf.b2 + Distance.To.Road.b2 + Distance.To.MAB.b2 +
    elevation.t2 + slope.t2 + aspect.t2 + plan.curve.t2 + profile.curve.t2 +
    twi.t2 + mrvbf.t2 + mrrtf.t2 + Distance.To.Road.t2 + Distance.To.MAB.t2
  
  
  ## Fit using  INLA
  # Define riors
  hyper_norm1 <- list(hyper = list(prec = list(prior = 'pc.prec',
                                               param = c(sd(data$var1$y.data)*2, 0.995))))
  hyper_norm2 <- list(hyper = list(prec = list(prior = 'pc.prec',
                                               param = c((sd(y_t1$y.data,na.rm=T)*0.2), 0.999))))
  hyper_norm3 <- list(hyper = list(prec = list(prior = 'pc.prec',
                                               param = c(sd(data$var2$y.data)*2, 0.995))))
  hyper_norm4 <- list(hyper = list(prec = list(prior = 'pc.prec',
                                               param = c((sd(y_t2$y.data,na.rm = T)*2), 0.999))))
  
  
  
  # Fit
  inla.bi<-inla(formula = form,
                data=inla.stack.data(monster.stack),
                family=c("gaussian","gaussian",
                         "gaussian","gaussian"),
                control.family=list(hyper_norm1, hyper_norm2,
                                    hyper_norm3, hyper_norm4),
                control.predictor=list(A=inla.stack.A(monster.stack), compute=TRUE),
                control.compute=list(cpo=TRUE,dic=TRUE,config=TRUE))
  
  dists_list<-c(inla.bi$marginals.fixed,
                inla.bi$marginals.hyperpar)
  names(dists_list)<-c(names(inla.bi$marginals.fixed),
                       names(inla.bi$marginals.hyperpar))
  
  # Save output
  
  if(final){
    return(list(model=inla.bi,
                mlik = inla.bi$mlik[[1]],
                dists = dists_list,
                indeces=list(b1=inla.stack.index(monster.stack,"y_b1.data")$data,
                             t1=inla.stack.index(monster.stack,"y_t1.data")$data,
                             b2=inla.stack.index(monster.stack,"y_b2.data")$data,
                             t2=inla.stack.index(monster.stack,"y_t2.data")$data)))
    
  }else{
    return(list(model=inla.bi,
                mlik = inla.bi$mlik[[1]],
                dists = dists_list))
  }
}

# # # Fit inla function
# # data: data[loc_1,loc_2,data]
# # idx.groups: list length 2 of indeces
# # cov: list length 2 of covariates
# # p: probabilities (p1,p2)
# # pars: hyperparameters
# # data_pred: prediction locations and covariates


fit.inlaQPRED<-function(data,idx.groups,covs,final=F,p,pars,data_pred){
  
  # Membership (body or tail)
  idx.groups1<-idx.groups[[1]]
  idx.groups2<-idx.groups[[2]]
  
  # Initialise prediction objects
  pred.preds<-data.frame(X=data_pred$X,
                         Y=data_pred$Y,
                         y.data=NA)
  pred.covs<-pred_data[,names(covs)]
  
  ##  Var 1
  # y1_body
  y_b1<-data[[1]]
  y_b1[idx.groups1==2,"y.data"]<-NA
  y_b1<-rbind(y_b1,pred.preds)# Add prediction Values
  cov_b1<-data.frame(covs)
  cov_b1[idx.groups1==2,]<-NA
  cov_b1<-rbind(cov_b1,pred.covs)
  
  # y1_tail
  y_t1<-data[[1]]
  qs.t1<-pgpd(y_t1$y.data[which(idx.groups1==2)],
              loc=pars[[1]]$gpd[1],
              scale=pars[[1]]$gpd[2],
              shape=pars[[1]]$gpd[3])
  qs.t1<-qnorm(qs.t1)
  y_t1$y.data<-NA
  y_t1[idx.groups1==2,"y.data"]<-qs.t1
  y_t1<-rbind(y_t1,pred.preds)
  cov_t1<-data.frame(covs)
  cov_t1[idx.groups1==1,]<-NA
  cov_t1<-rbind(cov_t1,pred.covs)
  
  ## Var 2
  # y2_body
  y_b2<-data[[2]]
  y_b2[idx.groups2==2,"y.data"]<-NA
  y_b2<-rbind(y_b2,pred.preds)
  cov_b2<-data.frame(covs)
  cov_b2[idx.groups2==2,]<-NA
  cov_b2<-rbind(cov_b2,pred.covs)
  
  # y2_tail
  y_t2<-data[[2]]
  qs.t2<-pgpd(y_t2$y.data[which(idx.groups2==2)],
              loc=pars[[2]]$gpd[1],
              scale=pars[[2]]$gpd[2],
              shape=pars[[2]]$gpd[3])
  qs.t2<-qnorm(qs.t2)
  y_t2$y.data<-NA
  y_t2[idx.groups2==2,"y.data"]<-qs.t2
  y_t2<-rbind(y_t2,pred.preds)
  cov_t2<-data.frame(covs)
  cov_t2[idx.groups2==1,]<-NA
  cov_t2<-rbind(cov_t2,pred.covs)
  
  # All covariates 
  all_pred_covs<-rbind(covs,pred.covs)
  
  ## Mesh
  mesh_all<-inla.mesh.2d(loc=data[[1]][,c("X","Y")],
                         loc.domain = matrix(data=c(min(y_b1[,c("X")]),
                                                    min(y_b1[,c("Y")]),
                                                    min(y_b1[,c("X")]),
                                                    max(y_b1[,c("Y")]),
                                                    max(y_b1[,c("X")]),
                                                    min(y_b1[,c("Y")]),
                                                    max(y_b1[,c("X")]),
                                                    max(y_b1[,c("Y")])),
                                             ncol=2,byrow=T),
                         max.edge=c(50000,50000))
  
  ## SPDE
  spde_b1<-inla.spde2.pcmatern(mesh=mesh_all,
                               alpha=2,
                               prior.range=c(20000,0.05),
                               prior.sigma=c(1,0.01))
  spde_t1<-inla.spde2.pcmatern(mesh=mesh_all,
                               alpha=2,
                               prior.range=c(30000,0.05),
                               prior.sigma=c(1,0.95))
  
  spde_b2<-inla.spde2.pcmatern(mesh=mesh_all,
                               alpha=2,
                               prior.range=c(50000,0.5),
                               prior.sigma=c(1,0.95))
  spde_t2<-inla.spde2.pcmatern(mesh=mesh_all,
                               alpha=2,
                               prior.range=c(30000,0.05),
                               prior.sigma=c(1,0.95))
  spde_all<-inla.spde2.pcmatern(mesh=mesh_all,
                                alpha=2,
                                prior.range=c(50000,0.05),
                                prior.sigma=c(1,0.95))
  
  
  # Projection Matrices
  # Var 1
  A_b1<- inla.spde.make.A(mesh=mesh_all,
                          loc=as.matrix(y_b1[,c("X","Y")]))
  A_t1<- inla.spde.make.A(mesh=mesh_all,
                          loc=as.matrix(y_t1[,c("X","Y")])) # same mesh using all coordinates and a single mesh
  # Var 2
  A_b2<- inla.spde.make.A(mesh=mesh_all,
                          loc=as.matrix(y_b2[,c("X","Y")]))
  A_t2<- inla.spde.make.A(mesh=mesh_all,
                          loc=as.matrix(y_t2[,c("X","Y")]))
  A_all<- inla.spde.make.A(mesh=mesh_all,
                           loc=as.matrix(y_t2[,c("X","Y")]))
  
  ## Stacks
  # Var 1
  # Stack 1 - body
  stack.covb1<-lapply(1:ncol(cov_b1),function(x) cov_b1[,x])
  names(stack.covb1)<-names(cov_b1)
  stack_b1<-inla.stack(data=list(y=cbind(as.vector(y_b1$y.data),NA,NA,NA)),
                       A=list(A_b1,1),
                       # A=list(1),
                       effects = list(list(z1=1:spde_b1$n.spde),
                                      c(intercept1=1,stack.covb1)), tag="y_b1.data")
  
  
  # Stack 2 - tail
  stack.covt1<-lapply(1:ncol(cov_t1),function(x) cov_t1[,x])
  names(stack.covt1)<-names(cov_t1)
  stack_t1<-inla.stack(data=list(y=cbind(NA,as.vector(y_t1$y.data),NA,NA)),
                       A=list(A_t1,1),
                       effects = list(list(z2=1:spde_t1$n.spde),
                                      # z42=1:spde_t2$n.spde),
                                      c(intercept2=1,
                                        stack.covt1)), tag="y_t1.data")
  # Var 2
  # Stack 3 - body
  stack.covb2<-lapply(1:ncol(cov_b2),function(x) cov_b2[,x])
  names(stack.covb2)<-names(cov_b2)
  stack_b2<-inla.stack(data=list(y=cbind(NA,NA,as.vector(y_b2$y.data),NA)),
                       A=list(A_b2,1),
                       # A=list(1),
                       effects = list(list(z3=1:spde_b2$n.spde),
                                      c(intercept3=1,stack.covb2)), tag="y_b2.data")
  
  
  # Stack 4 - tail
  stack.covt2<-lapply(1:ncol(cov_t2),function(x) cov_t2[,x])
  names(stack.covt2)<-names(cov_t2)
  stack_t2<-inla.stack(data=list(y=cbind(NA,NA,NA,as.vector(y_t2$y.data))),
                       A=list(A_t2,1),
                       effects = list(list(z24=1:spde_t1$n.spde,
                                           z4=1:spde_t2$n.spde),
                                      c(intercept4=1,
                                        stack.covt2)), tag="y_t2.data")
  
  # Stack all - all data
  all.covs<-lapply(1:ncol(all_pred_covs),function(x) all_pred_covs[,x])
  names(all.covs)<-names(all_pred_covs)
  stack_all<-inla.stack(data=list(y=cbind(rep(NA,nrow(y_b1)),
                                          NA,NA,NA)),
                        A=list(A_all,1),
                        effects=list(list(z1=1:spde_all$n.spde,
                                          z2=1:spde_all$n.spde,
                                          z3=1:spde_all$n.spde,
                                          z24=1:spde_all$n.spde,
                                          z4=1:spde_all$n.spde),
                                     c(intercept1=1,
                                       intercept2=1,
                                       intercept3=1,
                                       intercept4=1,
                                       all.covs)), tag="y_all.data")
  
  # Merge stacks
  monster.stack<-inla.stack(stack_b1,
                            stack_t1,
                            stack_b2,
                            stack_t2,
                            stack_all)
  ## Formula
  form<-y ~ -1 + intercept1 + intercept2 + intercept3 + intercept4+
    f(z1,model=spde_b1) + f(z3,model=spde_b2) +
    f(z2,model=spde_t1)+ f(z4,model=spde_t2) + 
    f(z24,copy='z2', fixed=F) + elevation + slope +
    aspect + plan.curve + profile.curve + twi + mrvbf +
    mrrtf  + Distance.To.Road + Distance.To.MAB
  
  ## Fit using INLA
  # Define priors
  hyper_norm1 <- list(hyper = list(prec = list(prior = 'pc.prec',
                                               param = c(sd(data$var1$y.data)*2, 0.995))))
  hyper_norm2 <- list(hyper = list(prec = list(prior = 'pc.prec',
                                               param = c((sd(y_t1$y.data,na.rm=T)*2), 0.999))))
  hyper_norm3 <- list(hyper = list(prec = list(prior = 'pc.prec',
                                               param = c(sd(data$var2$y.data)*2, 0.995))))
  
  hyper_norm4 <- list(hyper = list(prec = list(prior = 'pc.prec',
                                               param = c((sd(y_t2$y.data,na.rm = T)*2), 0.999))))
  
  # Fit
  inla.bi<-inla(formula = form,
                data=inla.stack.data(monster.stack),
                family=c("gaussian","gaussian",
                         "gaussian","gaussian"),
                control.family=list(hyper_norm1, hyper_norm2,
                                    hyper_norm3, hyper_norm4),
                control.predictor=list(A=inla.stack.A(monster.stack), compute=TRUE),
                control.compute=list(config=TRUE, return.marginals.predictor=TRUE))
  
  # Save output
  if(final){
    return(list(model=inla.bi,
                mlik = inla.bi$mlik[[1]],
                indeces=list(b1=inla.stack.index(monster.stack,"y_b1.data")$data,
                             t1=inla.stack.index(monster.stack,"y_t1.data")$data,
                             b2=inla.stack.index(monster.stack,"y_b2.data")$data,
                             t2=inla.stack.index(monster.stack,"y_t2.data")$data)))
  }else{
    return(list(model=summary(inla.bi),
                mlik = inla.bi$mlik[[1]],
                indeces=list(b1=inla.stack.index(monster.stack,"y_b1.data")$data,
                             t1=inla.stack.index(monster.stack,"y_t1.data")$data,
                             b2=inla.stack.index(monster.stack,"y_b2.data")$data,
                             t2=inla.stack.index(monster.stack,"y_t2.data")$data)))
  }
}

# Standardise covariates
stand.covs<-function(x){
  return((x-mean(x))/sd(x))
}

