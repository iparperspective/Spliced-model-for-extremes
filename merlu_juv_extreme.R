
########################################
### final5 merlu juv extreme





require(INLA)
require(rgdal)
require(ggplot2)

server="~/Dropbox/servidor/atlas/"
hp="C:/Users/ip30/Dropbox/Marie Curie/projects/MC objectives/atlas/"
link=hp

if(link==server){INLA:::inla.dynload.workaround()}

save.link = paste0(link,"Merlu_extreme")

load(file=paste0(link,"data/merluccius/merlu_juv.Rdata"))
bathy<-readGDAL(paste0(link,"data/bathy.tif"))
bathy$band1[which(bathy$band1==0)]<-NA

datos=merlu_juv
datos$occur=ifelse(datos$weight==0,0,1)
datos$cpue=datos$weight/datos$INICIO.VIRADO
datos$numbers[is.na(datos$numbers)]=0
datos$npue=round(datos$numbers/datos$INICIO.VIRADO)
datos$year2=datos$year-min(datos$year)+1

### remove baleares
idx_baleares = which(datos$LONGITUD>0.7 & datos$LATITUD<39.5)
datos = datos[-idx_baleares,]
plot(datos$LONGITUD ,datos$LATITUD)
head(datos)


#### utm
points=SpatialPoints(coords=cbind(datos$LONGITUD,datos$LATITUD),
                     proj4string=CRS("+proj=longlat"))
utm=spTransform(points, CRS=CRS(proj4string(bathy)))
datos$x=coordinates(utm)[,1]
datos$y=coordinates(utm)[,2]
urtebat=which(datos$year==2001)
datos$bathy = over(utm,bathy)$band1

##### mesh ####
spat.scale=100000
d<-matrix(c(datos$x/spat.scale,datos$y/spat.scale),ncol=2)
#boundary <- inla.nonconvex.hull(points=d,convex=-.06)
#max.edge <- c(max(dist(coordinates(d)))*.024,max(dist(coordinates(d)))*.1)
boundary <- inla.nonconvex.hull(points=d,convex=-.04)
#max.edge <- c(max(dist(coordinates(d)))*.015,max(dist(coordinates(d)))*.1)
max.edge <- c(max(dist(coordinates(d)))*.010,max(dist(coordinates(d)))*.1) ## FINE
mesh <- inla.mesh.2d( boundary=boundary,max.edge=max.edge)
plot(mesh,xlab=paste("num knots =",mesh$n))
points(d[urtebat,],pch=16,cex=.4,col="red")

save(mesh,file=paste0(save.link,"/final5/mesh_extreme_fine.Rdata"))
# save(mesh,file=paste0("C:/Users/ip30/Dropbox (Personal)/Marie Curie/projects/MC objectives/Merlu_extreme/final5/plotting/plot_data/mesh_extreme.Rdata"))
# save(d,file=paste0("C:/Users/ip30/Dropbox (Personal)/Marie Curie/projects/MC objectives/Merlu_extreme/final5/plotting/plot_data/d.Rdata"))
# save(datos,file=paste0("C:/Users/ip30/Dropbox (Personal)/Marie Curie/projects/MC objectives/Merlu_extreme/final5/plotting/plot_data/datos.Rdata"))
load(file=paste0(save.link,"/final5/mesh_extreme_fine.Rdata"))

####
mesh.points = SpatialPoints(coords=mesh$loc[,1:2]*spat.scale,
                            proj4string=CRS(proj4string(bathy)))

#### covariate RW2
bathy_model = c(datos$PROFUNDIDAD,over(mesh.points,bathy)$band1)
bathy_model[which(bathy_model>max(datos$PROFUNDIDAD))]=max(datos$PROFUNDIDAD) ## set maximum depth
bathy_model[which(bathy_model<min(datos$PROFUNDIDAD))]=min(datos$PROFUNDIDAD) ## set minimum depth
bathy_model=inla.group(bathy_model,n=30)

datos$bathy_model=bathy_model[1:length(datos[,1])]

### covariate 1D matern
mesh1d <- inla.mesh.1d(seq(min(datos$PROFUNDIDAD), max(datos$PROFUNDIDAD), by = 10)) 
save(mesh1d,file=paste0(save.link,"/final5/plotting/plot_data/mesh1d.Rdata"))
A1 <- inla.spde.make.A(mesh1d, datos$PROFUNDIDAD)
spde1 <- inla.spde2.pcmatern(mesh1d, prior.range=c(40, 0.9), prior.sigma=c(1, 0.2))
spde1_NA <- inla.spde2.pcmatern(mesh1d, prior.range=c(150, NA), prior.sigma=c(1, 0.2))
spde1.idx <- inla.spde.make.index("bathy", n.spde = spde1$n.spde)

## pred
pred_bathy = over(mesh.points,bathy)$band1
pred_bathy[which(pred_bathy>max(datos$PROFUNDIDAD))]=max(datos$PROFUNDIDAD)
pred_bathy_NA = pred_bathy 
pred_bathy_NA[which(is.na(pred_bathy_NA))] = max(datos$PROFUNDIDAD)
A1_pred <- inla.spde.make.A(mesh1d,pred_bathy_NA)


bound=Polygon(boundary$loc*spat.scale)
poly=SpatialPolygons(list(Polygons(list(bound),'0')))
ea=point.in.polygon(coordinates(mesh.points)[,1], coordinates(mesh.points)[,2], 
                    coordinates(bound)[,1], coordinates(bound)[,2], mode.checked=FALSE)
idx.mesh=which(ea%in%c(0))
plot(mesh$loc[-idx.mesh,1:2]*spat.scale)

#### check bathy
datos$Pr=over(utm,bathy)$band1
idx=which(datos$Pr>820)

# 1) variables selection
prior_occur= list(theta = list(prior="pc.prec", param=c(.5,0.01)))
prior_npue= list(theta = list(prior="pc.prec", param=c(sd(log(datos$npue)),0.01)))
spde_npue <- inla.spde2.pcmatern(mesh, prior.range=c(1.5, 0.15), prior.sigma=c(1, 0.2))




#############################
####### Find extremes ###########
#############################
method="laplace"

spat.index <- inla.spde.make.A(mesh, loc=d)
spat.pred <- inla.spde.make.A(mesh)

stack_stage1<-inla.stack(data=list(y=datos$numbers),
                        A=list(spat.index, 1,A1),
                        effects=list(spat=1:mesh$n,
                                     list(b0=1,
                                          prof=datos$PROFUNDIDAD,
                                          year=datos$year2,
                                          effort=datos$INICIO.VIRADO),
                                     spde1.idx),
                        tag='est')

pred_numbers_2000<-inla.stack(data=list(y=rep(NA,mesh$n)),
                                  A=list(spat.pred, 1,A1_pred),
                                  effects=list(spat.con=1:mesh$n,
                                               list(b0=1,
                                                    prof=pred_bathy,
                                                    year=1),
                                               spde1.idx),
                                  tag='pred.numbers.2000')

pred_numbers_2016<-inla.stack(data=list(y=rep(NA,mesh$n)),
                              A=list(spat.pred, 1,A1_pred),
                              effects=list(spat.con=1:mesh$n,
                                           list(b0=1,
                                                prof=pred_bathy,
                                                year=17),
                                           spde1.idx),
                              tag='pred.numbers.2016')

stack_both_ident=inla.stack(stack_stage1,pred_numbers_2000,pred_numbers_2016)

form_indetify<-y~-1+b0+
  f(bathy,model=spde1_NA)+
  f(year,model="rw2")+
  f(spat,model=spde_npue) + 
  offset(log(effort))

stage1_fit <- inla(form_indetify,
                          family=c("zeroinflatednbinomial2"),
                          #control.fixed=list(expand.factor.strategy="inla"),
                          data=inla.stack.data(stack_stage1), control.compute=list(config=T,waic=F),
                          control.predictor=list(A=inla.stack.A(stack_stage1), compute=TRUE,link=1),
                          verbose=T,num.threads=2,control.inla = list(strategy = method))

# ident_numbers_pred <- inla(form_indetify,
#                           family=c("nbinomial"),
#                           #control.fixed=list(expand.factor.strategy="inla"),
#                           data=inla.stack.data(stack_both_ident), control.compute=list(config=T,waic=T),
#                           control.predictor=list(A=inla.stack.A(stack_both_ident), compute=TRUE,link=1),
#                           verbose=T,num.threads=2,control.inla = list(strategy = method))

stage1_fit$summary.hyperpar
stage1_fit$waic$waic
stage1_fit2$waic$waic

fitted = inla.stack.index(stack_stage1,"est")$data
means=stage1_fit$summary.fitted.values[fitted,1]
plot(means,datos$numbers)
abline(0,1)

save(stage1_fit,file=paste0(save.link,"/final5/stage1/stage1_fit.Rdata"))
save(stack_stage1,file=paste0(save.link,"/final5/stage1/stack_stage1.Rdata"))

# save(ident_numbers_pred,file=paste0(save.link,"/final5/stage1/ident_numbers_pred.Rdata"))
# save(stack_both_ident,file=paste0(save.link,"/final5/stage1/stack_both_ident.Rdata"))

plot(stage1_fit$marginals.fitted.values[[fitted[1]]])
hist(rnbinom(10000,size=stage1_fit$summary.hyperpar[1,1],
             mu=stage1_fit$summary.fitted.values[fitted[1],1]),breaks=100)

###########################################
############ extremes 85 ###############
##########################################

load(file=paste0(save.link,"/final5/stage1/stage1_fit.Rdata"))
load(file=paste0(save.link,"/final5/stage1/stack_stage1.Rdata"))
model=stage1_fit
fitted = inla.stack.index(stack_stage1,"est")$data


means=model$summary.fitted.values[fitted,1]
perc_85 = unlist(lapply(means,function(x){qnbinom(p=0.85,size=model$summary.hyperpar[1,1],mu=x)}))

extremes_presence85 = which(datos$numbers>=perc_85 & datos$npue > mean(datos$npue))
datos$extremes85=0
datos$extremes85[extremes_presence85]=1

A1_CRPS <- inla.spde.make.A(mesh1d,datos$PROFUNDIDAD[extremes_presence85])
spat.CRPS <- inla.spde.make.A(mesh, loc=d[extremes_presence85,])

stack_stage3<-inla.stack(data=list(y=cbind(datos$extremes85)),
                            A=list(spat.index, 1,A1),
                            effects=list(spat.bin=1:mesh$n,
                                         list(b0=1,
                                              prof=datos$PROFUNDIDAD,
                                              prof2=datos$PROFUNDIDAD^2,
                                              year=datos$year2), ## change 
                                         spde1.idx),
                            tag='est')

pred.extreme.occur.2000<-inla.stack(data=list(y=rep(NA,mesh$n)),
                                  A=list(spat.pred, 1,A1_pred),
                                  effects=list(spat.con=1:mesh$n,
                                               list(b0=1,
                                                    prof=pred_bathy,
                                                    year=1),
                                               spde1.idx),
                                  tag='pred.2000')

pred.extreme.occur.2016<-inla.stack(data=list(y=rep(NA,mesh$n)),
                                  A=list(spat.pred, 1,A1_pred),
                                  effects=list(spat.con=1:mesh$n,
                                               list(b0=1,
                                                    prof=pred_bathy,
                                                    year=17),
                                               spde1.idx),
                                  tag='pred.2016')

pred.extreme.occur.CRPS<-inla.stack(data=list(y=rep(NA,length(extremes_presence85))),
                                  A=list(spat.CRPS, 1,A1_CRPS ),
                                  effects=list(spat.con=1:mesh$n,
                                               list(b0=1,
                                                    prof=pred_bathy[extremes_presence85],
                                                    year=datos$year2[extremes_presence85]),
                                               spde1.idx),
                                  tag='pred.CRPS')

stack_pred_extreme_occur = inla.stack(stack_stage2,pred.extreme.bin.2000,pred.extreme.bin.2016,pred.extreme.bin.CRPS)

#form_extreme1= y~-1 + b0 +  f(bathy,model=spde1_NA)
form_extreme1= y~-1 + b0 +  prof 

stage3_fit=inla(form_extreme1,
                        family=c('binomial'),control.inla = list(strategy = method),
                        data=inla.stack.data(stack_stage3), control.compute=list(config=T),
                        control.predictor=list(A=inla.stack.A(stack_stage3), compute=TRUE),
                        verbose=T,num.threads=2)

save(stack_stage3,file=paste0(save.link,"/final5/stage3/stack_stage3.Rdata"))
save(stage3_fit,file=paste0(save.link,"/final5/stage3/stage3_fit.Rdata"))

# occur_extreme_pred_fine=inla(form_extreme1,
#                         family=c('binomial'),control.inla = list(strategy = method),
#                         data=inla.stack.data(stack_pred_extreme_occur), control.compute=list(config=T),
#                         control.predictor=list(A=inla.stack.A(stack_pred_extreme_occur),link=1, compute=TRUE),
#                         verbose=T,num.threads=2)
# 
# save(stack_pred_extreme_occur,file=paste0(save.link,"/final5/stage3/stack_pred_extreme_occur_fine.Rdata"))
# save(occur_extreme_pred_fine,file=paste0(save.link,"/final5/stage3/occur_extreme_pred_fine.Rdata"))

# plot(stage2_fit$summary.random$bathy[,1:2])
# plot(stage2_fit$summary.random$year[,1:2])
stage3_fit$summary.fixed


##############################
####### extreme values PARETO
##############################
gdp_extr=which(datos$extremes85==1)
extreme_data=datos[gdp_extr,]
extreme_data$bathy_stage4 = inla.group(extreme_data$PROFUNDIDAD,n=15,method = "cut")

save(extreme_data,file=paste0(save.link,"/final5/stage4/extreme_data.Rdata"))

### covariate 1D matern
mesh1d_ext <- inla.mesh.1d(seq(min(extreme_data$PROFUNDIDAD), max(extreme_data$PROFUNDIDAD), by = 10)) 

A1_ext <- inla.spde.make.A(mesh1d_ext, extreme_data$PROFUNDIDAD)
spde1_ext <- inla.spde2.pcmatern(mesh1d_ext, prior.range=c(80, 0.15), prior.sigma=c(2, 0.2))
spde1_NA_ext <- inla.spde2.pcmatern(mesh1d_ext, prior.range=c(80, NA), prior.sigma=c(2, 0.2))
spde1.ext.idx <- inla.spde.make.index("bathy", n.spde = spde1_ext$n.spde)

pred_bathy_ext = pred_bathy_NA
pred_bathy_ext[which(pred_bathy_ext>max(extreme_data$PROFUNDIDAD))] = max(extreme_data$PROFUNDIDAD)
A1_pred_ext <- inla.spde.make.A(mesh1d_ext,pred_bathy_ext)

i.gdp_extr <- inla.spde.make.A(mesh, loc=d[gdp_extr,])

A1_CRPS <- inla.spde.make.A(mesh1d_ext,datos$PROFUNDIDAD[extremes_presence85])
spat.CRPS <- inla.spde.make.A(mesh, loc=d[extremes_presence85,])

stack_stage4<-inla.stack(data=list(y=extreme_data$numbers),
                         A=list(i.gdp_extr, 1,A1_ext),
                         effects=list(spat=1:mesh$n,
                                      list(b0=1,
                                           year=extreme_data$year2,
                                           #month=datos$MonthShot[gdp_extr],
                                           prof=extreme_data$bathy_stage4,
                                           effort=extreme_data$INICIO.VIRADO),
                                      spde1.ext.idx),
                         tag='est')

pred.gdp.2000<-inla.stack(data=list(y=rep(NA,mesh$n)),
                          A=list(spat.pred, 1,A1_pred_ext),
                          effects=list(spat=1:mesh$n,
                                       list(b0=1,
                                            prof=pred_bathy_ext,
                                            year=1,
                                            effort=0.5),
                                       spde1.ext.idx),
                          tag='pred_gdp.2000')

pred.gdp.2016<-inla.stack(data=list(y=rep(NA,mesh$n)),
                          A=list(spat.pred, 1,A1_pred_ext),
                          effects=list(spat=1:mesh$n,
                                       list(b0=1,
                                            prof=pred_bathy_ext,
                                            year=17,
                                            effort=0.5),
                                       spde1.ext.idx),
                          tag='pred_gdp.2016')

pred.gdp.CRPS<-inla.stack(data=list(y=rep(NA,length(extremes_presence85))),
                                  A=list(spat.CRPS, 1,A1_CRPS ),
                                  effects=list(spat.con=1:mesh$n,
                                               list(b0=1,
                                                    prof=datos$PROFUNDIDAD[extremes_presence85],
                                                    year=datos$year2[extremes_presence85]),
                                               spde1.ext.idx),
                                  tag='pred.gdp.CRPS')

stack_pred_gdp = inla.stack(est.gdp_extr,pred.gdp.2000,pred.gdp.2016,pred.gdp.CRPS)

alfa = 0.5 ## median
form_gdp5= y~-1 + b0 +  year + f(bathy,model=spde1_NA_ext)+ offset(log(effort)) 

stage4_fit=inla(form_gdp5,
                   family=c('dgp'),
                   control.inla = list(strategy = method),
                   control.family=list(control.link=list(quantile=alfa),
                                       hyper = list(tail = list(
                                         prior = "pc.gevtail",
                                         param = c(1, 0.0, 0.5)))),
                   data=inla.stack.data(stack_stage4),
                   control.compute=list(config=T),
                   control.predictor=list(A=inla.stack.A(stack_stage4), compute=TRUE),
                   verbose=T, num.threads = 2)

save(stack_stage4,file=paste0(save.link,"/final5/stage4/stack_stage4.Rdata"))
save(stage4_fit,file=paste0(save.link,"/final5/stage4/stage4_fit.Rdata"))
#save(extreme_data,file=paste0(save.link,"/final5/stage4/extreme_data.Rdata"))


gdp_extr_pred=inla(form_gdp5,
                   family=c('dgp'),
                   control.inla = list(strategy = method),
                   control.family=list(control.link=list(quantile=alfa),
                                       hyper = list(tail = list(
                                         prior = "pc.gevtail",
                                         param = c(1, 0.0, 0.5)))),
                   data=inla.stack.data(stack_pred_gdp), 
                   #control.mode = list(result = gdp_extr_fit5, restart = TRUE),
                   control.compute=list(config=T),
                   control.predictor=list(A=inla.stack.A(stack_pred_gdp),link=1, compute=TRUE),
                   verbose=T, num.threads = 2)

save(stack_pred_gdp,file=paste0(save.link,"/final5/stage4/stack_pred_gdp.Rdata"))
save(gdp_extr_pred,file=paste0(save.link,"/final5/stage4/gdp_extr_pred.Rdata"))
save(extreme_data,file=paste0(save.link,"/final5/stage4/extreme_data.Rdata"))

plot(gdp_extr_pred)

idx_est_gdp = inla.stack.index(stack_pred_gdp,tag="est.gdp_extr")$data
idx_pred_gdp = inla.stack.index(stack_pred_gdp,tag="pred_gdp")$data

### fit estimation model
hist(gdp_extr_fit5$summary.fitted.values[idx_est_gdp,1])
### fit prediction model
hist(gdp_extr_pred$summary.fitted.values[idx_est_gdp,1])
### predict prediction model
hist(gdp_extr_pred$summary.fitted.values[idx_pred_gdp,1])
plot(gdp_extr_pred$summary.fitted.values[idx_pred_gdp,1],ylim=c(50,1000))

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
### !!!!!! esta predicci?n NO ESTA BIEN !!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








# ##################################################
# ####### bulk fit and occurrence ##########
# ####################################################
# no_extreme=which(datos$extremes85==0)
# no_extreme_data=datos[no_extreme,]
# no_extreme_data$bathy_no_stage4 = inla.group(no_extreme_data$PROFUNDIDAD,n=25,method = "cut")
# 
# ### covariate 1D matern
# mesh1d_NOext <- inla.mesh.1d(seq(min(no_extreme_data$PROFUNDIDAD), max(no_extreme_data$PROFUNDIDAD), by = 10)) 
# A1_NOext <- inla.spde.make.A(mesh1d_NOext, no_extreme_data$PROFUNDIDAD)
# 
# spde1_NOext <- inla.spde2.pcmatern(mesh1d_NOext, prior.range=c(80, 0.15), prior.sigma=c(2, 0.2))
# spde1_NA_NOext <- inla.spde2.pcmatern(mesh1d_NOext, prior.range=c(120, NA), prior.sigma=c(2, 0.2))
# spde1.NOext.idx <- inla.spde.make.index("bathy", n.spde = spde1_NOext$n.spde)
# 
# pred_bathy_NOext = pred_bathy_NA
# pred_bathy_NOext[which(pred_bathy_NOext>max(no_extreme_data$PROFUNDIDAD))] = max(no_extreme_data$PROFUNDIDAD)
# A1_pred_NOext <- inla.spde.make.A(mesh1d_NOext,pred_bathy_NOext)
# 
# A1_CRPS <- inla.spde.make.A(mesh1d_NOext,datos$PROFUNDIDAD[extremes_presence85])
# spat.CRPS <- inla.spde.make.A(mesh, loc=d[extremes_presence85,])
# 
# A.no_extreme <- inla.spde.make.A(mesh, loc=d[no_extreme,])
# spat.pred <- inla.spde.make.A(mesh)
# 
# 
# 
# est.numbers<-inla.stack(data=list(y=no_extreme_data$numbers),
#                         A=list(A.no_extreme, 1,A1_NOext),
#                         effects=list(spat.con=1:mesh$n,
#                                      list(b0=1,
#                                           prof=no_extreme_data$PROFUNDIDAD,
#                                           year=no_extreme_data$year2,
#                                           effort=no_extreme_data$INICIO.VIRADO),
#                                      spde1.NOext.idx),
#                         tag='est_numbers')
# 
# pred.numbers.2000<-inla.stack(data=list(y=rep(NA,mesh$n)),
#                               A=list(spat.pred, 1,A1_pred_NOext),
#                               effects=list(spat.con=1:mesh$n,
#                                            list(b0=1,
#                                                 prof=pred_bathy_NOext,
#                                                 year=1,
#                                                 effort=0.5),
#                                            spde1.NOext.idx),
#                               tag='pred_numbers.2000')
# 
# pred.numbers.2016<-inla.stack(data=list(y=rep(NA,mesh$n)),
#                               A=list(spat.pred, 1,A1_pred_NOext),
#                               effects=list(spat.con=1:mesh$n,
#                                            list(b0=1,
#                                                 prof=pred_bathy_NOext,
#                                                 year=17,
#                                                 effort=0.5),
#                                            spde1.NOext.idx),
#                               tag='pred_numbers.2016')
# 
# pred.numbers.CRPS<-inla.stack(data=list(y=rep(NA,length(extremes_presence85))),
#                           A=list(spat.CRPS, 1,A1_CRPS ),
#                           effects=list(spat.con=1:mesh$n,
#                                        list(b0=1,
#                                             prof=datos$PROFUNDIDAD[extremes_presence85],
#                                             year=datos$year2[extremes_presence85]),
#                                        spde1.NOext.idx),
#                           tag='pred.numbers.CRPS')
# 
# 
# stack_bulk_pred = inla.stack(est.numbers,pred.numbers.2000,pred.numbers.2016,pred.numbers.CRPS)
# # stack_occur_pred = inla.stack(est.occur,pred.occur)
# 
# # prof_no_extreme_pred = c()
# # mesh_bathy =over(mesh.points,bathy)$band1
# # group_bathy_sort = sort(unique(no_extreme_data$bathy_no_stage4))
# # for(i in 1:mesh$n){
# #   #if(mesh_bathy[i]>max(no_extreme_data$PROFUNDIDAD) | is.na(mesh_bathy[i])){
# #   if( is.na(mesh_bathy[i])){
# #     prof_no_extreme_pred[i] = NA
# #   }else{
# #     idx=findInterval(mesh_bathy[i],sort(unique(no_extreme_data$bathy_no_stage4)))+1
# #     prof_no_extreme_pred[i] = sort(unique(no_extreme_data$bathy_no_stage4))[idx]
# #   }
# # }
# 
# 
# 
# 
# #### INDEP MODEL
# # form_occur<-y~-1+b0+
# #   f(bathy,model=spde1)+
# #   f(year,model="rw2",hyper = prior_occur)+
# #   f(spat.bin,model=spde_occur)
# 
# form_numbers<-y~-1+b0+
#   f(bathy,model=spde1_NA_NOext)+
#   f(year,model = "rw2")+
#   f(spat.con,model=spde_npue) + 
#   offset(log(effort))
# 
# 
# 
# numbers_fit<-inla(form_numbers,
#                   family=c('zeroinflatednbinomial2'),
#                   #control.fixed=list(expand.factor.strategy="inla"),
#                   data=inla.stack.data(est.numbers), control.compute=list(config=T),
#                   control.predictor=list(A=inla.stack.A(est.numbers), compute=TRUE),
#                   verbose=T,num.threads=1,control.inla = list(strategy = method))
# 
# 
# 
# save(numbers_fit,file=paste0(save.link,"/final5/bulk_occur_model/numbers_fit.Rdata"))
# save(est.numbers,file=paste0(save.link,"/final5/bulk_occur_model/est.numbers.stack.Rdata"))
# 
# ### its giving me errors. Maybe prediction stacks related
# numbers_pred<-inla(form_numbers,
#                    family=c('zeroinflatednbinomial2'),
#                    #control.fixed=list(expand.factor.strategy="inla"),
#                    #control.inla=list(strategy='adaptive',int.strategy='eb'), ### Memory error 
#                    control.inla = list(strategy = method),
#                    data=inla.stack.data(stack_bulk_pred), control.compute=list(config=T),
#                    control.predictor=list(A=inla.stack.A(stack_bulk_pred),link=1, compute=TRUE),
#                    verbose=T,num.threads=1)
# 
# save(numbers_pred,file=paste0(save.link,"/final5/bulk_occur_model/numbers_pred.Rdata"))
# save(stack_bulk_pred,file=paste0(save.link,"/final5/bulk_occur_model/stack_bulk_pred.Rdata"))
# 
# 
# numbers_pred$summary.fixed
# 
# 
# plot(numbers_pred$summary.random$bathy[,1:2])
# 
# 
# 
# ## remove files for space
# rm(numbers_pred)
# rm(numbers_fit)
# rm(gdp_extr_fit5)
# rm(gdp_extr_pred)
# rm(occur_extreme_fitX)
# rm(occur_extreme_pred)
# 



########################################################################################################################
############################################## Combined temporal trends ##############################################
#####################################################################################################################

est.numbers<-inla.stack(data=list(y=no_extreme_data$numbers),
                        A=list(A.no_extreme, 1,A1_NOext),
                        effects=list(spat.con=1:mesh$n,
                                     list(b0=1,
                                          prof=no_extreme_data$PROFUNDIDAD,
                                          year=no_extreme_data$year2,
                                          effort=no_extreme_data$INICIO.VIRADO),
                                     spde1.NOext.idx),
                        tag='est_numbers')

est.extreme.bin<-inla.stack(data=list(y=cbind(datos$extremes85)),
                            A=list(spat.index, 1,A1),
                            effects=list(spat.bin=1:mesh$n,
                                         list(b0=1,
                                              prof=datos$PROFUNDIDAD,
                                              year=datos$year2), ## change 
                                         spde1.idx),
                            tag='est.extreme.bin')

est.gdp_extr<-inla.stack(data=list(y=extreme_data$numbers),
                         A=list(i.gdp_extr, 1,A1_ext),
                         effects=list(spat=1:mesh$n,
                                      list(b0=1,
                                           year=extreme_data$year2,
                                           #month=datos$MonthShot[gdp_extr],
                                           prof=extreme_data$bathy_stage4,
                                           effort=extreme_data$INICIO.VIRADO),
                                      spde1.ext.idx),
                         tag='est.gdp_extr')

load(file=paste0(save.link,"/final5/bulk_occur_model/est.numbers.stack.Rdata"))
load(file=paste0(save.link,"/final5/stage4/est.gdp_extr.Rdata"))
load(file=paste0(save.link,"/final5/stage3/est.extreme.bin.Rdata"))


### 3 different bathymetries

idx1 = which(pred_bathy>70 & pred_bathy<80)
idx2 = which(pred_bathy>145 & pred_bathy<155)
idx3 = which(pred_bathy>220 & pred_bathy<230)

pred_bathy[idx1[7]]
pred_bathy[idx2[2]]
pred_bathy[idx3[2]]

IDX=c(idx1[7],idx2[2],idx3[2])




pred_frame = data.frame(lon  =rep(mesh$loc[IDX,1],each=17),
                        lat = rep(mesh$loc[IDX,2],each=17),
                        prof = rep(pred_bathy[IDX],each=17),
                        year = rep(1:17,3),
                        y=rep(NA,17*3))    
spat.pred.IDX = inla.spde.make.A(mesh, loc=cbind(pred_frame$lon,pred_frame$lat))

A1_pred_NOext.IDX <- inla.spde.make.A(mesh1d,pred_frame$prof)
pred.numbers_IDX<-inla.stack(data=list(y=pred_frame$y),
                             A=list(spat.pred.IDX, 1,A1_pred_NOext.IDX),
                             effects=list(spat.con=(1:mesh$n),
                                          list(b0=1,
                                               prof=pred_frame$prof,
                                               year=pred_frame$year,
                                               effort=0.5),
                                          spde1.NOext.idx),
                             tag='pred_numbers')

A1_pred_ext.IDX <- inla.spde.make.A(mesh1d_ext,pred_frame$prof)
pred.gdp_IDX<-inla.stack(data=list(y=pred_frame$y),
                         A=list(spat.pred.IDX, 1,A1_pred_ext.IDX),
                         effects=list(spat=1:mesh$n,
                                      list(b0=1,
                                           prof=pred_frame$prof,
                                           year=pred_frame$year,
                                           effort=0.5),
                                      spde1.ext.idx),
                         tag='pred_gdp')

A1_pred.IDX <- inla.spde.make.A(mesh1d,pred_frame$prof)
pred.extreme.bin_IDX<-inla.stack(data=list(y=pred_frame$y),
                                 A=list(spat.pred.IDX, 1,A1_pred.IDX),
                                 effects=list(spat.con=1:mesh$n,
                                              list(b0=1,
                                                   prof=pred_frame$prof,
                                                   year=pred_frame$year),
                                              spde1.idx),
                                 tag='pred.extreme.bin')



trend_stack_bulk = inla.stack(est.numbers,pred.numbers_IDX)
trend_stack_ECE_bin = inla.stack(est.extreme.bin,pred.extreme.bin_IDX)
trend_stack_ECE = inla.stack(est.gdp_extr,pred.gdp_IDX)


form_numbers<-y~-1+b0+
  f(bathy,model=spde1)+
  f(year,model = "rw2")+
  f(spat.con,model=spde_npue) + 
  offset(log(effort))
form_extreme1= y~-1 + b0 +  f(bathy,model=spde1)
form_gdp5= y~-1 + b0 +  year + f(bathy,model=spde1_ext)+ offset(log(effort)) 
alfa = 0.5 ## median

bulk_trend<-inla(form_numbers,
                 family=c('zeroinflatednbinomial2'),
                 #control.fixed=list(expand.factor.strategy="inla"),
                 #control.inla=list(strategy='adaptive',int.strategy='eb'), ### Memory error 
                 data=inla.stack.data(trend_stack_bulk), control.compute=list(config=T),
                 control.predictor=list(A=inla.stack.A(trend_stack_bulk), compute=TRUE),
                 verbose=T,num.threads=1)

bulk_trend_idx = inla.stack.index(trend_stack_bulk,tag="pred_numbers")$data



ECE_bin_trend<-inla(form_extreme1,
                    family=c('binomial'),
                    #control.fixed=list(expand.factor.strategy="inla"),
                    #control.inla=list(strategy='adaptive',int.strategy='eb'), ### Memory error 
                    data=inla.stack.data(trend_stack_ECE_bin), control.compute=list(config=T),
                    control.predictor=list(A=inla.stack.A(trend_stack_ECE_bin), compute=TRUE),
                    verbose=T,num.threads=1)

ECE_bin_trend_idx = inla.stack.index(trend_stack_ECE_bin,tag="pred.extreme.bin")$data


ECE_trend=inla(form_gdp5,
               family=c('dgp'),
               control.inla = list(strategy = method),
               control.family=list(control.link=list(quantile=alfa),
                                   hyper = list(tail = list(
                                     prior = "pc.gevtail",
                                     param = c(1, 0.0, 0.5)))),
               data=inla.stack.data(trend_stack_ECE), 
               #control.mode = list(result = gdp_extr_fit5, restart = TRUE),
               control.compute=list(config=T),
               control.predictor=list(A=inla.stack.A(trend_stack_ECE), compute=TRUE),
               verbose=T, num.threads = 1)

ECE_trend_idx = inla.stack.index(trend_stack_ECE,tag="pred_gdp")$data


save(ECE_trend_idx,file=paste0(save.link,"/final5/bulk_occur_model/ECE_trend_idx.Rdata"))
save(ECE_trend,file=paste0(save.link,"/final5/bulk_occur_model/ECE_trend.Rdata"))
save(ECE_bin_trend_idx,file=paste0(save.link,"/final5/bulk_occur_model/ECE_bin_trend_idx.Rdata"))
save(ECE_bin_trend,file=paste0(save.link,"/final5/bulk_occur_model/ECE_bin_trend.Rdata"))
save(bulk_trend_idx,file=paste0(save.link,"/final5/bulk_occur_model/bulk_trend_idx.Rdata"))
save(bulk_trend,file=paste0(save.link,"/final5/bulk_occur_model/bulk_trend.Rdata"))



### combine ####

n = 100000
inv.logit = function(x){exp(x)/(1+exp(x))}

trend_post1 = trend_post2 = trend_post3 = list()
trend_mean1 = trend_mean2 = trend_mean3 = list()

for(i in 1:(length(IDX)*17)){
  
  #### 1- occurrence #### 
  pred.idx = bulk_trend_idx
  
  p0 = function(mu,alpha){1-(mu/(1+mu))^alpha}
  p1 = function(mu,alpha){(mu/(1+mu))^alpha}
  alpha=bulk_trend$summary.hyperpar[2,1]
  mu = inla.emarginal(exp,bulk_trend$marginals.fitted.values[[pred.idx[i]]])
  #mu = mean(exp(inla.rmarginal(n,bulk_trend$marginals.fitted.values[[pred.idx[i]]])))
  r_mu = exp(inla.rmarginal(n,bulk_trend$marginals.fitted.values[[pred.idx[i]]]))
  r_prob = sapply(r_mu,p1,alpha)
  prob = p1(mu=mu,alpha=alpha)
  prob0 = p0(mu=mu,alpha=alpha)
  
  #prob = mean((inla.rmarginal(n,occur_pred$marginals.fitted.values[[pred.idx[i]]])))
  #prob = occur_fit$summary.fitted.values[pred.idx[i],1]
  
  occurrence = rbinom(n=n,size=1,prob=prob)
  zeros = rep(0,n - sum(occurrence))
  n_zeros = length(zeros)
  n_presence = sum(occurrence)
  
  n == n_zeros + n_presence
  
  #### 2- extreme events occurrence ####
  pred.idx = ECE_bin_trend_idx
  
  #plot(ECE_bin_trend$marginals.fitted.values[[pred.idx[i]]])
  prob = mean(inv.logit(inla.rmarginal(n,ECE_bin_trend$marginals.fitted.values[[pred.idx[i]]])))
  #prob = ECE_bin_trend$summary.fitted.values[pred.idx[i],1]
  
  extreme_occurrence = rbinom(n=n_presence,size=1,prob=prob)
  extreme_occurrence_mean = rbinom(n=n,size=1,prob=prob)
  
  n_bulk = n_presence - sum(extreme_occurrence)
  n_bulk_mean = n - sum(extreme_occurrence_mean)
  
  n_extreme = sum(extreme_occurrence)
  n_extreme_mean = sum(extreme_occurrence_mean)
  
  n == n_zeros + n_bulk + n_extreme
  
  #### 3- bulk abundance ####
  #est.idx = inla.stack.index(stack_bulk_pred,tag = "pred_numbers")$data
  pred.idx = bulk_trend_idx
  
  #plot(bulk_trend$marginals.fitted.values[[pred.idx[i]]])
  #kk = which(bulk_trend$marginals.fitted.values[[pred.idx[i]]][,2]<0.001) ## remove those marginal values with extremely small density
  #plot(bulk_trend$marginals.fitted.values[[pred.idx[i]]][-kk,])
  bulk = exp(inla.rmarginal(n_bulk,bulk_trend$marginals.fitted.values[[pred.idx[i]]]))
  bulk_mean = exp(inla.rmarginal(n_bulk_mean,bulk_trend$marginals.fitted.values[[pred.idx[i]]]))
  #bulk = inla.rmarginal(n_bulk,bulk_trend$marginals.fitted.values[[pred.idx[i]]])
  
  
  #### 4- extreme ####
  pred.idx = ECE_trend_idx
  # 
  #plot(ECE_trend$marginals.fitted.values[[pred.idx[i]]])
  extremes = exp(inla.rmarginal(n_extreme,ECE_trend$marginals.fitted.values[[pred.idx[i]]]))
  extremes_mean = exp(inla.rmarginal(n_extreme_mean,ECE_trend$marginals.fitted.values[[pred.idx[i]]]))
  # 
  #extremes = sample(predict_gdp_2000[[i]],n_extreme,replace = T)
  
  #### #### #### ####
  ### posterior ####
  if(i<=17){
    trend_post1[[i]]=c(zeros,bulk,extremes)
    trend_mean1[[i]]=r_prob*c(bulk_mean,extremes_mean)
  }
  if(i>17 & i<=34){
    trend_post2[[i-17]]=c(zeros,bulk,extremes)
    trend_mean2[[i-17]]=r_prob*c(bulk_mean,extremes_mean)
  }
  if(i>34 & i<=51){
    trend_post3[[i-34]]=c(zeros,bulk,extremes)
    trend_mean3[[i-34]]=r_prob*c(bulk_mean,extremes_mean)
  }
  
  #spliced_cumulative[[i]] = cumsum(sort(c(zeros,bulk,extremes)))
  #post = data.frame(values = c(zeros,bulk,extremes),cumulative = cumsum(sort(c(zeros,bulk,extremes))),x=1:length(c(zeros,bulk,extremes)))
  
  print(i)
}

mean_abu1=unlist(lapply(trend_mean1,FUN=mean))
q025_abu1=unlist(lapply(trend_mean1,FUN=function(x){quantile(x,probs=.025,na.rm=T)}))
q975_abu1=unlist(lapply(trend_mean1,FUN=function(x){quantile(x,probs=.975,na.rm=T)}))

trend_abundance1=data.frame(mean_abu1,q025_abu1,q975_abu1,
                            year=1:17,
                            loc=IDX[1])

t1 = ggplot(trend_abundance1, aes(x=year, y=mean_abu1)) + 
  geom_errorbar(aes(ymin=q025_abu1, ymax=q975_abu1), colour="black", width=.1) +
  geom_line() +
  geom_point(size=3) + ggtitle(paste0("mesh loc =",IDX[1]))


mean_abu2=unlist(lapply(trend_mean2,FUN=mean))
q025_abu2=unlist(lapply(trend_mean2,FUN=function(x){quantile(x,probs=.025,na.rm=T)}))
q975_abu2=unlist(lapply(trend_mean2,FUN=function(x){quantile(x,probs=.975,na.rm=T)}))

trend_abundance2=data.frame(mean_abu2,q025_abu2,q975_abu2,
                            year=1:17,
                            loc=IDX[2])

t2 = ggplot(trend_abundance2, aes(x=year, y=mean_abu2)) + 
  geom_errorbar(aes(ymin=q025_abu2, ymax=q975_abu2), colour="black", width=.1) +
  geom_line() +
  geom_point(size=3) + ggtitle(paste0("mesh loc =",IDX[2]))

mean_abu3=unlist(lapply(trend_mean3,FUN=mean))
q025_abu3=unlist(lapply(trend_mean3,FUN=function(x){quantile(x,probs=.025,na.rm=T)}))
q975_abu3=unlist(lapply(trend_mean3,FUN=function(x){quantile(x,probs=.975,na.rm=T)}))

trend_abundance3=data.frame(mean_abu3,q025_abu3,q975_abu3,
                            year=1:17,
                            loc=IDX[3])

t3 = ggplot(trend_abundance3, aes(x=year, y=mean_abu3)) + 
  geom_errorbar(aes(ymin=q025_abu3, ymax=q975_abu3), colour="black", width=.1) +
  geom_line() +
  geom_point(size=3) + ggtitle(paste0("mesh loc =",IDX[3]))


gridExtra::grid.arrange(t1,t2,t3,ncol=1)

png(filename = "C:/Users/ip30/Dropbox (Personal)/Marie Curie/projects/MC objectives/Merlu_extreme/final5/report/trends.png",width=600,height = 1100)
gridExtra::grid.arrange(t1,t2,t3,ncol=1)
dev.off()


save(trend_abundance1,file = "C:/Users/ip30/Dropbox (Personal)/Marie Curie/projects/MC objectives/Merlu_extreme/final5/report/trend_abundance1.Rdata")
save(trend_abundance2,file = "C:/Users/ip30/Dropbox (Personal)/Marie Curie/projects/MC objectives/Merlu_extreme/final5/report/trend_abundance2.Rdata")
save(trend_abundance3,file = "C:/Users/ip30/Dropbox (Personal)/Marie Curie/projects/MC objectives/Merlu_extreme/final5/report/trend_abundance3.Rdata")

