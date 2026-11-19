Sys.setenv(STLMM_REFERENCE_ONLY='1')
source('tests/native/latent-reference.R')
Sys.unsetenv('STLMM_REFERENCE_ONLY')
adj <- matrix(0,5,5,dimnames=list(letters[1:5],letters[1:5]))
adj[cbind(1:4,2:5)] <- 1;adj <- adj+t(adj);g <- car_graph(adj)
set.seed(73322)
d <- expand.grid(area=letters[1:4],time=c(1,2,4),replicate=1:2)
d$parent <- rep(c('A','B'),length.out=nrow(d));d$group <- factor(d$area,levels=letters[1:5])
d$x <- rnorm(nrow(d));d$y <- rnorm(nrow(d));d$y[c(2,7)] <- NA_real_
d$v <- seq(.2,1.2,length.out=nrow(d));d$offset_value <- rnorm(nrow(d))
cases <- list(car='car(area,graph=g)',leroux='car(area,graph=g,car_model="leroux")',
  dagar='dagar(area,graph=g,ordering=c(3,1,5,2,4))',
  car_time='car_time(area,time,graph=g)',car_time_exp='car_time(area,time,graph=g,time_model="exp")',
  dagar_time='dagar_time(area,time,graph=g,ordering=c(3,1,5,2,4))',
  dagar_time_exp='dagar_time(area,time,graph=g,time_model="exp",ordering=c(3,1,5,2,4))')
results <- list()
for(id in names(cases)) for(p0 in c(FALSE,TRUE)) {
  nm <- if(id=='leroux') 'car' else sub('_exp$','',id)
  time <- grepl('_time',id);exponential <- grepl('_exp$',id)
  ctrl <- list(sigma_sq=fixed(.4),rho=fixed(.6))
  if(time) ctrl[[if(exponential) 'lambda' else 'phi']] <- fixed(.3)
  form <- as.formula(paste('y ~',if(p0) '0' else 'x',
     '+ iid(group,~1+x,by=parent) +',cases[[id]],
     '+ offset(offset_value) + resid(model="fixed",variance=v)'))
  fit <- suppressWarnings(stLMM(form,d,starting=setNames(list(ctrl),paste0(nm,'_1')),
    priors=list(iid_1=list(covariance=iw(4,diag(2)))),n_samples=1,verbose=FALSE))
  be <- fit$backend
  for(b in seq_along(be$iid_covariance)) be$iid_covariance[[b]]$starting <- matrix(c(.7,-.2,-.2,.4),2)
  ref <- reference(be)
  set.seed(73322+length(results))
  out <- .Call('stLMM_latent_probe',be,6000L,PACKAGE='stLMM')
  error <- max(abs(out$H-ref$H),abs(out$M-ref$M),abs(out$Vinv-ref$Vi),
    abs(out$logs[1:3]-c(ref$logdet,ref$quad,ref$logdet_iid)),abs(out$mean-ref$mu),
    abs(out$Xprec-ref$Xprec),abs(out$Xrhs-ref$Xrhs))
  zz <- sweep(out$draws,2,ref$mu)%*%t(chol(ref$M))
  zmean <- max(abs(colMeans(zz)))*sqrt(nrow(zz))
  cv <- cov(zz); zcov <- abs(cv-diag(ncol(cv)))*sqrt(nrow(zz)-1);diag(zcov)<-diag(zcov)/sqrt(2)
  stopifnot(error<1e-8,zmean<6.5,max(zcov)<6.5,isTRUE(out$restored),abs(out$logs[6]-out$logs[7])<1e-10)
  results[[paste(id,p0)]] <- c(error=error,zmean=zmean,zcov=max(zcov),check_known_gaussian_chain(be,ref))
}
saveRDS(results,args[2]);print(do.call(rbind,results));cat(length(results),'areal reference cases passed\n')
