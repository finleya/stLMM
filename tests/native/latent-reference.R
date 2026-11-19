# Known-covariance checks against independent dense observation-space formulas.
# Uses a validation-only C++ probe in an isolated installation, not the public API.
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
.libPaths(c(args[1], .libPaths()))
suppressPackageStartupMessages(library(stLMM))
suppressPackageStartupMessages(library(Matrix))
set.seed(73301)
n <- 24L
d <- data.frame(y = rnorm(n), x = rnorm(n), z = 0,
  group = factor(rep(letters[1:4], each = 6), levels = letters[1:5]),
  parent = rep(c('A','B'), each = 12), child = rep(1:3, 8),
  crossed = rep(1:6, 4), time = rep(1:8, each = 3),
  lon = rep(runif(8), each = 3), lat = rep(runif(8), each = 3),
  v = seq(.2, 1.3, length.out = n))
d$nested <- interaction(d$parent, d$child, drop = TRUE)
d$xx <- d$x
d$y[c(2,9)] <- NA_real_
forms <- c(scalar = 'iid(group)', independent = 'iid(group, ~ 1 + x, cov_model="independent")',
  legacy = 'iid(group) + x:iid(group)', covariance = 'iid(group, ~ 1 + x)',
  crossed = 'iid(group, ~ 1 + x) + iid(crossed)',
  nested = 'iid(parent) + iid(nested, ~ 1 + x)',
  by = 'iid(child, ~ 1 + x, by=parent)', by_scalar = 'iid(child, by=parent)',
  collinear = 'iid(group, ~ 1 + x + xx + z)',
  scalar_covariance = 'iid(group, ~ 0 + x)')
get_priors <- function(kind) switch(kind,
  scalar = list(iid_1=list(sigma_sq=ig(3,1))),
  independent = list(iid_1=list(sigma_sq=list('(Intercept)'=ig(3,1),x=ig(4,.5)))),
  legacy = list(iid_1=list(sigma_sq=ig(3,1)),iid_2=list(sigma_sq=ig(4,.5))),
  crossed = list(iid_1=list(covariance=iw(4,diag(2))),iid_2=list(sigma_sq=ig(3,1))),
  nested = list(iid_1=list(sigma_sq=ig(3,1)),iid_2=list(covariance=iw(4,diag(2)))),
  by_scalar = list(iid_1=list(sigma_sq=ig(3,1))),
  collinear = list(iid_1=list(covariance=iw(6,diag(4)))),
  scalar_covariance = list(iid_1=list(covariance=iw(3,matrix(1)))),
  list(iid_1=list(covariance=iw(4,diag(2)))))

reference <- function(be) {
  Z <- as.matrix(be$Z_obs)
  D <- matrix(0,ncol(Z),ncol(Z))
  for(b in seq_along(be$iid_covariance)) {
    idx <- which(be$re$block_id == b)
    cb <- be$iid_covariance[[b]]
    if(is.null(cb)) D[cbind(idx,idx)] <- be$sigma_sq_re_starting[b] else {
      for(g in seq_len(length(idx)/cb$dim)) {
        jj <- idx[(g-1)*cb$dim+seq_len(cb$dim)]
        D[jj,jj] <- cb$starting
      }
    }
  }
  As <- Cs <- list()
  for(t in be$process_terms_obs) {
    gr <- be$graphs[[t$graph_index]]
    k <- t$n_node
    A <- matrix(0,nrow(Z),k)
    A[cbind(seq_len(nrow(Z)), t$map)] <- if(is.null(t$x)) 1 else t$x
    if(gr$graph_type == 'ar1') {
      C <- t$sigma_sq_starting * t$theta_starting^abs(outer(seq_len(k),seq_len(k),'-'))
    } else if(gr$graph_type %in% c('car','car_time','dagar','dagar_time')) {
      ns <- if(is.null(gr$n_space)) k else gr$n_space
      rho <- t$theta_starting[1]
      if(grepl('^car',gr$graph_type)) {
        adj <- matrix(0,ns,ns)
        adj[cbind(gr$edge_i,gr$edge_j)] <- gr$edge_w
        adj <- adj+t(adj)
        Qs <- if(identical(gr$car_model,'leroux'))
          (1-rho)*diag(ns)+rho*(diag(gr$degree)-adj) else diag(gr$degree)-rho*adj
      } else {
        B <- diag(ns); prec <- numeric(ns)
        for(i in seq_len(ns)) {
          ni <- gr$parent_count[i]
          jj <- if(ni) gr$parent_index[gr$parent_start[i]+seq_len(ni)]+1L else integer()
          den <- 1+(ni-1)*rho^2
          if(ni) B[i,jj] <- -rho/den
          prec[i] <- den/(1-rho^2)
        }
        Qs <- crossprod(B,B*prec)
      }
      C <- solve(Qs)
      if(grepl('_time$',gr$graph_type)) {
        dt <- abs(outer(gr$time_support,gr$time_support,'-'))
        Ct <- if(gr$time_model=='ar1')
          t$theta_starting[2]^abs(outer(seq_len(gr$n_time),seq_len(gr$n_time),'-')) else
          exp(-t$theta_starting[2]*dt)
        C <- kronecker(C,Ct)
      }
      C <- t$sigma_sq_starting*C
    } else {
      coords <- if(gr$graph_type=='nngp') gr$coords_ord else gr$coords
      C <- t$sigma_sq_starting * exp(-t$theta_starting*as.matrix(dist(coords)))
      if(gr$graph_type == 'nngp') {
        B <- diag(k); F <- numeric(k)
        for(i in seq_len(k)) {
          ni <- gr$nnIndxLU[k+i]
          jj <- if(ni) gr$nnIndx[gr$nnIndxLU[i]+seq_len(ni)]+1L else integer()
          if(ni) B[i,jj] <- -solve(C[jj,jj,drop=FALSE],C[jj,i])
          F[i] <- C[i,i] + sum(B[i,jj]*C[jj,i])
        }
        C <- solve(crossprod(B,B/F))
      }
    }
    As[[length(As)+1L]] <- A; Cs[[length(Cs)+1L]] <- C
  }
  H <- do.call(cbind,c(As,list(Z)))
  K <- as.matrix(bdiag(c(Cs,list(D))))
  var <- if(be$residual_model$type=='fixed_variance') be$residual_model$variance else rep(be$tau_sq_starting,nrow(Z))
  # The backend stores fixed precision on the observed rows.
  if(be$residual_model$type=='fixed_variance') var <- 1/be$residual_model$obs_precision_obs
  W <- 1/var
  V <- diag(var) + H %*% K %*% t(H)
  Vi <- solve(V)
  M <- solve(K) + crossprod(H,H*W)
  cond <- solve(M)
  list(H=H,K=K,W=W,M=M,Vi=Vi,mu=drop(cond%*%crossprod(H,be$y_model_obs*W)),cond=cond,
       logdet=as.numeric(determinant(V,logarithm=TRUE)$modulus),
       quad=drop(crossprod(be$y_model_obs,Vi%*%be$y_model_obs)),
       logdet_iid=-as.numeric(determinant(D,logarithm=TRUE)$modulus),
       Xprec=crossprod(be$X_obs,Vi%*%be$X_obs),Xrhs=drop(crossprod(be$X_obs,Vi%*%be$y_model_obs)))
}
check_known_gaussian_chain <- function(be,ref) {
  p <- be$p; q <- be$q; nw <- ncol(ref$H)-q
  perm <- c(nw+seq_len(q),seq_len(nw))
  design <- cbind(be$X_obs,ref$H[,perm,drop=FALSE])
  J <- crossprod(design,design*ref$W)
  ii <- p+seq_len(ncol(ref$H))
  J[ii,ii] <- J[ii,ii]+solve(ref$K)[perm,perm]
  mu <- drop(solve(J,crossprod(design,be$y_model_obs*ref$W)))
  draws <- .Call('stLMM_known_covariance_chain',be,3000L,TRUE,PACKAGE='stLMM')
  zz <- sweep(draws,2,mu)%*%t(chol(J))
  zmean <- max(abs(colMeans(zz)))*sqrt(nrow(zz))
  cv <- cov(zz); zcov <- abs(cv-diag(ncol(cv)))*sqrt(nrow(zz)-1)
  diag(zcov) <- diag(zcov)/sqrt(2)
  stopifnot(zmean<6.5,max(zcov)<6.5)
  c(chain_mean=zmean,chain_cov=max(zcov))
}
if (!identical(Sys.getenv("STLMM_REFERENCE_ONLY"), "1")) {
results <- list(); case <- 0L
for(kind in names(forms)) for(structure in c('none','ar1','nngp','gp')) for(weighted in c(FALSE,TRUE)) {
  case <- case+1L
  structured <- switch(structure, none='', ar1='+ ar1(time) + x:ar1(time)',
    nngp='+ nngp(lon, lat, m=3, cov_model="exp")',gp='+ gp(lon, lat, cov_model="exp")')
  st <- switch(structure, none=list(),
    ar1=list(ar1_1=list(sigma_sq=fixed(.4),phi=fixed(.6)),ar1_2=list(sigma_sq=fixed(.25),phi=fixed(-.3))),
    nngp=list(nngp_1=list(sigma_sq=fixed(.4),phi=fixed(2))),
    gp=list(gp_1=list(sigma_sq=fixed(.4),phi=fixed(2))))
  form <- as.formula(paste('y ~ x +',forms[[kind]],structured,
      if(weighted) '+ resid(model="fixed", variance=v)' else ''))
  fit <- suppressWarnings(stLMM(form,d,starting=st,priors=get_priors(kind),n_samples=1,verbose=FALSE))
  be <- fit$backend
  # Dense, non-diagonal, unequal covariance blocks and scalar variances.
  for(b in seq_along(be$iid_covariance)) {
    cb <- be$iid_covariance[[b]]
    if(is.null(cb)) be$sigma_sq_re_starting[b] <- .2 + .13*b else {
      L <- matrix(seq_len(cb$dim^2)/20,cb$dim); L[upper.tri(L)] <- 0
      diag(L) <- seq(.6,1.3,length.out=cb$dim)+b/10
      be$iid_covariance[[b]]$starting <- tcrossprod(L)
    }
  }
  ref <- reference(be)
  set.seed(73301+case)
  tm <- system.time(out <- .Call('stLMM_latent_probe',be,6000L,PACKAGE='stLMM'))['elapsed']
  err <- c(H=max(abs(out$H-ref$H)), M=max(abs(out$M-ref$M)),
    Vi=max(abs(out$Vinv-ref$Vi)),mu=max(abs(out$mean-ref$mu)),
    logs=max(abs(out$logs[1:3]-c(ref$logdet,ref$quad,ref$logdet_iid))),
    Xprec=max(abs(out$Xprec-ref$Xprec)),Xrhs=max(abs(out$Xrhs-ref$Xrhs)))
  stopifnot(max(err)<1e-8,isTRUE(out$restored),abs(out$logs[6]-out$logs[7])<1e-10)
  # Conditional draws whitened by independently computed posterior precision.
  # R chol(M)=U, so row draws must multiply U' to have unit covariance.
  zz <- sweep(out$draws,2,ref$mu) %*% t(chol(ref$M))
  zmean <- max(abs(colMeans(zz)))*sqrt(nrow(zz))
  cv <- cov(zz); zcov <- abs(cv-diag(ncol(cv)))*sqrt(nrow(zz)-1)
  diag(zcov) <- diag(zcov)/sqrt(2)
  stopifnot(zmean<6.5,max(zcov)<6.5)
  id <- paste(kind,structure,weighted,sep='/')
  results[[id]] <- c(err,zmean=zmean,zcov=max(zcov),check_known_gaussian_chain(be,ref),
                    elapsed=tm,nnzM=out$logs[4],nnzL=out$logs[5])
  if(case %% 10L==0L) cat(case,'reference cases passed\n')
}
saveRDS(results,args[2])
print(apply(do.call(rbind,results),2,max))
cat(length(results),'dense-reference and conditional-draw checks passed\n')

}
