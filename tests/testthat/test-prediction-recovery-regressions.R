test_that("Matern prediction correlations agree with closed forms", {
  for(nu in c(.5,1.5,2.5)) for(h in c(0,1e-8,.01,.3,1,10,100)){
    x <- .8*h
    expected <- exp(-x)*switch(as.character(nu),"0.5"=1,"1.5"=1+x,"2.5"=1+x+x*x/3)
    expect_equal(stLMM:::gp_correlation("matern",c(phi=.8,nu=nu),h,0),expected,tolerance=1e-13)
  }
})

test_that("recovery summaries distinguish w draws from w_node_order metadata", {
  d <- data.frame(y=c(.1,.3,-.2,.4,.2),time=1:5)
  f <- stLMM(y~ar1(time),d,n_samples=10,verbose=FALSE,
    starting=list(resid=list(tau_sq=fixed(.3)),ar1_1=list(sigma_sq=fixed(.2),phi=fixed(.4))))
  r <- recover(f)
  expect_output(print(summary(r)),"Recovery:")
  expect_output(print(summary(r,include_w=TRUE)),"w:")
})

test_that("saved fits recover in a fresh R session without loading Matrix explicitly", {
  d <- data.frame(y=c(.1,.3,-.2,.4,.2),time=1:5)
  f <- stLMM(y~ar1(time),d,n_samples=5,verbose=FALSE,
    starting=list(resid=list(tau_sq=fixed(.3)),ar1_1=list(sigma_sq=fixed(.2),phi=fixed(.4))))
  saved <- tempfile(fileext=".rds"); script <- tempfile(fileext=".R")
  on.exit(unlink(c(saved,script)),add=TRUE)
  saveRDS(f,saved)
  quote_r <- function(x) paste(deparse(x),collapse="")
  writeLines(c(paste0("library(stLMM,lib.loc=",quote_r(dirname(find.package("stLMM"))),")"),
               paste0("f <- readRDS(",quote_r(saved),")"),
               "r <- recover(f)",
               "stopifnot(nrow(r$w_samples$ar1_1)==5L)",
               "cat('Fresh-session recovery succeeded\\n')"),script)
  output <- system2(file.path(R.home("bin"),"Rscript"),c("--vanilla",shQuote(script)),stdout=TRUE,stderr=TRUE)
  expect_null(attr(output,"status"))
  expect_true(any(grepl("Fresh-session recovery succeeded",output,fixed=TRUE)))
})
