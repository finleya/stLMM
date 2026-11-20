# Run from a source checkout: Rscript tests/native/run-latent-checks.R [package-path]
# Build an isolated package with native probes; the installed public API is untouched.
args <- commandArgs(TRUE)
package <- normalizePath(if(length(args)) args[1] else '.', mustWork=TRUE)
stopifnot(file.exists(file.path(package,'DESCRIPTION')))
work <- tempfile('stlmm-latent-checks-',tmpdir=dirname(tempdir()))
dir.create(work)
work <- normalizePath(work)
dir.create(file.path(work,'lib'))
run <- function(executable, arguments, log) {
  status <- system2(executable, arguments, stdout=log, stderr=log)
  if(status != 0L) {
    cat(tail(readLines(log,warn=FALSE),30),sep='\n')
    stop('Native validation failed; see ',log,call.=FALSE)
  }
}
cat('Validation directory:',work,'\n')
setwd(work)
run(file.path(R.home('bin'),'R'),
    c('CMD','build','--no-build-vignettes','--no-manual',shQuote(package)),
    file.path(work,'build.log'))
archive <- list.files(work,pattern='^stLMM_.*[.]tar[.]gz$',full.names=TRUE)
stopifnot(length(archive)==1L)
untar(archive,exdir=work)
src <- file.path(work,'stLMM','src')
stopifnot(file.copy(file.path(package,'tests/native/latent-probe.cpp'),src))
init <- file.path(src,'init.c')
code <- paste(readLines(init),collapse='\n')
code <- sub('static const R_CallMethodDef CallEntries[] = {',paste0(
  'extern SEXP stLMM_latent_probe(SEXP, SEXP);\n',
  'extern SEXP stLMM_known_covariance_chain(SEXP, SEXP, SEXP);\n',
  'extern SEXP stLMM_iid_mh_probe(SEXP);\n',
  'extern SEXP stLMM_iid_density_probe(SEXP, SEXP);\n',
  'static const R_CallMethodDef CallEntries[] = {\n',
  '    {"stLMM_latent_probe", (DL_FUNC) &stLMM_latent_probe, 2},\n',
  '    {"stLMM_iid_mh_probe", (DL_FUNC) &stLMM_iid_mh_probe, 1},\n',
  '    {"stLMM_iid_density_probe", (DL_FUNC) &stLMM_iid_density_probe, 2},\n',
  '    {"stLMM_known_covariance_chain", (DL_FUNC) &stLMM_known_covariance_chain, 3},'),
  code,fixed=TRUE)
writeLines(code,init)
run(file.path(R.home('bin'),'R'),
    c('CMD','INSTALL','--preclean',shQuote(paste0('--library=',file.path(work,'lib'))),
      shQuote(file.path(work,'stLMM'))),file.path(work,'install.log'))
setwd(package)
for(test in c('latent-reference','latent-areal-reference','iid-scalar-reference',
             'iid-proposal-reference','iid-iw-reference','iid-posterior-reference',
             'iid-augmented-reference','iid-augmented-schedules',
             'iid-augmented-posterior-reference')) {
  log <- file.path(work,paste0(test,'.log'))
  run(file.path(R.home('bin'),'Rscript'),
      shQuote(c(file.path(package,'tests/native',paste0(test,'.R')),
                file.path(work,'lib'),file.path(work,paste0(test,'.rds')))),log)
  cat(tail(readLines(log),1),'\n')
}
cat('Native latent-system checks passed. Results:',work,'\n')
