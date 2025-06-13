
rm graphicalExtremes_0.3.4.tar.gz

rm -rf graphicalExtremes.Rcheck

export RCMD=R

# %RCMD% CMD build --no-build-vignettes .
# $RCMD -e "install.packages(c('remotes','ggplot2','glasso','V8'))"
# $RCMD -e "remotes::install_deps(dependencies=TRUE)"

$RCMD CMD build .

# %RCMD% CMD check --ignore-vignettes --as-cran graphicalExtremes_0.2.0.tar.gz
$RCMD CMD check --as-cran graphicalExtremes_0.3.4.tar.gz
