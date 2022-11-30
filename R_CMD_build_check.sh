
rm graphicalExtremes_0.1.0.9000.tar.gz

rm -rf graphicalExtremes.Rcheck

export RCMD=R

# %RCMD% CMD build --no-build-vignettes .
$RCMD CMD build .

# %RCMD% CMD check --ignore-vignettes --as-cran graphicalExtremes_0.2.0.tar.gz
$RCMD CMD check --as-cran graphicalExtremes_0.2.0.tar.gz
