
del graphicalExtremes_0.1.0.9000.tar.gz

rmdir graphicalExtremes.Rcheck /S /Q

R CMD build --no-build-vignettes .

R CMD check --ignore-vignettes --as-cran graphicalExtremes_0.1.0.9000.tar.gz
