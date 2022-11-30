
del graphicalExtremes_0.1.0.9000.tar.gz

rmdir graphicalExtremes.Rcheck /S /Q

set RCMD=R

@REM %RCMD% CMD build --no-build-vignettes .
%RCMD% CMD build .

@REM %RCMD% CMD check --ignore-vignettes --as-cran graphicalExtremes_0.2.0.tar.gz
%RCMD% CMD check --as-cran graphicalExtremes_0.2.0.tar.gz
