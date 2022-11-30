
del graphicalExtremes_0.2.0.tar.gz

rmdir graphicalExtremes.Rcheck /S /Q

@REM set RCMD=R
set RCMD="C:\Program Files\R\R-devel\bin\R.exe"

@REM %RCMD% CMD build --no-build-vignettes .
%RCMD% CMD build .

@REM %RCMD% CMD check --ignore-vignettes --as-cran graphicalExtremes_0.2.0.tar.gz
%RCMD% CMD check --as-cran graphicalExtremes_0.2.0.tar.gz
