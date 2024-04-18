
del graphicalExtremes_0.3.2.tar.gz

rmdir graphicalExtremes.Rcheck /S /Q

@REM Make sure dependencies are installed before!
@REM set RCMD=R
set RCMD="C:\Program Files\R\R-devel\bin\R.exe"
@REM set RCMD="%USERPROFILE%\AppData\Local\Programs\R\R-devel\bin\R.exe"

@REM %RCMD% -e "unlink(Sys.getenv('R_LIBS_USER'),recursive=TRUE,force=TRUE)"
@REM %RCMD% -e "dir.create(Sys.getenv('R_LIBS_USER'),recursive=TRUE)"
@REM %RCMD% -e "install.packages(c('remotes','ggplot2','glasso','V8'))"
@REM %RCMD% -e "remotes::install_deps(dependencies=TRUE)"

@REM %RCMD% CMD build --no-build-vignettes .
%RCMD% CMD build .

@REM %RCMD% CMD check --ignore-vignettes --as-cran graphicalExtremes_0.2.0.tar.gz
%RCMD% CMD check --as-cran graphicalExtremes_0.3.2.tar.gz
