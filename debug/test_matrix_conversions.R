
devtools::load_all()

set.seed(1)

options(graphicalExtremes.tol.small = 1e-10)

d <- 10
Gamma <- generate_random_Gamma(d)
# m <- generate_random_model(d)
# Gamma <- m$Gamma

reprGamma <- list(
    name = 'Gamma',
    k = NULL,
    full = FALSE
)

repr0 <- reprGamma
repr0[['M']] <- Gamma

generate_random_representation <- function(d=1){
    repNames <- c('Gamma', 'Sigma', 'Theta')
    repName <- sample(repNames, 1, prob = c(0.2, 0.4, 0.4))
    if(repName == 'Gamma'){
        return(list(
            name = repName,
            k = NULL,
            full = FALSE
        ))
    }
    kIsNull <- sample(c(TRUE, FALSE), 1, prob = c(0.3, 0.7))
    if(kIsNull){
        return(list(
            name = repName,
            k = NULL,
            full = FALSE
        ))
    }
    k <- sample(seq(d), 1)
    full <- sample(c(TRUE, FALSE), 1)
    return(list(
        name = repName,
        k = k,
        full = full
    ))
}

convertRepresentation <- function(repr1, repr2, M = NULL, check=TRUE){
    if(!is.null(M)){
        repr1[['M']] <- M
    }
    M2 <- matrix2matrix(
        repr1$M,
        repr1$name,
        repr2$name,
        repr1$k,
        repr2$k,
        repr1$full,
        repr2$full,
        check = check
    )
    repr2[['M']] <- M2
    return(repr2)
}

convertToRandomRepresentation <- function(repr1){
    d <- computeD(repr1$M, repr1$k, repr1$full)
    repr2 <- generate_random_representation(d)
    repr2b <- convertRepresentation(repr1, repr2)
    return(repr2b)
}

checkReprGamma <- function(repr, check = TRUE){
    repr2 <- convertRepresentation(repr, reprGamma, check = check)
    return(max(abs(repr2$M - repr0$M)))
}

compareRepr <- function(repr1, repr2, check = TRUE){
    repr2b <- convertRepresentation(repr2, repr1, check = check)
    return(max(abs(repr1$M - repr2b$M)))
}

repr2string <- function(repr){
    ret <- paste0(
        repr$name,
        ' --- k=',
        format(repr$k),
        ' --- full=',
        repr$full
    )
}

maxErr <- 0
nIter <- 5000
repr <- repr0
reprList <- rep(list(NULL), nIter)
for(i in seq_len(nIter)){
    repr2 <- generate_random_representation(d)
    cat(i, ': ', repr2string(repr2), ':  ', sep='')
    repr <- convertRepresentation(repr, repr2, check = TRUE)
    d <- nrow(repr$M)
    err <- checkReprGamma(repr, FALSE)
    maxErr <- max(maxErr, err)
    cat(err, ' --- d=', d, '\n', sep='')
    reprList[[i]] <- repr
    if(maxErr > 1){
        stop('Big max err!')
    }
}

maxCompErr <- 0
nComps <- 500
cat('Comparing random pairs of representations...\n')
for(i in seq_len(nComps)){
    repr1 <- sample(reprList, 1)[[1]]
    repr2 <- sample(reprList, 1)[[1]]
    err <- compareRepr(repr1, repr2, check = TRUE)
    maxCompErr <- max(maxCompErr, err)
    cat(
        i,
        ': ',
        repr2string(repr1),
        '  <-->  ',
        repr2string(repr2),
        ': ',
        err,
        '\n',
        sep = ''
    )
}

cat('Max err:', maxErr, '\n')

cat('Max comp err:', maxCompErr, '\n')

cat('Done.\n')


