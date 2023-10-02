

# This script converts some (text) data in ./data-raw to rds objects in ./inst/extdata.
# The ./data-raw is ignored from the package build.

basenames <- c(
    'Texas_cluster_days_train_data',
    'Texas_cluster_days_test_data',
    'Texas_cluster_IATAS_all',
    'Texas_cluster_IATAS_cluster'
)

for(bn in basenames){
    txtPath <- file.path('data-raw', paste0(bn, '.txt'))
    rdsPath <- file.path('inst', 'extdata', paste0(bn, '.rds'))

    cat(txtPath, '->', rdsPath, '...\n')

    dat <- read.table(txtPath)[[1]]
    saveRDS(dat, file=rdsPath)

    cat('Ok.\n')
}

cat('Done.\n')

