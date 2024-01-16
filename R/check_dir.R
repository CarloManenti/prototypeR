check_dir <- function(directory2check.dir){
    ### Description ####
    # Checks if the directory exists;
    # if not, it creates the directory

    # example usage
    # check_dir('~/Documents/PhD')



    if(!dir.exists(directory2check.dir)){
        message(paste0('--- new directory defined at : ',
                       directory2check.dir , '---\n'))
        dir.create(directory2check.dir)
    }
}
