
loadfunctions = function(
    projectname,
    directorypattern=NULL,
    functionname=NULL,
    RcodeDirectory="R",
    keydirectories=c("r", "\\.r", "\\_r", "rfunctions", "\\.rfunctions", "\\_rfunctions" ),
    toignore = c("retired", "_archive", "archive", "orphan", "request", "example" ),
    filepattern="\\.r$",
    directory=NULL ) {

    # used to load local functions conveniently
    # sequence slightly important ... modify with care
    # the  '\\' used in the above are escape sequences used with regex

    project.codedirectory = function(...) {
        ## this function is required to bootstrap the other project level functions
        if (!exists("code_root")) {
        if (is.null(directory)) {
            stop("Define 'code_root' in your Rprofile...")  
        }
        code_root = directory
        }
        sep = .Platform$file.sep
        dirs = paste0( c(...) , collapse=sep )
        pd = file.path( code_root, dirs )
        pd = gsub( paste( sep, "$", sep=""), "", pd) # strip last path element
        if (! file.exists( pd)) {
        stop (paste("Directory", pd, "not found. "))
        }
        return ( pd )
    }

    LoadFiles = function( ... ) {
        filelist = c(...)
        if ( length(filelist) > 0  ) for ( nm in c(...) ) source( file=nm )
        return( filelist )
    }



    fs = .Platform$file.sep
    keydirectories = paste( "\\", fs, keydirectories, "\\", fs, sep="") # make sure only match directories and not file names
    filestosource = NULL

    for (pn in projectname ) {

        projectdirectory = project.codedirectory( directory, pn )

        for (searchdirectory in c( file.path( projectdirectory, RcodeDirectory ), projectdirectory ) ) {  # first try in RcodeDirectory and then the project if not found in first pass

        projectfiles = NULL
        projectfiles = list.files( path=searchdirectory, pattern=filepattern, full.names=T, recursive=T,  ignore.case=T, include.dirs=F )

        # filter on directory pattern
        if (!is.null(directorypattern)) {
            keep = grep ( directorypattern, projectfiles, ignore.case =T )
            if (length(keep) == 0) {
            # not found, try an approximate match
            keep = agrep( directorypattern, projectfiles, ignore.case =T )
            }
            if (length(keep)>0) projectfiles = projectfiles[ keep ]
        }

        # functionname filter
        if (!is.null(functionname)) {
            keep = grep ( functionname, projectfiles, ignore.case =T )
            if (length(keep) == 0) {
            # if no match try approximate matching
            keep = agrep( functionname, projectfiles, ignore.case =T )
            }
            if (length(keep)>0) {
            for ( nm in keep ) source( file= projectfiles[nm] )
            return( projectfiles[keep] )  # this breaks out of the loops as it is a single file-search and source call
            }
        }


        # remove archived functions, etc. that are to be "ignored"
        rem  = NULL
        for (i in toignore) {
            rem0 = grep ( i, projectfiles,  ignore.case =T )
            if (length( rem0)>0 ) rem = c( rem, rem0)
        }
        if ( length(rem)>0 ) {
            rem = unique(rem)
            projectfiles = projectfiles[-rem]
        }

        # filter on key directories
        # determine correct directory .. first look for a set of key directories, if not then scan all files
        md = NULL
        for ( kd in keydirectories) {
            md = c( md, grep ( kd, projectfiles, ignore.case =T ) )
            md = unique( md )
        }
        if (length(md) > 0 ) {
            projectfiles = projectfiles[ md ]
            filestosource = c( filestosource, projectfiles )
            break()  # no need to continue with the loop
        }
        }

    }

    LoadFiles( filestosource )

    return( filestosource )

}




