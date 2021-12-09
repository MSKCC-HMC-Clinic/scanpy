#!/usr/bin/env Rscript

dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)  # create personal library
.libPaths(Sys.getenv("R_LIBS_USER"))  # add to the path

# test taking in and reading different argument types
args = commandArgs(trailingOnly=TRUE)
print(args)

# test if there is exactly 2 arguments: if not, return an error
# if (length(args)!=2) {
#   stop("Missing arguments", call.=FALSE)
# }

# args[1]: --vanilla

# args[1]: try to read in a filename and open it as a csv
filename = args[2]
df = read.csv(file=filename, header=FALSE)
print(df)

# args[2]: try to read in a list in string format, ie "a,b,c" and parse it
#  as list("a", "b", "c")

string_list = args[3]
list = as.list(strsplit(string_list, ",")[[1]])
test_list = list("a", "b", "c")
print(identical(test_list,list))

# args[3]: try to read in AnnData converted to robject??