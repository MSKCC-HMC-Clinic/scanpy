#!/usr/bin/env Rscript

# assume this is run using 'test_arguments()' in scanpy/tests/external/test_execute_r_function.py
dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)  # create personal library
.libPaths(Sys.getenv("R_LIBS_USER"))  # add to the path

args = commandArgs(trailingOnly=TRUE)

# default arguments in the python subprocess call in execute_r_script include:
# args[1]:  '--vanilla'
# args[2]:  'cache'
args1 = args[1]
cache_dir_path = args[2]

# optional arguments are strings that could be filenames, input parameters, etc 
filename = args[3] # 'test.txt': this will be the filename
content = args[4] # 'contents': this will be the text file's contents

# create .txt file under /_tmp directory and write args[4] to file
write_path = file.path(cache_dir_path, filename)
write_file = file(write_path)
writeLines(content, write_file)
close(write_file)
