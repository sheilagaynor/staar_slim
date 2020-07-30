# Description: Generate a null model using the STAAR package using a wrapper of the GMMAT package.
# Inputs:
# results files : comma separated list of files containing output from STAAR (string of .Rds)

# Parse arguments
args <- commandArgs(T)
print(args)
# File inputs
result_files <- unlist(strsplit(args[1],","))
print(result_files)

result_out <- c()
for (i in 1:length(result_files)){
  result <- readRDS(result_files[i])
  result_out <- rbind(result_out, result)
}

# Save output
saveRDS(result_out, 'staar_out.Rds')
