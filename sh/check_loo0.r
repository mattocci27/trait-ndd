# R values
output <- "loo.csv"
cat(paste("file", 
          "season",
          "hab",
          "trait",
          "n",
          "loo", sep = ","), file = output, append = FALSE, sep = "\n")
