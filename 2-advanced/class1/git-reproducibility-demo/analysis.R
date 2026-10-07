# Gene expression analysis
# Version 1

data <- read.csv("expression.csv")

# Calculate fold change
data$fold_change <- data$Treatment / data$Control

# Identify genes with at least 2-fold increase
results <- subset(data, fold_change >= 2)

# Save results
dir.create("out", showWarnings = FALSE)

write.table(results, "out/results.txt", sep="\t", row.names = FALSE)

print(results)
