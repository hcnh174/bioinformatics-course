# Gene expression analysis
# Version 2

data <- read.csv("expression.csv")

# Calculate fold change
data$fold_change <- data$Treatment / data$Control

# Identify genes with at least 2-fold increase or decrease
#results <- subset(data, fold_change >= 2)
results <- subset(data, fold_change >= 2 | fold_change <= 0.5)

# Save results
dir.create("out", showWarnings = FALSE)

write.table(results, "out/results.txt", sep="\t", row.names = FALSE)

print(results)
