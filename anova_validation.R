library(readxl)

excel_file <- "C:/Users/aburtnerabt/Documents/K-State/Experimental Design/Assignment 3 Data.xlsx"

data <- read_excel(excel_file, sheet = "hay_example")

str(data)

hay_anova = aov(y ~ block + m*n*p*, data = data)
summary(hay_anova)