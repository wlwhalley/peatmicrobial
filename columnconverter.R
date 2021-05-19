# this script was written by Will Burn wlwhalley@gmail.com or wlb501@york.ac.uk
# lots of problems occur when columns in your dataframe are the wrong class
# this is especially common in ecological data because R gets confused by NA
# before analysis, check all of your columns using
str(yourdata)
## in R you can select individual columns like this
yourdata$yourcolumn
# categorical factors should read out for example:
# natcode     : Factor w/ 5 levels "E","NP","P","S",
# converting individual columns to a factor is easy 
yourdata$categorical <- as.factor(yourdata$categorical)
# numeric variables are harder, they should look like
# doc: num  20.4 26.5 22.5 34.8 34.4 
# but converting them is a bit different for reasons I won't go into
# they need to be converted to a character first
yourdata$numeric <- as.numeric(as.character(yourdata$numeric))
# sometimes you may wish to bulk converter columns, for example
# you may want to make all integers to numeric
# you can use this package and code
### converts all columns of a specific class to numeric
install.packages("remotes")
remotes::install_github("gsk3/taRifx.geo")
library(taRifx)
yourdata <- japply(yourdata, which(sapply(yourdata, class)=="character"), as.numeric) # this line takes all the columns in your dataset that are characters and converts them to numeric. 
# always check again after doing so
str(yourdata)