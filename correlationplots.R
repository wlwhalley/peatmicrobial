########
# this code written by Will Burn wlb501@york.ac.uk or wlwhalley@gmail.com
# for more code visit: https://github.com/wlwhalley/peatmicrobial
# first we install all the packages required:
packages <- c("devtools","gridExtra", "remotes", "car","psych", "data.table","scales","dendextend","Hmisc", "rlang", "plyr", "Hmisc", "ggplot2", "dplyr", "tidyr", "vegan", "DT", "reshape2", "knitr", "lubridate", "pwr", "psy", "car", "doBy", "corrplot", "RcmdrMisc", "questionr", "vcd", "multcomp", "KappaGUI", "rcompanion", "gridExtra", "factoextra", "corrplot", "FSA", "MASS", "scales", "nlme", "psych", "ordinal", "lmtest", "ggpubr", "dslabs", "stringr", "assist", "ggstatsplot", "forcats", "styler", "remedy", "addinslist", "esquisse", "here", "summarytools", "magrittr", "tidyverse", "funModeling", "pander", "cluster", "abind")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])}
invisible(lapply(packages, library, character.only = TRUE)) # this bit takes all the packages in the above list, if they arent installed it installs them, if they are already it ust loads them
# first import your data
yourdata <- read.csv("C:/yourfolders/somefolders/yourdata", row.names=1)
str(yourdata) # check the structure of your data
yourdata$natcode <- as.factor(yourdata$natcode) # in my data, a lot of columns get imported as characters, this changes them to a factor
install.packages("remotes")
remotes::install_github("gsk3/taRifx.geo")
library(taRifx)
yourdata <- japply(yourdata, which(sapply(yourdata, class)=="character"), as.numeric) # this takes any column that is a character and converts it to a numeric variable
str(yourdata) # check the structure of your data again after transofrmations
# now centre and z score all your columns of interest
azenwater$lat<-scale(as.numeric(yourdata$lat), center=TRUE,scale=TRUE)
yourdata$long<-scale(as.numeric(yourdata$long), center=TRUE,scale=TRUE)
yourdata$slope<-scale(as.numeric(yourdata$slope), center=TRUE,scale=TRUE)
yourdata$aspect<-scale(as.numeric(yourdata$aspect), center=TRUE,scale=TRUE)
yourdata$elevationcm<-scale(as.numeric(yourdata$elevationcm), center=TRUE,scale=TRUE)
yourdata$peatdepth<-scale(as.numeric(yourdata$peatdepth), center=TRUE,scale=TRUE)
# etc etc - you may need to do many lines of this. Remove junk columns or ones that don't need to be correlated
res1 <- cor.mtest(as.matrix(yourdata[6:38]), conf.level = .95, na.rm = TRUE) # after yourdata, you need to put in the numbers of the columns you want to correlate
# for example, here I want to correlate columns 6 to 38
# to find the number of a specific column use which( colnames(yourdata)=="hazen" ) 
pmat<-as.matrix(res1$p)
lowCI<-as.matrix(res1$lowCI)
uppCI<-as.matrix(res1$uppCI)
cortest<-rcorr(as.matrix(yourdata[6:38])) # remember to adjust the numbers here
plot <- corrplot(as.matrix(cortest$r),p.mat=pmat,low = lowCI, upp = uppCI,
         rect.col = "navy", plotC = "rect", cl.pos = "n", sig.level=0.05)
plot
# now you can export your plot as a hi res figure (or a pdf) 
tiff("C:/somefiles/folder/imagename.tiff", units="in", width=5, height=5, res=1200) # units = inches here, width and height in inches, res is the resolution of the image in DPI
plot
dev.off()














