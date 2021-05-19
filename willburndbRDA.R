# this script was written by William Burn as an adaption from a script by Phil Brailey (!)
# wlwhalley@gmail.com or wlb501@york.ac.uk
##############################################################
# this is the code I use to install packages
packages <- c("plyr", "Hmisc", "ggplot2", "speedyseq", "dplyr", "tidyr", "vegan", "DT", "reshape2", "knitr", "lubridate", "pwr", "psy", "car", "doBy", "corrplot", "RcmdrMisc", "questionr", "vcd", "multcomp", "KappaGUI", "rcompanion", "gridExtra", "factoextra", "corrplot", "FSA", "MASS", "scales", "nlme", "psych", "ordinal", "lmtest", "ggpubr", "dslabs", "stringr", "assist", "ggstatsplot", "forcats", "styler", "remedy", "snakecaser", "addinslist", "esquisse", "here", "summarytools", "magrittr", "tidyverse", "funModeling", "pander", "cluster", "abind")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])}
# load the packages
invisible(lapply(packages, library, character.only = TRUE))
##############################################################
#distance based redundancy analysis
# first create a correlation matrix of measures
# see if values are sig related to each other
# and sense check whether to include them in the analysis
# first ensure all variables are *numeric* and convert to Z-scores
metadata$lat<-scale(as.numeric(metadata$lat), center=TRUE,scale=TRUE) # you will need to change all of these to your variables! 
metadata$long<-scale(as.numeric(metadata$long), center=TRUE,scale=TRUE)
metadata$slope<-scale(as.numeric(metadata$slope), center=TRUE,scale=TRUE)
metadata$aspect<-scale(as.numeric(metadata$aspect), center=TRUE,scale=TRUE)
metadata$elevationcm<-scale(as.numeric(metadata$elevationcm), center=TRUE,scale=TRUE)
metadata$peatdepth<-scale(as.numeric(metadata$peatdepth), center=TRUE,scale=TRUE)
metadata$calluna<-scale(as.numeric(metadata$calluna), center=TRUE,scale=TRUE)
metadata$moss<-scale(as.numeric(metadata$moss), center=TRUE,scale=TRUE)
metadata$sphag<-scale(as.numeric(metadata$sphag), center=TRUE,scale=TRUE)
metadata$sedge<-scale(as.numeric(metadata$sedge), center=TRUE,scale=TRUE)
metadata$grass<-scale(as.numeric(metadata$grass), center=TRUE,scale=TRUE)
metadata$bare<-scale(as.numeric(metadata$bare), center=TRUE,scale=TRUE)
metadata$doc<-scale(as.numeric(metadata$doc), center=TRUE,scale=TRUE)
metadata$n<-scale(as.numeric(metadata$n), center=TRUE,scale=TRUE)
metadata$X254<-scale(as.numeric(metadata$X254), center=TRUE,scale=TRUE)
metadata$X400<-scale(as.numeric(metadata$X400), center=TRUE,scale=TRUE)
metadata$X465<-scale(as.numeric(metadata$X465), center=TRUE,scale=TRUE)
metadata$X665<-scale(as.numeric(metadata$X665), center=TRUE,scale=TRUE)
metadata$suva<-scale(as.numeric(metadata$suva), center=TRUE,scale=TRUE)
metadata$moisturem<-scale(as.numeric(metadata$moisturem), center=TRUE,scale=TRUE)
metadata$hazen<-scale(as.numeric(metadata$hazen), center=TRUE,scale=TRUE)
metadata$e4e6<-scale(as.numeric(metadata$e4e6), center=TRUE,scale=TRUE)
metadata$dnamgg<-scale(as.numeric(metadata$dnamgg), center=TRUE,scale=TRUE)
metadata$abbymoisture<-scale(as.numeric(metadata$abbymoisture), center=TRUE,scale=TRUE)
metadata$bdgcm3<-scale(as.numeric(metadata$bdgcm3), center=TRUE,scale=TRUE)
metadata$lightmolm2<-scale(as.numeric(metadata$lightmolm2), center=TRUE,scale=TRUE)
metadata$rainfallmm<-scale(as.numeric(metadata$rainfallmm), center=TRUE,scale=TRUE)
metadata$tair<-scale(as.numeric(metadata$tair), center=TRUE,scale=TRUE)
metadata$tsoil<-scale(as.numeric(metadata$tsoil), center=TRUE,scale=TRUE)
metadata$lightcol<-scale(as.numeric(metadata$lightcol), center=TRUE,scale=TRUE)
metadata$rainfallcol<-scale(as.numeric(metadata$rainfallcol), center=TRUE,scale=TRUE)
metadata$taircol<-scale(as.numeric(metadata$taircol), center=TRUE,scale=TRUE)
metadata$tsoilcol<-scale(as.numeric(metadata$tsoilcol), center=TRUE,scale=TRUE)
metadata$ph<-scale(as.numeric(metadata$ph), center=TRUE,scale=TRUE)
# do not include variables you don't have or that are blank


##########################################################
# this creates a correlation matrix of variables so you can examine which ones are related
res1 <- cor.mtest(as.matrix(metadata[6:30]), # you will need to change this to the *numeric* variables in your dataset e.g if you data is called "data" and you want columns 6:9 use as.matrix(data[6:9]) here
                  conf.level = .95, na.rm = TRUE) # conf level sets the significance for a correlation at 0.05
pmat<-as.matrix(res1$p)
lowCI<-as.matrix(res1$lowCI) 
uppCI<-as.matrix(res1$uppCI)
cortest<-rcorr(as.matrix(metadata[6:30])) # does a correlation test on the data
cplot <- corrplot(as.matrix(cortest$r),p.mat=pmat,low = lowCI, upp = uppCI,
         rect.col = "navy", plotC = "rect", cl.pos = "n", sig.level=0.05) # makes a correlation plot. Sometimes if you have many variables these are quite hard to read
cplot
##########################################################
# export your plots as high resolution .tiff so you can zoom in one them! 
tiff("C:/Will/data/plots/filename.tiff", units="in", width=5, height=5, res=1200) # this can also do other files e.g pdf or jpeg but I like tiff. The width and height are the height of the image in inches - an A3 graph would be 11 x 16 inch ish. 
cplot
dev.off() # you can change the resolution of your exported image with the "res" above it is measured in DPI - 1200 is a very high quality graph. 600DPI is suitable for publications, 900DPI is what Nat Geo would want for a magazine

##########################################################
# do the dbRDA

#Create environment file of subsetted variables
env <- as.data.frame(cbind(metadata$ph, metadata$lat, metadata$long, metadata$slope, metadata$aspect, metadata$elevationcm, metadata$peatdepth, metadata$calluna, metadata$moss, metadata$sphag, metadata$sedge, metadata$grass, metadata$bare, metadata$doc, metadata$n, metadata$suva, metadata$hazen, metadata$moisturem, metadata$e4e6, metadata$dnamgg, metadata$abbymoisture, metadata$bdgcm3))
colnames(env)<-c("ph", "lat","long","slope","aspect","elevationcm","peatdepth", "calluna","moss","sphag","sedge","grass","bare","doc","n","suva","hazen","moisturem","e4e6","dnamgg","abbymoisture","bdgcm3")
# make a distance based RDA model

######## MAKE YOUR DIST OBJECT USING VEGAN HERE ######### 
dist <- vegdist(x, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE, ...)
##########
dbRDA1<-capscale(dist ~ lat+long+slope+aspect+elevationcm+peatdepth+ph+calluna+moss+sphag+grass+sedge+bare+doc+n+suva+hazen+moisturem+e4e6+dnamgg+bdgcm3, data=metadata, na.action = "na.omit")
summary(dbRDA1)
# vif checks model coefficients, if it is above 5 for a variable it is not ideal, above 10 is very not ideal - it means they are probably autocorrelating values. 
vif(dbRDA1)
#Test the overall model
test1<-anova(dbRDA1, perm = 9999)## this tests the significance of your dbRDA model.
test1
#Test which variables are significant
test2<-anova(dbRDA1, by="terms", permu=9999)
test2
#Test which axes are significant
test3<-anova(dbRDA1, by="axis", permu=9999)
test3
plot(dbRDA1) # gives a very quick plot of the dbRDA with all the variables (even the ones that are non-significant.)
### Plot a dbRDA
metadata2<-na.omit(metadata) #remove row with missing variable - might not be relevant for your work but it was for mine 
smry <- summary(dbRDA1)
df1  <- data.frame(smry$sites[,1:2])       # PC1 and PC2
df1$manage<-metadata2$manage # you can change which variable your dots are grouped by here
df1 <- df1[order(df1$manage),] 
find_hull <- function(df1) df1[chull(df1$CAP1, df1$CAP2), ]
hulls <- ddply(df1, "manage", find_hull)
cent<-aggregate(cbind(df1$CAP1,df1$CAP2) ~ manage, data = df1, FUN = mean)
segs<-merge(df1, setNames(cent, c('manage', 'V1','V2')), 
            by = 'manage', sort = TRUE)
df1$seg1<-segs$V1
df1$seg2<-segs$V2
df2  <- data.frame(smry$biplot[,1:2])     # loadings for PC1 and PC2
#subset df2 with only significant variables
test2<-na.omit(test2)
df2$p<-as.numeric(test2$`Pr(>F)`)
df2<-subset(df2, p < 0.05)
#Plot with sig variables and no ellipses / hulls etc
rda.plot <- ggplot(df1, aes(x=CAP1, y=CAP2, color = metadata2$manage)) + 
  geom_point(size=3) +
  scale_color_manual(values=c("red", "green", "orange","yellow", "blue", "pink",
                              "#AD6F3B", "purple"), name="Management",
                     #breaks=c("E", "P", "S", "Y"), labels=c("Exmoor","Peak District", "Scotland", "Yorkshire (Grouse Moor)")) + # I use this bit of code when grouping by national site an just swap them out
                     breaks=c("F", "M", "U", "MH", "10Y", "5R", "D", "I"), labels=c("Rotationally Burnt (Grouse Moor)","Mown (Grouse Moor)", "Unmanaged (uncut Grouse Moor)", "80y post burn (ex-Grouse Moor)","10y post restoration", "5y post restoration", "Degraded",  "Intact")) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  theme_classic() # I liked the classic theme, but you can change it here. There's a 1998 Excel theme somewhere if you want to feel retro. 
rda.plot
rda.biplot <- rda.plot +
  geom_segment(data=df2, aes(x=0, xend=CAP1, y=0, yend=CAP2), 
               color="black", arrow=arrow(length=unit(0.01,"npc"))) + # this bit is where you can change size and colour of the variable arrows. 
  geom_text(data=df2, 
            aes(x=CAP1,y=CAP2,label=rownames(df2),
                hjust=0.5*(1-sign(CAP1)),vjust=0.5*(1-sign(CAP2))), 
            color="black", size=4.5) # change the size, colour of the variable text here
rda.biplot 

##### export it as a nice image
tiff("C:/Will/data/plots/filename.tiff", units="in", width=5, height=5, res=1200) 
rda.biplot
dev.off() 