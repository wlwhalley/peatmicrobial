# this code exports a higher resolution plot for use in publication etc
tiff("C:/Will/data/plots/plot.tiff", # the file extension where you want your plot saved
     units="in", #this is the units for the width. I use inches for some reason
     width=5, height=5, #this is the width and height of the plot in inches
     res=1200) # this is the resolution of the plot in DPI (dots per inch)
# insert ggplot code here for the full graph. The two lines should surround the graph. 
dev.off()

# Note: DPI is the plot resolution. 1200 is high. 600 is okay. 300 is blurry:
# the default DPI is 96!

# here is an example with a plot in it:
tiff("C:/Will/data/plots/plot.tiff", units="in", width=5, height=5, res=1200) 
ggplot(the_data, aes(x, y)) +
  geom_line() +
  xlab("stuff") +
  ylab("measured stuff")
dev.off()
