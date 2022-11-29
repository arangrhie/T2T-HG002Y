#========== for the dot plot ==========

plotFilenameTemplate = "hg38Y_onto_hg002Y.dotplot.pdf"

addTitle = F
addLegend = F
delineateChromosome = F

#~~~~~ read the hg002 and hg38 lengths

scaffolds.hg002 = read_scaffold_lengths("chrY_hg002_v2.7.lengths")
xLength = sum(scaffolds.hg002$length)

scaffolds.hg38 = read_scaffold_lengths("hg38.chrY.lengths",sortByLength=F)
yLength = sum(scaffolds.hg38$length)

#~~~~~ read the dotplot

colClasses=c("character","integer","character","integer","integer")
dots = read.table("hg38Y_onto_hg002Y.dotplot.gz",header=T,comment.ch="",colClasses=colClasses)

#~~~~~ set up hg002 region boundaries
# these names and colors manually taken from Monika_HG002_classes.v5.dat

regionNames = c("PAR","X-DEG","XTR","AMPL","SAT","CEN","DYZ","HET","OTHER")

regionColor = c(rgb(151,203,153,max=255), # PAR
                rgb(255,239,87,max=255),  # X-DEG
                rgb(238,169,186,max=255), # XTR
                rgb(136,192,234,max=255), # AMPL
                rgb(106,71,0,max=255),    # SAT
                rgb(176,32,38,max=255),   # CEN
                rgb(106,71,0,max=255),    # DYZ
                rgb(119,119,119,max=255), # HET
                rgb(217,216,216,max=255)) # OTHER
names(regionColor) = regionNames

colClasses=c("character","numeric","numeric","character","numeric","character","numeric","numeric","character")
regions = read.table("Monika_HG002_classes.v5.dat",header=T,comment.ch="",colClasses=colClasses)

regions$class = gsub("^PAR.*$","PAR",regions$class)
regions$class = gsub("^XTR.*$","XTR",regions$class)
regions$class = gsub("^DYZ.*$","DYZ",regions$class)

legText  = regionNames
legColor = regionColor[legText]

#~~~~~ plot

xSpan = xLength
ySpan = yLength
xlim  = c(0,xLength)
ylim  = c(0,yLength)

windowSize = 12
aspect = 1.069 # this makes the plot area ≈ square if genomes are same length
aspect = aspect*ySpan/xSpan
width  = ifelse(aspect<=1,windowSize,windowSize/aspect);
height = width*aspect;
height = height*0.9909365    # special adjustment to match horizontal axis

plotFilename = NULL
if (!is.null(plotFilenameTemplate))
        plotFilename = plotFilenameTemplate

turnDeviceOff = F

if (is.null(plotFilename)) {
        quartz(width=width,height=height)
} else {
        print(paste("drawing to",plotFilename))
        pdf(file=plotFilename,width=width,height=height,pointsize=9)
        turnDeviceOff = T
        }

title = "lastz alignment, HG002 Y vs GRCh38 Y"
xlab  = "HG002 Y"
ylab  = "GRCh38 Y"

# create window
options(scipen=10)
par(mar=c(4,4.1,2.5,0.2)+0.1);    # BLTR
if (addTitle)
        {
        plot(NA,xlim=xlim,ylim=ylim,main=title,xlab=xlab,ylab=ylab)
} else {
        plot(NA,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab)
        }

# draw background regions
for (ix in 1:(nrow(regions)))
        {
        xLeft   = regions$start[ix]
        xRight  = regions$end[ix]
        color   = regionColor[regions$class[ix]]
        rect(xLeft,0,xRight,yLength,border=NA,col=color)
        }

# plot the dotplot
lines(dots$pos1,dots$pos2,lwd=1,col="black")

# delineate the chromosome
if (delineateChromosome)
        {
        lines(c(0,xLength),c(0,0),col="red")
        lines(c(0,xLength),c(yLength,yLength),col="red")
        lines(c(0,0),c(0,yLength),col="red")
        lines(c(xLength,xLength),c(0,yLength),col="red")
        }

if (addLegend)
        legend("bottomright",bg="white",cex=1.0,legText,fill=legColor)

if (turnDeviceOff) dev.off()




========== for the identity plot ==========

plotFilenameTemplate = "hg38Y_onto_hg002Y.hg002_covered.lastz.best{identity}.pdf"

addTitle = F

assembly1Name   = "HG002 Y"
assembly2Name   = "GRCh38 Y"
lengthsFilename = "chrY_hg002_v2.7.lengths"
bestsFilename   = "hg38Y_onto_hg002Y.hg002Y_best_identity.dat"
identityTag     = "identity"
yAnnotation     = "identity m/(m+mm)%"

#~~~~~ read the hg002 lengths

lengths = read.table(lengthsFilename,header=F)
nameToLength = lengths$V2
names(nameToLength) = lengths$V1
xLength = sum(nameToLength)

#~~~~~ read the best-identity table

colClasses=c("character","integer","integer","numeric")
bb = read.table(bestsFilename,header=F,comment.ch="",colClasses=colClasses)
colnames(bb) = c("chrom","start","end","identity")

#~~~~~ set up hg002 region boundaries
# these names and colors manually taken from Monika_HG002_classes.v5.dat

regionNames = c("PAR","X-DEG","XTR","AMPL","SAT","CEN","DYZ","HET","OTHER")

regionColor = c(rgb(151,203,153,max=255), # PAR
                rgb(255,239,87,max=255),  # X-DEG
                rgb(238,169,186,max=255), # XTR
                rgb(136,192,234,max=255), # AMPL
                rgb(106,71,0,max=255),    # SAT
                rgb(176,32,38,max=255),   # CEN
                rgb(106,71,0,max=255),    # DYZ
                rgb(119,119,119,max=255), # HET
                rgb(217,216,216,max=255)) # OTHER
names(regionColor) = regionNames

colClasses=c("character","numeric","numeric","character","numeric","character","numeric","numeric","character")
regions = read.table("Monika_HG002_classes.v5.dat",header=T,comment.ch="",colClasses=colClasses)

regions$class = gsub("^PAR.*$","PAR",regions$class)
regions$class = gsub("^XTR.*$","XTR",regions$class)
regions$class = gsub("^DYZ.*$","DYZ",regions$class)

#~~~~~ plot

width  = 12*1.017198               # to match horizontal scale of the dot plot
height = 5
turnDeviceOff = F

plotFilename = NULL
if (!is.null(plotFilenameTemplate))
        {
        plotFilename = plotFilenameTemplate
        plotFilename = gsub("[{]identity[}]","",plotFilename)
        plotFilename = gsub("[{]intervalName[}]",paste(".",identityTag,sep=""),plotFilename)
        }

if (is.null(plotFilename)) {
        quartz(width=width,height=height)
} else {
        print(paste("drawing to",plotFilename))
        pdf(file=plotFilename,width=width,height=height,pointsize=9)
        turnDeviceOff = T
        }

xMax = nameToLength["chrY_hg002"]

title = paste("highest alignment identities across ",assembly1Name,
              "\n(for lastz alignments of ",assembly2Name," to ",assembly1Name,")",sep="")
xlab = paste("position on ",assembly1Name,sep="")
ylab = yAnnotation
xlim = c(0,xMax)
ylim = c(-5,100)
gridColor = "blue"

# create window
options(scipen=10)
par(mar=c(4,4,2.4,1.7)+0.1)     # BLTR
if (addTitle)
        {
        plot(NA,xlim=xlim,ylim=ylim,main=title,xlab=xlab,ylab=ylab)
} else {
        plot(NA,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab)
        }

# draw background regions
for (ix in 1:(nrow(regions)))
        {
        xLeft   = regions$start[ix]
        xRight  = regions$end[ix]
        color   = regionColor[regions$class[ix]]
    rect(xLeft,0,xRight,100,border=NA,col=color)
        }

gridSpacing = 5*1000*1000
for (x in gridSpacing*(0:floor(xMax/gridSpacing)))
        lines(c(x,x),c(0,100),col=gridColor,lwd=0.25);
for (y in c(100,95,90,85,80))
        lines(c(0,xMax),c(y,y),col=gridColor,lwd=0.25);

#~~~~~ plot best-identity
bb.cs <- cityscape(bb[,2:4])
lines(bb.cs,type="l")

if (turnDeviceOff) dev.off()





========== support code ==========

read_scaffold_lengths <- function(lengthsFilename,scaffoldsOfInterest=NULL,sortByLength=T)
        #
        # Read a file containing scaffold names and lengths. Result is either
        # (a) reduced to scaffolds of interest, in order or (b) sorted by
        # decreasing length. A column is added, giving offsets for which the
        # scaffolds can be arranged along a number line.
        #
        # typical input:
        #       scaffold_100_arrow_ctg1 29829
        #       scaffold_101_arrow_ctg1 29808
        #       scaffold_102_arrow_ctg1 28782
        #       scaffold_103_arrow_ctg1 27584
        #       scaffold_104_arrow_ctg1 27050
        #    ...
        #
        # returns, e.g.
        #                        name    length    offset
        #            super_scaffold_1 215872496         0
        #            super_scaffold_2 165051286 215872496
        #       scaffold_2_arrow_ctg1 125510139 380923782
        #            super_scaffold_Z  84526827 506433921
        #       scaffold_5_arrow_ctg1  82438829 590960748
        #    ...
        #
        {
        scaffolds = read.table(lengthsFilename,header=F,colClasses=c("character","integer"))
        colnames(scaffolds) <- c("name","length")

        if (!is.null(scaffoldsOfInterest))
                {
                # pick ordered subset
                scaffolds = scaffolds[scaffolds$name %in% scaffoldsOfInterest,]
                scaffoldsToNumber = 1:length(scaffoldsOfInterest)
                names(scaffoldsToNumber) = scaffoldsOfInterest
                scaffolds = scaffolds[order(scaffoldsToNumber[scaffolds$name]),]
                }
        else if (sortByLength)
                {
                # sort by decreasing length
                scaffolds = scaffolds[order(-scaffolds$length),]
                }

        scaffolds[,"offset"] = rev(sum(as.numeric(scaffolds$length))-cumsum(rev(as.numeric(scaffolds$length))))

        scaffolds;
        }

cityscape <- function(df)
        {
        # expects df to be 2-dimensional with three columns
        #  col 1 is start of interval, col 2 is end (inclusive)
        #  col3 is the 'value' over the interval
        # result will have two columns
        #  col1 is an interval endpoint       (NA for gap between intervals)
        #  col2 is the value at that endpoint (NA in gap)

        # create an array with two rows for every interval and 1 for every gap
        rows <- 2*dim(df)[1] + gapcount(df)
        cs   <- array(NA,c(rows,2))

        # fill in the array

        nextStart <- df[1,1]

        csRow <- 0
        for (dfRow in  1:dim(df)[1])
                {
                if (df[dfRow,1]!=nextStart)
                        {
                        csRow <- csRow + 1      # leave a gap
                        }
                csRow <- csRow + 1
                cs[csRow,1] <- df[dfRow,1] # copy interval start
                cs[csRow,2] <- df[dfRow,3]
                csRow <- csRow + 1
                cs[csRow,1] <- df[dfRow,2] # copy interval end
                cs[csRow,2] <- df[dfRow,3]
                nextStart   <- df[dfRow,2] + 1
                }

        cs
        }
