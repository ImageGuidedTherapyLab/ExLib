#  Usage:
#   > source('plotDataSummary.R')
#library( ANTsR )
#library( randomForest )
#library( snowfall )
library( MASS )
#library( rlecuyer )
#rm(list=ls())

options("width"=200)
stopQuietly <- function(...)
  {
  blankMsg <- sprintf( "\r%s\r", paste( rep(" ", getOption( "width" ) - 1L ), collapse = " ") );
  stop( simpleError( blankMsg ) );
  } # stopQuietly()


# load data
modelData <-  read.csv( "./analysissummary.csv" )
modelData <-  subset(modelData,  !is.na(DCEAvg) )
modelData$EntropyPimo[ is.na(modelData$EntropyPimo) ] <- 0
modelData$HaralickPimo[is.na(modelData$HaralickPimo)] <- 0
modelData$Status = modelData$Status + 1




# lexical scope on color code and symbol code
# http://www.johndcook.com/R_language_for_programmers.html#scope
## put histograms on the diagonal
panel.summaryhist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 2.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col='cyan', ...)
    meanvalue <- mean(x)
    localdigits <- 4
    avgtxt <- paste("avg=", format( meanvalue , digits = localdigits ))
    stdtxt <- paste("std=", format( sd(x)     , digits = localdigits ))
    mintxt <- paste("min=", format( min(x)    , digits = localdigits ))
    maxtxt <- paste("max=", format( max(x)    , digits = localdigits ))
    text(meanvalue , 1.2, mintxt )
    text(meanvalue , 1.4, avgtxt )
    text(meanvalue , 1.6, stdtxt )
    text(meanvalue , 1.8, maxtxt )
    #if (do.legend) legend(43,2.5,c("GR","S2"), col = c(1,3), pch = c(1,3), bg="white")
    do.legend <<- FALSE
  
}

# lexical scope on color code and symbol code
# http://www.johndcook.com/R_language_for_programmers.html#scope
panel.summarylinear <- function(x, y)
{
  points(x[modelData$Status==1],y[modelData$Status==1],col=colcode[1] )
  points(x[modelData$Status==2],y[modelData$Status==2],col=colcode[2] )
  abline(lm(y~x), col='blue',lty='dashed')
  #TODO l1 regression
  #abline(rlm(x~y,abs(y),method = c("MM")), col='red')
       
  #abline(h=CrossOver,v=CrossOver,col='blue',lty=2)
}

# lexical scope on color code and symbol code
# http://www.johndcook.com/R_language_for_programmers.html#scope
panel.sorafenibsummarylinear <- function(x, y)
{
  #points(x,y,col=colcode,pch=symcode)
  points(   x[SorafenibData$idmethod=='manual' ],y[SorafenibData$idmethod=='manual' ] ,col="blue")
  points(   x[SorafenibData$idmethod=='rfmodel'],y[SorafenibData$idmethod=='rfmodel'] ,col="green")
  abline(lm(y[SorafenibData$idmethod=='manual' ]~x[SorafenibData$idmethod=='manual' ]),col='blue' ,lty='dashed')
  #abline(rlm(y[SorafenibData$idmethod=='manual' ]~x[SorafenibData$idmethod=='manual' ]),psi  = psi.huber(u, k = 0, deriv = 0),col='blue' )


  abline(lm(y[SorafenibData$idmethod=='rfmodel']~x[SorafenibData$idmethod=='rfmodel']),col='green',lty='dashed')
  #TODO l1 regression
  #abline(rlm(x~y,abs(y),method = c("MM")), col='red')
       
  #abline(h=CrossOver,v=CrossOver,col='blue',lty=2)
}

# lexical scope on color code and symbol code
# http://www.johndcook.com/R_language_for_programmers.html#scope
panel.tacesummarylinear <- function(x, y)
{
  #points(x,y,col=colcode,pch=symcode)
  points(   x[TACEData$idmethod=='manual' ],y[TACEData$idmethod=='manual' ] ,col="blue" ) 
  points(   x[TACEData$idmethod=='rfmodel'],y[TACEData$idmethod=='rfmodel'] ,col="green")
  abline(lm(y[TACEData$idmethod=='manual' ]~x[TACEData$idmethod=='manual' ]),col='blue' ,lty='dashed')
  abline(lm(y[TACEData$idmethod=='rfmodel']~x[TACEData$idmethod=='rfmodel']),col='green',lty='dashed')
  #TODO l1 regression
  #abline(rlm(x~y,abs(y),method = c("MM")), col='red')
       
  #abline(h=CrossOver,v=CrossOver,col='blue',lty=2)
}



# view first 10 lines
print(head(modelData ,n=10))

# print columns
print(names(modelData ))


#color code
numberOfUniqueLabels = 2
colcode <- character(numberOfUniqueLabels)
colcode[1] <- "red"
colcode[2] <- "green"
# lexical scope on color code and symbol code
# http://www.johndcook.com/R_language_for_programmers.html#scope
## put histograms on the diagonal
## put histograms on the diagonal
panel.summaryhist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}
# lexical scope on color code and symbol code
# http://www.johndcook.com/R_language_for_programmers.html#scope
## put histograms on the diagonal
panel.hist <- function(x, ...)
{
    ## Calculate and plot the two histograms
    usr <- par("usr"); on.exit(par(usr))

    maxden = 0
    for (iii in 1:numberOfUniqueLabels) {
       densitylist <- density(x[modelData$Status==(iii)])
       maxden = max(maxden,max(densitylist$y))
    }

    par(usr = c(usr[1:2], 0, maxden  ))
    for (iii in 1:numberOfUniqueLabels) {
       densitylist <- density(x[modelData$Status==(iii)])
       lines(densitylist,col=colcode[iii])
    }

}
## put (absolute) correlations on the upper panels,
## with size proportional to the correlations.
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y)
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste("cor=" , txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    #text(0.5, 0.5, txt, cex = cex.cor * r)
    text(0.5, 0.5, txt, cex = 0.6)
}
## put (absolute) correlations on the upper panels,
## with size proportional to the correlations.
panel.sorafenibcor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    rmanual  <- cor(x[SorafenibData$idmethod=='manual' ], y[SorafenibData$idmethod=='manual' ])
    rrfmodel <- cor(x[SorafenibData$idmethod=='rfmodel'], y[SorafenibData$idmethod=='rfmodel'])
    txtmanual  <- format(c(rmanual , 0.123456789), digits = digits)[1]
    txtrfmodel <- format(c(rrfmodel, 0.123456789), digits = digits)[1]
    if(rmanual  <  0.0){ txtmanual  <- paste("cor=-", txtmanual )}
    if(rmanual  >= 0.0){ txtmanual  <- paste("cor=" , txtmanual )}
    if(rrfmodel <  0.0){ txtrfmodel <- paste("cor=-", txtrfmodel)}
    if(rrfmodel >= 0.0){ txtrfmodel <- paste("cor=" , txtrfmodel)}
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txtmanual )
    #text(0.5, 0.5, txt, cex = cex.cor * r)
    text(0.5, 0.5, txtmanual  , cex = 1.0,col="blue" )
    text(0.5, 0.4, txtrfmodel , cex = 1.0,col="green")
}

## put (absolute) correlations on the upper panels,
## with size proportional to the correlations.
panel.tacecor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    rmanual  <- cor(x[TACEData$idmethod=='manual' ], y[TACEData$idmethod=='manual' ])
    rrfmodel <- cor(x[TACEData$idmethod=='rfmodel'], y[TACEData$idmethod=='rfmodel'])
    txtmanual  <- format(c(rmanual , 0.123456789), digits = digits)[1]
    txtrfmodel <- format(c(rrfmodel, 0.123456789), digits = digits)[1]
    if(rmanual  <  0.0){ txtmanual  <- paste("cor=-", txtmanual )}
    if(rmanual  >= 0.0){ txtmanual  <- paste("cor=" , txtmanual )}
    if(rrfmodel <  0.0){ txtrfmodel <- paste("cor=-", txtrfmodel)}
    if(rrfmodel >= 0.0){ txtrfmodel <- paste("cor=" , txtrfmodel)}
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txtmanual )
    #text(0.5, 0.5, txt, cex = cex.cor * r)
    text(0.5, 0.5, txtmanual  , cex = 1.0,col="blue" )
    text(0.5, 0.4, txtrfmodel , cex = 1.0,col="green")
}


modelDataControl = subset(modelData,  Status == 1  )
modelDataTreat   = subset(modelData,  Status == 2  )
print(head(modelDataControl ,n=10))
print(head(modelDataTreat   ,n=10))

# select image features
featureNames  = c("EntropyHE","HaralickHE","DistViableHE","DistNecrosisHE","HEPIMOOverlap","EntropyPimo","HaralickPimo","DistO2Pimo","DCEAvg","T2Abs","MedAir","OxyT2","T2Pct","T1Pre","T1Pst")

numtests= 2
treatcontrolpvalue     <- matrix( NA, nrow = numtests, ncol = length( featureNames ) )
colnames( treatcontrolpvalue    ) <- c( featureNames)
for (iii in 1:length(featureNames))
{
  print(featureNames[iii] )
  treatcontroltest     =      t.test( modelDataControl [[featureNames[iii]]], modelDataTreat   [[featureNames[iii]]])
  print( treatcontroltest    )
  treatcontrolpvalue[   1,iii]  =   2*pt( -abs(    treatcontroltest$statistic), df=    treatcontroltest$parameter,log=TRUE)
}
treatcontrolpvalue    =  treatcontrolpvalue[   ,order(treatcontrolpvalue[   1,])]

## # plot volume change 
## do.legend <- TRUE
pdf('DataSummary.pdf')
 pairs(~EntropyHE+HaralickHE+DistViableHE+DistNecrosisHE+HEPIMOOverlap+EntropyPimo+HaralickPimo+DistO2Pimo+DCEAvg+T2Abs+MedAir+OxyT2+T2Pct,
         data=modelData,
        diag.panel  = panel.hist, 
        lower.panel = panel.summarylinear,
        upper.panel = panel.cor
      )
dev.off()


