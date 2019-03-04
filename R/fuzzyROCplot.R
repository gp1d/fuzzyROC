fuzzyROCplot = function(fR,
                        mainTitle="",
                        set=TRUE, # if TRUE plots the upper and lower curves and shades the aea between them
                        segm=FALSE # if TRUE plots the segments joing starred and daggered points
                        )
  # fR is the output of the fuzzyROC function
{
  plot(c(0,1),c(0,1),
       type="n",
       xlim=c(0,1),ylim=c(0,1),
       main=mainTitle,
       xlab="1-Specificity",ylab="Sensitivity")
  points(c(1-fR$F0star,rev(1-fR$F0dagger)),
         c(1-fR$F1star,rev(1-fR$F1dagger)),
         col=gray(.85),
         lwd=6)
  if(set==TRUE){
  polygon(c(1,1-fR$F0star,0,rev(1-fR$F0dagger)),
          c(1,1-fR$F1star,0,rev(1-fR$F1dagger)),
          col=gray(0.75),
          lwd=0.001)
  }
  if (segm==TRUE) {
    for (jj in 1:length(fR$cL)) {
      seg.col = ifelse(set,gray(0.85),gray(.65))
      seg.lwd = ifelse(set,2,4)
      segments(1-fR$F0star[jj],1-fR$F1star[jj],1-fR$F0dagger[jj],1-fR$F1dagger[jj],
               col=seg.col,lwd=seg.lwd)
  }
  abline(0,1)
}
}
