#======================================================
#this function plots shaded polygons, even when there
#are NAs in the sequence.
#this ideal for representing SEs or SDs around means in an
#timeseries with gaps
#======================================================
#NAs in lowerCI and upperCI need to match
shade.confidence.interval <- function(xvals,lowerCI,upperCI, band.col='dark grey', border.col=NA,transparent=0.5)
{
  i <- min(which(!is.na(lowerCI)))
  inds <- NULL
  while (i <= length(xvals))
  {   #cat(i,"\n")
      if ((is.na(lowerCI[i])) | (i==length(xvals)))
        {
            if (length(inds)>1)
            {
              if (i==length(xvals)){inds <- c(inds,i)}
              #cat("about to plot ",i)
              polygon(x=c(xvals[inds],xvals[inds[order(inds,decreasing=T)]]),
                      y=c(lowerCI[inds],upperCI[inds[order(inds,decreasing=T)]]),col=adjustcolor(band.col, alpha.f=transparent),border=border.col)
            inds <- NULL
            }
        }else{
            inds <- c(inds,i)
        }


      i <- i+1
  }
}

