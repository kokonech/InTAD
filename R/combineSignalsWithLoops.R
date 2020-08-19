
#' Preparation for correlation analysis via loops
#'
#' This function combines signals and genes based on the usage
#' of loops obtained from HiC data analysis
#' @param object InTADSig object
#' @param loopsInitDf Data frame with loops. By default 6-column format
#' \emph{(chr1,start1,end1,chr2,start2,pos2)} is expected.
#' @param fragmentLength In case the input format is 4-column \emph{(chr1,middlePos1,
#' chr2, middlePos2)} fragment length should be provided to extend the
#' corresponding loci for loop start and end positions.
#' @param tssWidth The transcription start site width is used to control overlaps
#' with loop anchor. Default is 2000 base pairs.
#' @param extSize The loop endings can be extended upstream and downstream
#' with provided corresponding increase size in base pairs.
#' @return Updated InTADSig object containing genes connected to signals
#' via loops
#'
#' @details
#' The expected input is the loops data.frame applied to find
#' connections of signals to genes. This data.frame could be
#' in two formats: either \emph{(chr1,start1,end1,chr2,start2,end2)} or
#' \emph{(chr1,middlePos1,chr2,middlePos2)} with fragment size.
#'
#' @importMethodsFrom GenomicRanges intersect
#' @importFrom BiocGenerics start end
#' @importFrom S4Vectors queryHits subjectHits
#'
#' @export
#'
combineWithLoops <- function( object, loopsInitDf,
                              fragmentLength = 0,
                              tssWidth = 2000,
                              extSize = 0) {

    if (!is(object, "InTADSig"))
        stop("Object must be an InTADSig!")

    if (fragmentLength == 0) {
      message("NOTE: 6-column loops format is assumed.")
      if (ncol(loopsInitDf) < 6) {
        stop("Loops data.frame should have at leaset 6 columns!")
      }
      selLoopsDf <-loopsInitDf[,1:6]
    } else {
      message("NOTE: binSize is provded, 4-column loops format is assumed.")
      if (ncol(loopsInitDf) < 4) {
        stop("Loops data.frame should have at least 6 columns!")
      }
      selLoopsDf <- data.frame( chr1 = loopsInitDf[,1],
                                x1 = loopsInitDf[,2] - fragmentLength/2,
                                y1 = loopsInitDf[,2] + fragmentLength/2,
                                chr2 = loopsInitDf[,3],
                                x2 = loopsInitDf[,4] - fragmentLength/2,
                                y2 = loopsInitDf[,4] + fragmentLength/2)
    }
    # clear enhancers with no connection
    loops <- paste0("loop",1:nrow(selLoopsDf))
    rownames(selLoopsDf) <- loops
    object@loopsDf <- selLoopsDf

    # form regions
    regions1 <- GRanges(seqnames = selLoopsDf[,1],
                        ranges = IRanges(start = selLoopsDf[,2], end = selLoopsDf[,3]))
    regions2 <- GRanges(seqnames = selLoopsDf[,4],
                        ranges = IRanges(start = selLoopsDf[,5], end = selLoopsDf[,6]))

    if (extSize > 0) {
      regions1 <- shift(regions1, -extSize)
      regions1 <- resize(regions1, width(regions1) + extSize*2)
      regions2 <- shift(regions2, -extSize)
      regions2 <- resize(regions2, width(regions2) + extSize*2)
    }
    geneGR <- rowRanges(object@sigMAE[["exprs"]])
    sigGR <- rowRanges(object@sigMAE[["signals"]])

    tss <- GRanges( seqnames=as.character(seqnames(geneGR)),
                    IRanges(start=ifelse(as.character(strand(geneGR))=="+",
                                         start(geneGR),end(geneGR)), width=tssWidth,
                            names=values(geneGR)$gene_name),
                    strand=strand(geneGR), geneid=values(geneGR)$gene_id )
    tss <- shift(tss, -tssWidth/2)

    # overlaps with enhancers
    enhUpstream <- findOverlaps(regions1,sigGR)
    enhDownstream <- findOverlaps(regions2,sigGR)

    # overlap with genes
    genesUpstream <- findOverlaps(regions1,tss)
    genesDownstream <- findOverlaps(regions2,tss)

    enhUpDf <- data.frame(Loop=loops[queryHits(enhUpstream)],
                          Enhancer=as.character(sigGR[subjectHits(enhUpstream)]),
                          stringsAsFactors = F)
    genesDownDf <- data.frame(Loop=loops[queryHits(genesDownstream)],
                              GeneId=tss[subjectHits(genesDownstream)]$geneid,
                              GeneName=names(tss[subjectHits(genesDownstream)]),
                              stringsAsFactors = F)
    res1 <- merge(enhUpDf, genesDownDf, by = "Loop")


    enhDownDf <- data.frame(Loop=loops[queryHits(enhDownstream)],
                            Enhancer=as.character(sigGR[subjectHits(enhDownstream)]),
                            stringsAsFactors = F)
    genesUpDf <- data.frame(Loop=loops[queryHits(genesUpstream)],
                            GeneId=tss[subjectHits(genesUpstream)]$geneid,
                            GeneName=names(tss[subjectHits(genesUpstream)]),
                            stringsAsFactors = F)
    res2 <- merge(enhDownDf, genesUpDf, by = "Loop")

    combRes <- rbind(res1,res2)
    object@loopConnections <- split(combRes, 1:nrow(combRes))

    message(paste("Combined",sum(sapply(object@loopConnections, nrow)),
                    "signal-gene pairs with loops" ))
    return(object)

}


findLoopPairCorrelation <- function(inData, signalVals, countVals, corMethod) {

  countSel <- as.numeric( countVals[inData$GeneId,])
  if ( sum(countSel) == 0) {
    return(NULL)
  }
  sigSel = as.numeric(signalVals[inData$Enhancer,])
  cors <- suppressWarnings(cor.test( sigSel, countSel, method=corMethod))
  euDis = as.numeric(dist(rbind ( sigSel, countSel)))
  corDF <- data.frame(peak=inData$Enhancer,loop=inData$Loop,
                      gene=inData$GeneId, name=inData$GeneName,
                      cor=cors$estimate,
                      pvalue=cors$p.value,
                      eucDist= euDis,
                      stringsAsFactors =FALSE)
  return(corDF)
}


#' Function to perfrom correlation analysis via loops.
#'
#' This function combines genes and signals using obtained loop connections.
#' @param object InTADSig object with signals and genes combined via loops
#' @param method Correlation method: "pearson" (default), "kendall", "spearman"
#' @param adj.pval Perform p-value adjsutment and include q-values in result
#' @return A table with correlation values for signal-gene pairs including
#' correlation p-value and euclidian distance.
#' @importFrom stats na.omit cor.test dist
#' @import qvalue
#' @export
#'
findCorFromLoops <- function(object, method = "pearson", adj.pval = FALSE ) {

  if (!is(object, "InTADSig"))
    stop("Object must be an InTADSig!")

  if (length(object@loopConnections)  == 0)
    stop("No signals and genes are connected via loops!
         Use combineWithLoops()  function.")

  combList <- object@loopConnections
  selLoopsDf <- object@loopsDf

  if (object@ncore > 1 && requireNamespace("parallel")) {
    message("Running in parallel...")
    allRes <- parallel::mclapply(combList, findLoopPairCorrelation,
                                 assay(object@sigMAE[["signals"]]),
                                 assay(object@sigMAE[["exprs"]]),
                                 method)
  } else {
    allRes <- lapply(combList, findLoopPairCorrelation,
                     assay(object@sigMAE[["signals"]]),
                     assay(object@sigMAE[["exprs"]]),
                     method
    )
  }

  allRes <- allRes[lengths(allRes) != 0]
  allCortab <- do.call(rbind, lapply(allRes, function(x) x))

  loopAnn <- selLoopsDf[ allCortab$loop,]
  allCortab$loopStart = paste0(loopAnn[,1],":",loopAnn[,2],"-",loopAnn[,3])
  allCortab$loopEnd = paste0(loopAnn[,4],":",loopAnn[,5],"-",loopAnn[,6])

  # re-order, discard IDs of loops
  allCortab <- allCortab[, c(1,8,9,3:7)]

  if (adj.pval) {
    qvals <- tryCatch({
      qvalue(allCortab$pvalue)
    }, error = function(e) {
      qvalue(allCortab$pvalue, pi0=1)
    })
    allCortab <- cbind(allCortab, qvals$qvalues)
    colnames(allCortab)[ncol(allCortab)] <- "qvalue"
  }


  return(allCortab)
}


