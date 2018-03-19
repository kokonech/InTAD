#' Function to plot correlation
#'
#' This function creates a plot of selected pair signal-gene
#' @param obj InTADSig object with signals and genes combined in TADS
#' @param sId Signal id based on genomic cooridantes i.e. "chr:start-end"
#' @param geneName Gene name to select. Based on "gene_name" attribute.
#' @param xLabel The label to mark signal X-axis. Default: "Gene expression"
#' @param yLabel The label to mark signal Y-axis. Default: "Signal enrichment"
#' @param colByPhenotype The pheno data column i.e. tumour type
#' that can be use for colour
#' @param corMethod Correlation method. Default: Pearson
#' @importFrom stats cor
#' @import ggpubr
#' @import ggplot2
#' @return A \code{ggplot} object for visualization or customization.
#' @examples
#' inTadSig <- newSigInTAD(enhSel, enhSelGR, rpkmCountsSel, txsSel)
#' inTadSig <- combineInTAD(inTadSig, tadGR)
#' plotCorrelation(inTadSig, "chr15:26372163-26398073", "GABRA5")
#'
#' @export
plotCorrelation <- function( obj, sId, geneName,
                            xLabel = "Gene expression",
                            yLabel = "Signal enrichment",
                            colByPhenotype = "",
                            corMethod = "pearson") {

    if (!is(obj, "InTADSig"))
        stop("Object must be an InTADSig!")

    if (length(obj@signalConnections)  == 0)
        stop("No signals and genes are combined in TAD!
            Use combineInTAD() function")

    if (nchar(colByPhenotype) > 0) {
        if (!( colByPhenotype %in% colnames(colData(obj@sigMAE))) )
            stop(paste0("Phenotype ", colByPhenotype," is not found!"))
    }

    allIdGnX <- obj@signalConnections

    if (!sId %in% names(allIdGnX)  ) {
        stop("Signal is not found!")
    }

    selSig <- allIdGnX[[sId]]

    if (!geneName %in% selSig$names) {
        stop(sprintf("The gene %s is not in the
                    same TAD as signal %s", geneName, sId))
    }
    geneId <- selSig[selSig$names == geneName,]$geneid

    toPlot<-data.frame(
        sig = signals(obj)[row.names(signals(obj))==sId, ],
        exp = exprs(obj)[row.names(exprs(obj))==geneId, ]
    )

    selColor <- "black"
    if (nchar(colByPhenotype) > 0) {
        toPlot <- cbind( toPlot, colData(obj@sigMAE)[, colByPhenotype] )
        colnames(toPlot)[3] <- colByPhenotype
        selColor <- colByPhenotype
    }

    title <- sprintf("%s - %s (cor = %.3f)",sId, geneName,
                    cor(toPlot$exp, toPlot$sig,method = corMethod))

    res <- ggscatter(toPlot, x="exp", y="sig",
                    color=selColor,
                    #palette = group.colors,
                    size = 4)+
        xlab(xLabel)+
        ylab(yLabel)+
        ggtitle(title)

    res

}



#' Function to plot correlation across genome
#'
#' This function creates a plot of correlation strength
#' in target genomic region from the result table.
#' The X-coordinates represent signals, Y-coords represent genes, while
#' each dot represents -log10(P-value) from correlation test.
#' Additionallly all TAD boundaries can be visualized.
#' @param obj InTADSig object with signals and genes combined in TADS
#' @param corRes Correlation result table created by function findCorrelation()
#' @param targetRegion Target genomic region visualise.
#' @param showCorVals Use this option to visualize postive correlation values
#' instead of correlation strength
#' @param tads TAD regions to visualize. By default only TADs persent in
#' correlation result table are applied (NULL value).
#'
#' @import ggplot2
#' @return A \code{ggplot} object for visualization or customization.
#' @examples
#' inTadSig <- newSigInTAD(enhSel, enhSelGR, rpkmCountsSel, txsSel)
#' inTadSig <- combineInTAD(inTadSig, tadGR)
#' corData <- findCorrelation(inTadSig, method="pearson")
#' plotCorAcrossRef(inTadSig,corData,GRanges("chr15:25000000-28000000"))
#'
#' @export
plotCorAcrossRef <- function( obj, corRes, targetRegion,
                              showCorVals = FALSE,
                              tads = NULL) {

    if (!is(obj, "InTADSig"))
        stop("Object must be an InTADSig!")

    if (!is(targetRegion, "GRanges"))
        stop("Target region must be GRanges")

    if (sum( colnames(corRes)[seq_len(3)]  == c("peakid","tad","gene")) != 3)
        stop("Incorrect correlation table! Expected first 3 table
            column names: peakid, tad, gene")

    peaks <- GRanges(corRes$peakid)
    enhOverlap <- findOverlaps(sigCoords(obj), targetRegion)
    if (length(enhOverlap) == 0) {
        stop("No signals inside selected region detected.")
    }
    enhToCheck <- sigCoords(obj)[queryHits(enhOverlap)]

    # peaks within regions
    corSel <- corRes[ corRes$peakid %in% names(enhToCheck), ]
    genes <- geneCoords(obj)[unique(corRes$gene)]
    genesOverlap <- findOverlaps(genes, targetRegion,type = "within")
    genesToCheck <- genes[queryHits(genesOverlap)]
    # both peaks and genes within regions
    corSel2 <- corSel[ corSel$gene %in% genesToCheck$gene_id, ]

    sigSel <- GRanges(corSel2$peakid)
    geneSel <- genesToCheck[ corSel2$gene ]

    if (showCorVals) {
        res <- cbind(start(sigSel), end(sigSel),
                start(geneSel), end(geneSel),
                corSel2$cor)
    } else {
        res <- cbind(start(sigSel), end(sigSel),
                   start(geneSel), end(geneSel),
                   -log10(corSel2$pvalue))

    }

    # start position version
    dt <- as.data.frame(res[,c(1,3,5)])
    colnames(dt) <- c("enh","gene","cor")
    if (showCorVals) {
      dt <- dt[dt$cor > 0, ]
    }

    if (nrow(dt) == 0)
        stop("No correlation between signal and genes
            found inside selected region")

    if (showCorVals) {
        legendLab = "Cor"
    } else {
        legendLab = "-log10(P)"
    }

    sp <- ggplot( dt, aes_string(x="enh",y="gene",color="cor")) +
        geom_point() + labs(color=legendLab) +
        labs(title = paste("Region",as.character(targetRegion)),
            x = "Signal coords (bp)", y="Gene coords (bp)") +
        theme(plot.title = element_text(hjust = 0.5)) +
        expand_limits( x=c(start(targetRegion), end(targetRegion)),
                    y=c(start(targetRegion), end(targetRegion))) +
        scale_color_gradient(low = "white", high = "red")

    # add TAD borders

    if (!is.null(tads)) {
        if (is(tads, "GRanges")) {
            selTadGR <- tads
        }
    } else {
        selTadGR <-  GRanges(unique(corRes$tad))
    }

    tadOverlap <- selTadGR[queryHits(findOverlaps(selTadGR, targetRegion))]
    yStartLim = min( c(min(start(tadOverlap)),start(targetRegion) ))
    yEndLim = max( c(max(end(tadOverlap)),end(targetRegion) ))

    if (length(tadOverlap) > 0) {
        sp <- sp + geom_rect(data = as.data.frame(tadOverlap),
                            inherit.aes = FALSE,
                            aes(xmin = start, xmax=end,ymin=start, ymax=end),
                            fill=NA, color="black", linetype=3) +
            coord_cartesian(ylim=c(start(targetRegion), end(targetRegion)),
                        xlim=c(start(targetRegion), end(targetRegion)))
    }

    diagDf <- as.data.frame(rbind( rep(start(targetRegion),2),
                                 rep(end(targetRegion),2) ))
    colnames(diagDf) <- c("start","end")
    sp <- sp + geom_line(data = diagDf,inherit.aes = FALSE,
                         aes(x=start,y=end), col="gray",linetype=8)

    sp <- sp + scale_y_reverse(lim=c(yEndLim, yStartLim))
    sp

}
