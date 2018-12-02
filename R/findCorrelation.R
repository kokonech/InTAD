
findGeneCorrelation <- function(x, signalVals, countVals, corMethod) {
    genes <- x$geneid
    pId <- as.character(unique(x$peakid))
    #print(pId)

    corDat <- list()
    for (i in seq_len(length(genes))) {
        if (! (as.character(genes[i]) %in% row.names(countVals))) {
            #print(paste("NOT FOUND!", genes[i]))
            next
        }
        # log2 is taken into account in the generation of object
        countSel <- as.numeric( countVals[as.character(genes[i]), ])
        #print(as.character(genes[i]) )
        if ( sum(countSel) == 0) {
            next
        }
        sigSel = as.numeric(signalVals[pId,])
        cors <- suppressWarnings(cor.test( sigSel, countSel, method=corMethod))
        euDis = as.numeric(dist(rbind ( sigSel, countSel)))
        corDat[[i]] <- data.frame(gene=genes[i], eudis= euDis,
                                corr=cors$estimate,
                                pvalue=cors$p.value,
                                stringsAsFactors =FALSE)
    }
    corDF <- do.call(rbind, lapply(corDat, function(x) x))
    return(corDF)
}


#' Function to perfrom correlation analysis
#'
#' This function combines genes and signals in inside of TADs
#' @param object InTADSig object with signals and genes combined in TADS
#' @param method Correlation method: "pearson" (default), "kendall", "spearman"
#' @param adj.pval Perform p-value adjsutment and include q-values in result
#' @param plot.proportions Plot proportions of signals and genes in correlation
#' @return A table with correlation values for signal-gene pairs including
#' correlation p-value, euclidian distance and rank.
#' @importFrom stats na.omit cor.test dist
#' @import qvalue
#' @export
#' @examples
#' ## perform analysis on test data
#' inTadSig <- newSigInTAD(enhSel, enhSelGR, rpkmCountsSel, txsSel)
#' inTadSig <- filterGeneExpr(inTadSig, geneType = "protein_coding")
#' inTadSig <- combineInTAD(inTadSig, tadGR)
#' corData <- findCorrelation(inTadSig, method="pearson")
#'

findCorrelation <- function(object, method = "pearson", adj.pval = FALSE,
                            plot.proportions = FALSE ) {

    if (!is(object, "InTADSig"))
        stop("Object must be an InTADSig!")

    if (length(object@signalConnections)  == 0)
        stop("No signals and genes are combined in TAD!
            Use combineInTAD() function.")

    sigCons <- object@signalConnections
    txs <- rowRanges(object@sigMAE[["exprs"]])

    # this null filtering was already performed in previous step
    allIdGnX <- sigCons[!sapply(sigCons, is.null)]
    allIdGnY <- allIdGnX[sapply(allIdGnX, function(x) dim(x)[1]) > 0]

    if (object@ncore > 1 && requireNamespace("parallel")) {
        message("Running in parallel...")
        allRes <- parallel::mclapply(allIdGnY, findGeneCorrelation,
                assay(object@sigMAE[["signals"]]),
                assay(object@sigMAE[["exprs"]]),
                method)
    } else {
        allRes <- lapply(allIdGnY, findGeneCorrelation,
                    assay(object@sigMAE[["signals"]]),
                    assay(object@sigMAE[["exprs"]]),
                    method)
    }

    allRes <- lapply(allRes, function(x) na.omit(x))

    allX <- list()
    for(i in seq_len(length(allRes))) {
        gnname <- values(txs)$gene_name[match(allRes[[i]]$gene,
                                            values(txs)$gene_id)]
        rnk <- allRes[[i]]$gene[order(allRes[[i]]$corr, decreasing=TRUE)]
        corRank <- match(allRes[[i]]$gene, rnk)
        #cluster <- cltab$cluster[match(names(allRes)[i], cltab$id)]
        peakName <- names(allRes)[i]
        tad <- unique(allIdGnX[[peakName]]$tad)
        #tad <- inTADpairs$tad[match(names(allRes)[i], inTADpairs$peakid)]
        allX[[i]] <- data.frame(peakid=rep(peakName, length(gnname)),
                            tad=rep(tad,length(gnname)),
                            gene=allRes[[i]]$gene, name=gnname,
                            cor=allRes[[i]]$corr,
                            pvalue=allRes[[i]]$pvalue,
                            eucDist = allRes[[i]]$eudis,
                            corRank=corRank,
                            stringsAsFactors =FALSE)
    }

    allCortab <- do.call(rbind, lapply(allX, function(x) x))

    if (adj.pval) {
        qvals <- tryCatch({
            qvalue(corData$pvalue)
        }, error = function(e) {
            qvalue(corData$pvalue, pi0=1)
        })
        allCortab <- cbind(allCortab, qvals$qvalues)
        colnames(allCortab)[9] <- "qvalue"
    }

    if (plot.proportions) {
        corDataF <- allCortab[allCortab$pvalue < 0.05,]
        if (nrow(corDataF) == 0) {
            warning("No conections detected")
            return(allCortab)
        }
        # plot 1
        assignedGeneProportion <- summary(as.factor(corDataF$gene),
                                            maxsum = nrow(corDataF))
        agpSorted <- sort(assignedGeneProportion, decreasing = TRUE)

        h <- hist( agpSorted, plot =FALSE )
        h$density = h$counts/sum(h$counts)*100
        plot(h, freq=FALSE, xlab="Number of signals",
            ylab = "Assigned genes (%)",
            col="blue",cex.main = 0.7,cex.axis=0.8,cex.lab=0.7,
            main = "Genes regulated by signals (p-val 0.05)")

        # plot 2
        assignedEnhProportion <- summary(as.factor(corDataF$peakid),
                                            maxsum = nrow(corDataF))
        aepSorted <- sort(assignedEnhProportion, decreasing = TRUE)

        h <- hist( aepSorted, plot =FALSE )
        h$density = h$counts/sum(h$counts)*100
        plot(h, freq=FALSE, xlab="Number of genes",
            ylab = "Assigned signals (%)",
            col="blue", cex.main = 0.7,cex.axis=0.8,cex.lab=0.7,
            main = "Signals regulating genes (p-val 0.05)")

    }

    return(allCortab)
}
