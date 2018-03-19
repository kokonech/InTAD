#' Preparation for correlation analysis for a signal
#'
#' This function collects all genes for signal genomic region inside of
#' Topologically Associated Domains (TADs)
#' @param id Id of signal from the list
#' @param sigList List of signal GRs and their names
#' @param tadGR TAD genomic regions
#' @param tss Gene transcription start sites
#' @param pickMaxOvlp Use TAD with max overlap
#' @param nearestTad The table listing TADs nearest to each TSS
#'#'
#' @details
#' The signal is checked if it is lying inside of TAD.
#' Then all genes in this TAD are collected.
#'
#' @return Data.frame containing genes connected to signal
#'
fnSE <- function(id, sigList, tadGR, tss,  pickMaxOvlp, nearestTad) {

    sGR <- sigList[[id]]
    sName <- names(sigList)[id]
    #message(paste("Checking peak...", sName, sGR))
    ov1 <- findOverlaps(query=sGR, subject=tadGR)
    if (length(subjectHits(ov1)) == 0) {
        return(NULL)
    }

    if ( (length(subjectHits(ov1))>1) & pickMaxOvlp) {
        #message(paste("Overlapping:",sGR) )
        temp <- tadGR[subjectHits(ov1),]
        intsW <- NULL
        for(k in seq_len(start(temp))) {
            ints <- intersect(sGR,temp[k,])
            intsW[k] <- width(ints)
        }
        ovTad <- temp[which.max(intsW),]
    } else {
        ovTad <- tadGR[subjectHits(ov1),]
    }

    if (is.null(nearestTad)) {
        ovTss <- findOverlaps(query=tss, subject=ovTad, type="within")
        dattab <- data.frame(peakid=rep(sName,length(queryHits(ovTss))),
            geneid=values(tss)$geneid[queryHits(ovTss)],
            names=names(tss)[queryHits(ovTss)],
            tad=names(ovTad[subjectHits(ovTss)]),
            stringsAsFactors = FALSE)
    } else {

        closestGenes <- nearestTad[nearestTad$TAD %in% names(ovTad),,drop=FALSE]
        geneSel <- tss [ closestGenes$Gene ]

        dattab <- data.frame(peakid=rep(sName,nrow(closestGenes)),
            geneid=geneSel$geneid,names=names(geneSel),tad=closestGenes$TAD,
            stringsAsFactors = FALSE)

    }

    dattab2 <- unique(dattab)
    return(dattab2)
}

#' Preparation for correlation analysis
#'
#' This function combines signals and genes in inside of
#' Topologically Associated Domains (TADs)
#' @param object InTADSig object
#' @param tadGR TAD genomic regions
#' @param selMaxTadOvlp If a signal overlaps 2 or more TADs by default only
#' single TAD with max overlap is selected.All overlaps can be included by
#' deactivating this option.
#' @param closestGene By default closest to TAD genes are selected based
#' on TSS location. Deactivate this option to use genes only lying within TAD.
#'
#' @details
#' Each signal is checked if it is lying inside of TAD. Signals out of TADs
#' are ignored. The genomic regions reprenting gene coordiantes are converted
#' to TSS. By default, the closest genes are assigned belonging to TAD.
#' If this option deactivated, only those lying with TAD are collected.
#' Result is a list of signals connected to tables with gene details.
#'
#' @importMethodsFrom GenomicRanges intersect
#' @importFrom BiocGenerics start end
#' @importFrom S4Vectors queryHits subjectHits getListElement
#' @importFrom IRanges IRanges
#' @import graphics
#'
#' @return Updated InTADSig object containing genes connected to eash signal
#'
#' @export
#' @examples
#' # create sigInTAD object
#' inTadSig <- newSigInTAD(enhSel, enhSelGR, rpkmCountsSel, txsSel)
#' # combine signals and genes in TAD
#' inTadSig <- combineInTAD(inTadSig, tadGR)
#'
#'
combineInTAD <- function( object, tadGR, selMaxTadOvlp = TRUE,
                            closestGene = TRUE) {

    if (!is(object, "InTADSig"))
        stop("Object must be an InTADSig!")

    # make sure TAD regions have names
    if (is.null(names(tadGR))) {
        names(tadGR) <- as.character(tadGR)
    }
    # coords requrired
    geneGR <- rowRanges(object@sigMAE[["exprs"]])
    tss <- GRanges( seqnames=as.character(seqnames(geneGR)),
                IRanges(start=ifelse(as.character(strand(geneGR))=="+",
                start(geneGR),end(geneGR)), width=1,
                names=values(geneGR)$gene_name),
                strand=strand(geneGR), geneid=values(geneGR)$gene_id )

    nearestTad <- NULL
    if (closestGene) {
        nearestTss <- nearest(tss,tadGR)
        validTss <- !(is.na(nearestTss))
        nearestTad <- as.data.frame( cbind( names(tss)[validTss],
                                    names(tadGR[nearestTss[validTss]])),
                                    stringsAsFactors=FALSE )
        colnames(nearestTad) <- c("Gene", "TAD")
    }

    # fix issue with list
    sigGR <- rowRanges(object@sigMAE[["signals"]])
    sigList <- as( sigGR, "GRangesList")

    if (object@ncore > 1 && requireNamespace("parallel") ) {
        message("Running in parallel...")
        object@signalConnections <-
            parallel::mclapply( seq_along(sigList), fnSE,
                                sigList, tadGR, tss, selMaxTadOvlp,nearestTad )
    } else {
        object@signalConnections <-
            lapply( seq_along(sigList),fnSE,
                    sigList, tadGR, tss, selMaxTadOvlp, nearestTad)
    }
    names(object@signalConnections) <- names(object@sigMAE[["signals"]])

    # clear enhancers with no connection
    object@signalConnections <-
        object@signalConnections[!sapply(object@signalConnections, is.null)]
    message(paste("Combined",sum(sapply(object@signalConnections, nrow)),
                    "signal-gene pairs in TADs" ))
    return(object)

}

