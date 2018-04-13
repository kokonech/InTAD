# Classes

#' The InTADSig Class
#'
#' The InTADSig object stores signals and gene expression data for the samples.
#'
#' It uses MultiAssayExperiment object to store information. Key slots to access
#' are listed below.
#'
#'
#'@section Slots:
#'    \describe{
#'        \item{\code{sigMAE}:}{\code{"MultiAssayExperiment"},
#'        MultiAssayExperiment object containg signals and gene counts }
#'        \item{\code{signalConnections}:}{\code{"list"},
#'        The list of signals representing gene data frames in the same TAD }
#'        \item{\code{ncore}:}{\code{"numeric"},
#'        Number of cores to use for parallel computing }#'
#'}
#' @name InTADSig
#' @rdname InTADSig
#' @aliases InTADSig-class
#' @exportClass InTADSig
#'
InTADSig <- setClass("InTADSig",
    slots = c(
        sigMAE = "MultiAssayExperiment",
        signalConnections = "list",
        ncore = "numeric"
    )
)


################################################################################
### Validity check for InTADSig class object

setValidity("InTADSig", function(object) {

    msg <- NULL
    valid <- TRUE

    numSignalSamples <- ncol(signals(object))
    numExprSamples <- ncol(exprs(object))

    if ( (numSignalSamples != 0) && (numSignalSamples != numExprSamples) )  {
        valid <- FALSE
        msg <- c(msg, "Number of signal samples doesn't
                match number of RNA-seq samples")
    }

    if (  (length(rowRanges(object@sigMAE[["signals"]])) != 0) &&
        !identical(colnames(assay(object@sigMAE[["signals"]])),
                    colnames(assay(object@sigMAE[["signals"]]))) ) {
        valid <- FALSE
        msg <- c(msg, "Signal sample names must match RNA-seq sample names")
    }

    geneGR <- rowRanges(object@sigMAE[["exprs"]])

    if (!"gene_id" %in% colnames(values(geneGR)))
        stop("Annotation parameter \"gene_id\" is missing!")

    if (!"gene_name" %in% colnames(values(geneGR)))
        stop("Annotation parameter \"gene_name\" is missing!")

    numSameIds = sum(rownames(exprs(object)) == geneGR$gene_id )
    if ( numSameIds != length(geneGR))
        stop("The \"gene_id\" values in GRanges are not expr matrix row names!")

    if (valid) TRUE else msg
})


setMethod('show', signature='InTADSig', definition=function(object) {
    cat("S4 InTADSig object\n")
    signals <- signals(object)
    cat(paste("Num samples:", ncol(signals),"\n"))
    cat(paste("Num signals:", nrow(signals), "\n"))
    cat(paste("Num genes: ", nrow(exprs(object)), "\n"))
    #if (is.na(tads(object)))
    #    cat(paste("TAD regions are not set"))
    #else:
    #    cat(paste("TAD regions: ", length(tads(object), '\n'))
    }
)


#### Methods for extracting slots ###

#' @rdname signals
#' @export
setGeneric('signals', function(object) standardGeneric('signals'))

#' Signal values table
#'
#' This funcion returns the signal values table
#' @param object InTADSig object with signals and genes
#' @return Signals table
#' @rdname signals
#' @examples
#' inTadSig <- newSigInTAD(enhSel, enhSelGR, rpkmCountsSel, txsSel)
#' head(signals(inTadSig))
#' @export
setMethod('signals', 'InTADSig',
            function(object) assay(object@sigMAE[["signals"]]))


#' @rdname sigCoords
#' @export
setGeneric('sigCoords',
            function(object) standardGeneric('sigCoords'))

#' Signal coords GRanges
#'
#' This funcion returns the signal GRanges
#' @param object InTADSig object with signals and genes
#' @return Signal GRanges
#' @rdname sigCoords
#' @examples
#' inTadSig <- newSigInTAD(enhSel, enhSelGR, rpkmCountsSel, txsSel)
#' head(sigCoords(inTadSig))
#' @export
setMethod('sigCoords', 'InTADSig',
            function(object) rowRanges(object@sigMAE[["signals"]]))


#' Gene expression counts table
#'
#' This funcion returns gene expression counts table
#' @param object InTADSig object with signals and genes
#' @return Gene expression table
#' @rdname exprs
#' @examples
#' inTadSig <- newSigInTAD(enhSel, enhSelGR, rpkmCountsSel, txsSel)
#' head(exprs(inTadSig))
#' @export
setMethod('exprs', signature='InTADSig',
          definition=function(object) assay(object@sigMAE[["exprs"]]))

# get pairs in TADs
setGeneric('sigLinks', function(object) standardGeneric('sigLinks'))
setMethod('sigLinks', 'InTADSig', function(object) object@signalConnections)

#' @rdname geneCoords
#' @export
setGeneric('geneCoords', function(object) standardGeneric('geneCoords'))

#' Gene coords GRanges
#'
#' This funcion returns the gene GRanges
#' @param object InTADSig object with signals and genes
#' @return Gene GRanges
#' @rdname geneCoords
#' @examples
#' inTadSig <- newSigInTAD(enhSel, enhSelGR, rpkmCountsSel, txsSel)
#' head(geneCoords(inTadSig))
#' @export
setMethod('geneCoords', 'InTADSig',
            function(object) rowRanges(object@sigMAE[["exprs"]]))












