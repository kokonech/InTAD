### Methods for the InTADSig class

#' Create InTADSig object
#'
#' The fuction generates an object that contains the signals and gene expression
#' data.frames along with their genomic coordinates for further processing.
#'
#' @param signalData data frame containing signals
#' @param signalRegions genomic regions of the signals
#' @param countsData data matrix containing count expression values
#' @param geneRegions gene coordiantes
#' @param sampleInfo data frame containing additional sample info
#' @param performLog  Perform log2 convertion of expression values.
#' Default: TRUE.
#' @param logExprsOffset Offset x for log2 gene exrpression
#' i.e. log2(value + x). Default: 1
#' @param ncores Number of cores to use for parallel computing
#' @return Novel InTADSig object
#
#' @details
#' InTADSig object stores matrices of signals and gene expression values along
#' with coordinates. The order of samples and names of columns should match in both
#' datasets. For gene coordinates GRanges "gene_id" and "gene_name" are required
#' in metadata. These are typical markers of genes in GTF anntotation format.
#'
#' @import methods
#' @import SummarizedExperiment
#' @import MultiAssayExperiment
#' @importFrom S4Vectors DataFrame SimpleList
#' @export
#' @examples
#' ## create sigInTAD object
#' inTadSig <- newSigInTAD(enhSel, enhSelGR, rpkmCountsSel, txsSel)
#'
#'


newSigInTAD <- function(signalData = NULL,
                        signalRegions = NULL,
                        countsData = NULL,
                        geneRegions = NULL,
                        sampleInfo = NULL,
                        performLog = TRUE,
                        logExprsOffset = 1,
                        ncores = 1)
{
    ## Check that we have some expression data
    if ( is.null(signalData) || is.null(countsData))
        stop("Require signal and gene expression input data.frames.")

    if (!is(signalRegions, "GRanges"))
        stop("signalRegions must be GRanges object!")

    if (!is(geneRegions, "GRanges"))
        stop("geneRegions must be GRanges object!")

    if (is.null(sampleInfo)) {
        sampleInfo <- DataFrame(Type=paste0("Sample",seq_len(ncol(signalData))),
                                row.names=colnames(signalData))
    }

    # GRanges names should be ignored, only row names of tables are valid
    signalSE <- SummarizedExperiment(
        assays=SimpleList(counts=as.matrix(signalData, rownames.force=TRUE)),
        rowRanges=signalRegions, colData=sampleInfo)

    #names(rowRanges(signalSE)) <- rownames(assay(signalSE))
    exprSE <- SummarizedExperiment(
        assays=SimpleList(counts=as.matrix(countsData)),
        rowRanges=geneRegions, colData=sampleInfo)


    inTadMAE <- MultiAssayExperiment(list("signals"=signalSE,"exprs"=exprSE),
                                        sampleInfo)

    # Generate new object
    inTadSig <- new( "InTADSig", sigMAE = inTadMAE )



    # Check validity of object
    validObject(inTadSig)

    # Assign names for signal coords; is this required for MAE?
    # names(inTadSig@sigCoords) <- rownames(inTadSig@signals)

    # Log2 gene expression
    if (performLog) {
        assay(inTadSig@sigMAE[["exprs"]]) <-
            log2(assay(inTadSig@sigMAE[["exprs"]]) + logExprsOffset)
    }

    message(paste("Created signals and genes object for",
                    ncol(signalData) ,"samples"))
    inTadSig@ncore <- ncores
    if (ncores > 1) {
        pkgAvail <- requireNamespace("parallel")
        if (pkgAvail) {
            message(paste("Activate multicore computation. Num cores:",ncores))
            options(mc.cores=ncores)
        } else {
            stop("Package parallel is required for multi-core anlaysis!")
        }
    }

    inTadSig
}


#' Load InTADSig object from text files
#'
#' The fuction loads the data tables to create an object that contains the
#' signals and gene expression data.frames along with their genomic coordinates
#' for further processing.
#'
#' @param signalsFile Tab-seprated data table containg signals and their
#' coordinates as row.names
#' @param countsFile Tab-seprated counts table
#' @param gtfFile GTF file containing all gene coordinates
#' @param annFile Tab-delimited phenotype annotation of samples
#' @param performLog  Perform log2 convertion of expression values.
#' Default: TRUE.
#' @param logExprsOffset Offset x for log2 gene exrpression
#' i.e. log2(value + x). Default: 1
#' @param ncores Number of cores to use for parallel computing
#' @return Novel InTADSig object
#'
#' @importFrom rtracklayer import.gff
#' @importFrom utils read.delim
#' @import GenomicRanges
#'
#' @details
#' The function loads data from input files and creates object that stores
#' matrices of signals and gene expression values along with coordiantes.
#' The samples order and names of columns should match in both tables.
#' It is expected that gene ids are applied in the validation of counts table.
#'
#' @export
#' @examples
#' # create sigInTAD object
#' inTadSig <- newSigInTAD(enhSel, enhSelGR, rpkmCountsSel, txsSel)
#'
#'


loadSigInTAD <- function(signalsFile,
                        countsFile,
                        gtfFile,
                        annFile = "",
                        performLog = TRUE,
                        logExprsOffset = 1,
                        ncores = 1 )
{

    message("Loading signals...")
    sigtab <- read.delim(signalsFile)
    if (sum(colnames(sigtab)[seq_len(3)] ==   c("chr","start","end")) != 3)
        stop("First 3 columns of signal table should be genomic coordiantes
            with names \"chr\",\"start\",\"end\"!")

    sigCoords = sigtab[,c(1,2,3)]
    sigGR <- GRanges(sigCoords)
    sigData <- sigtab[, 4:ncol(sigtab)]

    ids <- paste(sigtab$chr,":", sigtab$start, "-", sigtab$end, sep="")
    rownames(sigData) <- ids

    message("Loading gene expression...")
    geneCounts <- read.delim(countsFile)

    message("Loading gene coordinates...")
    annData <- import.gff(gtfFile)
    if ("type" %in%  names(mcols(annData))) {
        message("Selecting only gene coordinates...")
        annData <- annData[annData$type == "gene",]
    }

    if (!"gene_id" %in% colnames(values(annData)))
        stop("Anntotation parameter \"gene_id\" is missing! Check GTF file.")

    # check order
    if (!"gene_id" %in% colnames(values(annData)))
        stop("Anntotation parameter \"gene_id\" is missing! Check GTF file.")

    numSameIds = sum(rownames(geneCounts) == annData$gene_id )
    if ( numSameIds != length(annData)) {
        message("Fixing row order of gene counts to correspond annotation...")
        geneCounts <- geneCounts[ annData$gene_id, ]

    }

    sInfo <- NULL
    if (nchar(annFile) > 0 ) {
        message("Loading annotation...")
        sInfo <- read.delim(annFile, header = TRUE,stringsAsFactors = FALSE)
    }

    # create object
    inTadSig <- newSigInTAD(sigData, sigGR, geneCounts, annData,
                            sampleInfo = sInfo,
                            performLog = performLog,
                            logExprsOffset = logExprsOffset,
                            ncores = ncores )

    inTadSig

}
