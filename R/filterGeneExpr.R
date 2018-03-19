#' Function to filter gene expression
#'
#' This function performs filtering of gene expression counts based on various
#'  parameters
#' @param obj InTADSig object
#' @param cutVal Exclude genes that have max expression less or equal to
#' this value in all samples. Default: 0
#' @param geneType Type of gene to select for filtering i.e. "protein_coding".
#' Default:NA
#' @param checkExprDistr Adjust cutVal based on gene expression distribution
#' @param plotExprDistr Perform visualziation of the distribution
#' @return InTADSig object with filtered counts table
#' @details
#' The function allows to stabilize the functional activity of the genes. By
#' default all not expressed genes are filtered. It is also possible to set type
#' of gene to take into account i.e. "protein_coding" only. This option requires
#' additional metadata column "transcript_type".
#' Also, special filtering option based on mclust library allows to analyze
#' distribution of counts and adjust the cut value to exclude low expressed
#' genes.
#' @importFrom stats qnorm dnorm density
#' @export
#' @examples
#' ## perform analysis on test data
#' inTadSig <- newSigInTAD(enhSel, enhSelGR, rpkmCountsSel, txsSel)
#' ## default filtering
#' inTadSig <- filterGeneExpr(inTadSig)
#' ## filter based on gene type
#' inTadSig <- filterGeneExpr(inTadSig, geneType = "protein_coding")
#'
#'
filterGeneExpr <- function(obj, cutVal = 0, geneType = NA,
                        checkExprDistr = FALSE, plotExprDistr = FALSE ) {

    if (!is(obj, "InTADSig"))
        stop("Object must be an InTADSig!")

    if (length(obj@signalConnections)  != 0)
        warning("Signals and genes are already combined in TAD,
            rerun combineInTAD() function to update the result.")

    # TODO: create an object to keep initial gene expression

    RNA <- assay(obj@sigMAE[["exprs"]])
    txs <- rowRanges(obj@sigMAE[["exprs"]])

    message(paste("Initial result:", nrow(RNA), "genes"))

    if (checkExprDistr) {
        vals <- as.vector(RNA)
        # P value <--> log2 ratio
        # pval2ratio <- function(pval,  background)
        # {qnorm(pval,mean=background[1], sd=background[2], lower.tail=FALSE)}

        ### define positive/negative cut-off (by fitting normal mixture model)
        pcut <- 0.01

        d1 <- density(vals,na.rm=TRUE)
        bg1 <- get.enr.bg.normfit(vals)
        if (sum(is.na(bg1)) == 0) {
            #ratiocut <- pval2ratio(pcut, background=bg1)
            ratiocut <- qnorm(pcut,  mean=bg1[1], sd=bg1[2], lower.tail=FALSE)

            if (plotExprDistr) {
                bgProp1 = attr(bg1, "proportions")[1]
                bgProp2 = attr(bg1, "proportions")[2]
                plot(d1, col='gray', lwd=3, main="Normal Mixture Model",
                    ylim=c(0,2.5),xlab="Expression", ylab="Density")
                lines(x=d1$x,
                    y=dnorm(d1$x, mean=bg1[1], sd=bg1[2])*bgProp1,
                    lwd=3, col="blue")
                lines(x=d1$x,
                    y=dnorm(d1$x, mean=attr(bg1, "enr.mean"),
                    sd=attr(bg1, "enr.sd")) * bgProp2,
                    lwd=3, lty=2, col="green")
                abline(v=ratiocut, lwd=3, col="red")
                legend(x="topright", bty="n", lty=c(1,1,2,1,NA),
                    lwd=c(3,3,3,3,NA),
                    col=c("gray","blue","green","red","black"),
                    legend=c(sprintf("alldata (n=%d)", d1$n),
                        sprintf("background (%.1f%%)", 100*bgProp1),
                        sprintf("foreground (%.1f%%)", 100*bgProp2),
                        sprintf("P = %.3g (n_sig=%d)",pcut,
                            length(which(vals>=ratiocut)))))
            }

            cutVal <- ratiocut
        }


    }

    message(paste("Gene expression cut value:", cutVal))

    maxR <- apply(RNA, 1, function(x) max(x))
    RNA <- RNA[which(maxR>cutVal),]
    txs <- txs[ txs$gene_id %in% rownames(RNA)  ]

    if (!is.na(geneType)) {
        cols <- colnames(as.data.frame(txsSel))
        if ("gene_type" %in% cols) {
            txsSel <- txs[which(txs$gene_type==geneType),]
            if (length(txsSel) == 0) {
                stop(sprintf("The gene type \"%s\" is not found!", geneType))
            }
            geneIds <- unique(txsSel$gene_id)
            txs <- txsSel
            RNA <- RNA[  rownames(RNA) %in% geneIds , ]
        } else {
            warning("The metadata column \"gene_type\" is not present in
                gene coordinates. Skipping this step...")
        }

    }

    message(paste("Filtered result:", nrow(RNA), "genes"))
    #obj@exprs <- RNA
    #obj@geneCoords <- txs

    exprsInitial <- experiments(obj@sigMAE)[["exprs"]]
    experiments(obj@sigMAE)[["exprs"]] <- exprsInitial[rownames(RNA),]

    return(obj)

}

#' Function to estimate gene expression
#'
#' This function uses mclust package to analyze gene expression distribution
#' @param x Full gene expression vector
#' @return Distribution properties: mean and std
#' @import mclust
#' @details
#' The function adjust filtering cut value based on mclust library to exclude
#' low expressed genes. It is a part of filtering procedure.
get.enr.bg.normfit <- function(x) {

    result <- c(NA, NA)

    # NOTE: no idea how not to depend on mclust
    # requireNamespace() did not work properly

    if (requireNamespace("mclust")) {
        names(result) <- c("mean","sd")

        # fit two normal distributions to data
        # model <- Mclust(na.omit(x), G=2, modelNames="V")
        model <- mclust::Mclust(na.omit(x), G=2)

        # bg corresponds to first (lower mean) normal
        result[1] <- model$parameters$mean[1]
        result[2] <- sqrt(model$parameters$variance$sigmasq[1])

        # store additional information in attributes
        attr(result, "proportions") <- model$parameters$pro
        attr(result, "enr.mean") <- as.numeric(model$parameters$mean[2])
        attr(result, "enr.sd") <- sqrt(model$parameters$variance$sigmasq[2])

        # return results
    } else {
        message("Mclust library not detected, skipping cut value adjustment...")
    }
    result
}
