## Initial testing

context("test on inputs")

test_that("signals and exprs input data are present", {
    expect_that(newSigInTAD(),
                throws_error("Require signal and gene expression input data.frames."))
})

test_that("genomic regions of signals are present", {
    data("mbSelEnhSignals")
    data("mbSamplesRPKM")
    
    expect_that(newSigInTAD(enhSel, list(), rpkmCountsSel, list()),
                throws_error("signalRegions must be GRanges object!"))
})


test_that("genomic regions of gene expresssion are present", {
    data("mbSelEnhSignals")
    data("enhSelCoords")
    data("mbSamplesRPKM")
    
    expect_that(newSigInTAD(enhSel, enhSelGR, rpkmCountsSel, list()),
                throws_error("geneRegions must be GRanges object!"))
})

# Control 
# Test: less region coords than regions
# Test: less gene coords than genes
# Test: less samples 

test_that("genomic regions have gene_id", {
    data("mbSelEnhSignals")
    data("enhSelCoords")
    data("mbSamplesRPKM")
    data("txsSel")

    txsSel$gene_id <- NULL
    expect_that(newSigInTAD(enhSel, enhSelGR, rpkmCountsSel,txsSel ),
                throws_error("Annotation parameter \"gene_id\" is missing!"))
})


test_that("genomic regions have gene_name", {
    data("mbSelEnhSignals")
    data("enhSelCoords")
    data("mbSamplesRPKM")
    data("txsSel")

    txsSel$gene_name <- NULL
    expect_that(newSigInTAD(enhSel, enhSelGR, rpkmCountsSel,txsSel ),
                throws_error("Annotation parameter \"gene_name\" is missing!"))
})

test_that("genomic regions correspond to gene counts row names", {
    data("mbSelEnhSignals")
    data("enhSelCoords")
    data("mbSamplesRPKM")
    data("txsSel")

    rpkmCounts2 <- rpkmCountsSel
    rownames(rpkmCounts2)[1] <- "BOO!"

    expect_that(newSigInTAD(enhSel, enhSelGR, rpkmCounts2,txsSel ),
                throws_error("The \"gene_id\" values in GRanges are not expr matrix row names!"))
})




test_that("test log2 gene expression", {
    data("mbSelEnhSignals")
    data("enhSelCoords")
    data("mbSamplesRPKM")
    data("txsSel")

    gID <- "ENSG00000186297.7"
    
    expect_equal( rpkmCountsSel[gID,24], 5.598033 , tolerance=1e-6)

    inTadSig <- newSigInTAD(enhSel, enhSelGR, rpkmCountsSel, txsSel)

    normCounts <- assay(inTadSig@sigMAE[["exprs"]])

    # log2 normalization
    expect_equal( normCounts[gID, 24] , 2.722036,tolerance=1e-6)    



})



