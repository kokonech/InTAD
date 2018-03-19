## Gene expression filtering

context("test filtering of gene expression")

test_that("genes with no expression are filtered", {
    data("mbSelEnhSignals")
    data("enhSelCoords")
    data("mbSamplesRPKM")
    data("txsSel")

    inTadSig <- newSigInTAD(enhSel, enhSelGR, rpkmCountsSel, txsSel)

    inTadSig <- filterGeneExpr(inTadSig)
    expect_equal( nrow(assay(inTadSig@sigMAE[["exprs"]])) , 1898 )

})

test_that("transcritpt type missing mark", {
    data("mbSelEnhSignals")
    data("enhSelCoords")
    data("mbSamplesRPKM")
    data("txsSel")
    
    inTadSig <- newSigInTAD(enhSel, enhSelGR, rpkmCountsSel, txsSel)
    
    expect_error( filterGeneExpr(inTadSig, geneType = "boo"),  
                    "The gene type \"boo\" is not found!" )

})

test_that("transcritpt type missing mark", {
    data("mbSelEnhSignals")
    data("enhSelCoords")
    data("mbSamplesRPKM")
    data("txsSel")
    
    inTadSig <- newSigInTAD(enhSel, enhSelGR, rpkmCountsSel, txsSel)
    inTadSig <- filterGeneExpr(inTadSig, geneType = "protein_coding") 

    expect_equal(  nrow(assay(inTadSig@sigMAE[["exprs"]])) , 612 )

})

# problematic test due to issue with mclust

test_that("transcritpt type missing mark", {
    data("mbSelEnhSignals")
    data("enhSelCoords")
    data("mbSamplesRPKM")
    data("txsSel")
    
    inTadSig <- newSigInTAD(enhSel, enhSelGR, rpkmCountsSel, txsSel)
    inTadSig <- filterGeneExpr(inTadSig, geneType = "protein_coding", checkExprDistr  = TRUE) 

    # update 19.02 with gene_id: 
    # should be 429 by default for protein_coding only
    expect_equal(  nrow(assay(inTadSig@sigMAE[["exprs"]])) , 429 )

})




