## Correlation testing

context("test on correlation")

test_that("correlation computation is correct", {
    data("mbSelEnhSignals")
    data("enhSelCoords")
    data("mbSamplesRPKM")
    data("txsSel")

    inTadSig <- newSigInTAD(enhSel, enhSelGR, rpkmCountsSel, txsSel)
    inTadSig <- combineInTAD(inTadSig, tadGR)
    cID <- "chr15:26372163-26398073" #E2
 

    corData <- findCorrelation(inTadSig)
    selCorData <- corData[corData$peakid == cID, ]

    expect_equal( nrow(selCorData) , 30 )

    expect_equal( selCorData[ selCorData$name == "GABRA5","cor" ] , 0.878531,tolerance=1e-6)

    # this gene should be in since  option closest is active
    expect_equal( selCorData[ selCorData$name == "GABRB3","cor" ] , 0.6304704, tolerance=1e-6)

})

test_that("correlation properties are working correctly", {
    data("mbSelEnhSignals")
    data("enhSelCoords")
    data("mbSamplesRPKM")
    data("txsSel")

    inTadSig <- newSigInTAD(enhSel, enhSelGR, rpkmCountsSel, txsSel)
    inTadSig <- combineInTAD(inTadSig, tadGR)
    corData <- findCorrelation(inTadSig, method="spearman")
    cID <- "chr15:26372163-26398073" #E2
 
    selCorData <- corData[corData$peakid == cID, ]
    expect_equal( selCorData[ selCorData$name == "GABRA5","cor" ] , 0.8361538, tolerance=1e-6)

    corData <- findCorrelation(inTadSig, adj.pval = TRUE)
    expect_equal( ncol(corData) , 9)

    selCorData <- corData[corData$peakid == cID, ]

    expect_equal( selCorData[ selCorData$name == "GABRA5","qvalue" ] , 6.276324e-06, tolerance=1e-6)

    # novel option 18.11.2018
    expect_equal( selCorData[ selCorData$name == "GABRA5","eucDist" ] , 10.92154 , tolerance=1e-6)



})

test_that("combine intad options are working correctly", {
    data("mbSelEnhSignals")
    data("enhSelCoords")
    data("mbSamplesRPKM")
    data("txsSel")

    inTadSig <- newSigInTAD(enhSel, enhSelGR, rpkmCountsSel, txsSel)
    inTadSig <- combineInTAD(inTadSig, tadGR, closestGene = FALSE)
    corData <- findCorrelation(inTadSig)
    cID <- "chr15:26372163-26398073" #E2
 
    selCorData <- corData[corData$peakid == cID, ]

    expect_equal( nrow(selCorData) , 13)

    expect_equal( selCorData[ selCorData$name == "GABRA5","cor" ] , 0.878531,tolerance=1e-6)

})

test_that("check intad selMaxTadOvlp option", {
    data("mbSelEnhSignals")
    data("enhSelCoords")
    data("mbSamplesRPKM")
    data("txsSel")

    # novel enhnacer overlapping 2 tads
    # std tad chr15:25728907-27128907, close tad chr15:23608559-25368907
    novEnh <- "chr15:25367907-25738907"
    enhSel2 <- rbind(enhSel, c(rep(1, 20),rep(2,5)))
    rownames(enhSel2)[66] <-  novEnh
    enhSelGR2 <- c(enhSelGR,GRanges(novEnh))

    inTadSig <- newSigInTAD(enhSel2, enhSelGR2, rpkmCountsSel, txsSel)

    inTadSig <- filterGeneExpr(inTadSig, checkExprDistr = TRUE)

    # combine genes and signals in TAD: default largest overlap
    inTadSig <- combineInTAD(inTadSig, tadGR)
    expect_equal( length(inTadSig@signalConnections[[novEnh]]$tad) , 8)

    # combine enh with both TADs
    inTadSig <- combineInTAD(inTadSig, tadGR,selMaxTadOvlp = FALSE)
    expect_equal( length(inTadSig@signalConnections[[novEnh]]$tad) , 68)


})







