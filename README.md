## InTAD R package ##

The main target of the package is to detect correlation between gene expression obtained from RNA-seq and selected epigenomic signals i.e. enhancers from ChIP-seq either within topologically associated domains (TADs) or between anchors of chromatin contact loops. Coordinates of TADs/loops can be detected using HiC and other methods and should be provided as additional input.

The analysis within TADs starts from collecting signals and genes lying in the same TAD. Then the combined groups in a TAD are analyzed to detect correlations among them. 

Another analysis approach is focused on the usage of confident chromatin loops. It uses the loops between genomic locations as input and searches if one anchor of the loop is overlapping with gene TSS, while another with target signal peak.

Various analysis parameters can be further controlled. For example, the expression filters, correlation limits (positive,negative), methods (Pearson, Spearman), etc. Each step also provides specific visualization plots.

Package is avialable at BioConductor:
http://bioconductor.org/packages/release/bioc/html/InTAD.html 

For more details please refer to the main vignette and documentation.
