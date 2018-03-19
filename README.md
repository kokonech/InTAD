## InTAD R package ##

The main target of the package is to detect correlation between gene expression obtained from RNA-seq and selected epigenomic signals i.e. enhancers from ChIP-seq within topologically associated domains (TADs). Coordinates of TADs can be detected using HiC and other methods.

The analysis starts from collecting signals and genes lying in the same TAD. Then the combined groups in a TAD are analyzed to detect correlations among them. Various parameters can be further controlled. For example, the expression filters, correlation limits (positive,negative), methods (Pearson, Spearman), etc. Each step also provides specific visualization plots.

For more details please refer to the main vignette and documentation.
