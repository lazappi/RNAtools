## ----global_options, include = FALSE-------------------------------------
# Default RMarkdown options
# Changing these will change them for all chunks unless otherwise set
knitr::opts_chunk$set(
                        autodep        = TRUE,
                        cache          = TRUE,
                        cache.comments = TRUE,
                        collapse       = TRUE,
                        comment        = "##",
                        #dev            = "png",
                        echo           = TRUE,
                        error          = FALSE,
                        fig.width      = 7,
                        fig.height     = 5,
                        fig.align      = "center",
                        highlight      = TRUE,
                        include        = TRUE,
                        message        = FALSE,
                        prompt         = FALSE,
                        results        = "markup",
                        size           = "normalsize",
                        strip.white    = TRUE,
                        tidy           = FALSE,
                        tidy.opts      = NULL,
                        warning        = FALSE
               )

## ----libraries-----------------------------------------------------------
library("lazappiRNA")

## ----data----------------------------------------------------------------
library("HTSFilter")

data("sultan")

counts <- exprs(sultan)
groups <- pData(sultan)$cell.line

rm(sultan)

groups

head(counts)

## ----count-density-------------------------------------------------------
countDensity(counts)

## ----count-density-bw----------------------------------------------------
countDensity(counts) + ggplot2::theme_bw() + ggplot2::ggtitle("Count Densities")

## ----transform-----------------------------------------------------------
transformed <- transformCounts(counts,
                               methods = c("log", "vst", "rlog", "logCPM"))

names(transformed)

head(transformed$log)

## ----densities-----------------------------------------------------------
densities <- listDensity(transformed)

names(densities)

densities$combined

## ----boxplots------------------------------------------------------------
boxplots <- listBoxplots(transformed)

boxplots$combined

## ----ma-plots------------------------------------------------------------
ma.plots <- listCountMA(transformed)

ma.plots$combined

## ----single-ma-----------------------------------------------------------
countMA(transformed$log)

## ----heatmaps, fig.height = 7--------------------------------------------
heatmap <- countHeatmap(transformed$log, groups = groups)

names(heatmap)

showHeatmap(heatmap)

## ----PCA-----------------------------------------------------------------
PCA <- listPCA(transformed, top = 500, groups = groups)

PCA$combined

## ----MDS-----------------------------------------------------------------
MDS <- listMDS(transformed, top = 500, group = groups, selection = "pairwise")

MDS$combined

## ----objects-------------------------------------------------------------
objects <- counts2Objects(counts, groups, filter = TRUE,
                          methods = c("edgeR", "DESeq", "DESeq2", "voom"))

names(objects)

for(object in objects) {
    print(class(object))
}

## ----filtering-----------------------------------------------------------
nrow(counts)

nrow(objects$voom$counts)

## ----normalise-----------------------------------------------------------
normalised <- listNorm(objects)

## ----test----------------------------------------------------------------
group1 <- "HEK293T"
group2 <- "Ramos B cell"

tested <- listTest(normalised, group1, group2, filter = TRUE)

## ----p-values------------------------------------------------------------
pval.plots <- listPlotPvals(tested, alpha = 0.05)

pval.plots$combined

## ----results-MA----------------------------------------------------------
results.ma <- listResultsMA(tested, alpha = 0.05)

results.ma$combined

## ----plotSmear-----------------------------------------------------------
genes.de.names <- rownames(tested$edgeR)[tested$edgeR$FDR < 0.05]
edgeR::plotSmear(normalised$edgeR, de.tags = genes.de.names)

rm(genes.de.names)

## ----volcano-------------------------------------------------------------
volcanos <- listVolcano(tested)

volcanos$combined

## ----voom-volcano--------------------------------------------------------
volcanos$voom

## ----jaccard-------------------------------------------------------------
jaccardTable(tested, alpha = 0.05)

## ----venn, fig.height = 7------------------------------------------------
venn <- geneVenn(tested, alpha = 0.05)

plot.new()
grid::grid.draw(venn)
text(0.5, 1, "Comparison of Lists", vfont = c("serif", "bold"), cex = 2)
text(0.5, 0, "Number of Differentially Expressed Genes",
     vfont = c("serif", "plain"), cex = 1.5)

## ----venn-gene-----------------------------------------------------------
venn.genes <- vennGenes(tested, alpha = 0.05)

sapply(venn.genes, length)

## ----all-genes-----------------------------------------------------------
all.genes <- Reduce(union, venn.genes)

length(all.genes)

## ----pasilla, fig.height = 7---------------------------------------------
if (requireNamespace("pasilla", quietly = TRUE) &&
        requireNamespace("Biobase", quietly = TRUE)) {
    
    data("pasillaGenes", package = "pasilla")
    
    group1 <- "untreated"
    group2 <- "treated"
    
    counts <- counts(pasillaGenes)
    groups <- factor(pData(pasillaGenes)[, c("condition")])
    groups <- relevel(groups, ref = "untreated")
    
    objects <- counts2Objects(counts, groups = groups, 
                              methods = c("edgeR", "DESeq", "DESeq2", "voom"),
                              filter = TRUE)
    normalised <- listNorm(objects)
    tested <- listTest(normalised, group1 = group1, group2 = group2,
                       filter = TRUE)
    
    venn <- geneVenn(tested, alpha = 0.05)

    plot.new()
    grid::grid.draw(venn)
    text(0.5, 1, "Comparison of Lists", vfont = c("serif", "bold"), cex = 2)
    text(0.5, 0, "Number of Differentially Expressed Genes",
         vfont = c("serif", "plain"), cex = 1.5)
}

## ----sessionInfo---------------------------------------------------------
sessionInfo()

