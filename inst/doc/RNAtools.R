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

## ----regularise----------------------------------------------------------
results <- listRegularise(tested)

head(results[[1]])

## ----p-values------------------------------------------------------------
pval.plots <- listPlotPvals(results, alpha = 0.05)

pval.plots$combined

## ----results-MA----------------------------------------------------------
results.ma <- listResultsMA(results, alpha = 0.05)

results.ma$combined

## ----plotSmear-----------------------------------------------------------
genes.de.names <- rownames(tested$edgeR)[tested$edgeR$FDR < 0.05]
edgeR::plotSmear(normalised$edgeR, de.tags = genes.de.names)

## ----volcano-------------------------------------------------------------
volcanos <- listVolcano(results)

volcanos$combined

## ----voom-volcano--------------------------------------------------------
volcanos$voom

## ----jaccard-------------------------------------------------------------
jaccardTable(results, alpha = 0.05)

## ----venn, fig.height = 7------------------------------------------------
venn <- geneVenn(results, alpha = 0.05)

plot.new()
grid::grid.draw(venn)
text(0.5, 1, "Comparison of Lists", vfont = c("serif", "bold"), cex = 2)
text(0.5, 0, "Number of Differentially Expressed Genes",
     vfont = c("serif", "plain"), cex = 1.5)

## ----venn-gene-----------------------------------------------------------
venn.genes <- vennGenes(results, alpha = 0.05)

sapply(venn.genes, length)

## ----all-genes-----------------------------------------------------------
all.genes <- Reduce(union, venn.genes)

length(all.genes)

## ----gene-summary--------------------------------------------------------
gene.summ <- geneSummary(results, alpha = 0.05)

head(data.frame(gene.summ))

## ----gene-summ-selection-------------------------------------------------
gene.summ.de <- gene.summ[gene.summ$DECount >= 2, ]
gene.summ.de.up <- gene.summ.de[gene.summ.de$meanFC >= 0, ]

head(data.frame(gene.summ.de.up))

## ----gene-sets-----------------------------------------------------------
set.seed(1)
gene.sets <- list(Set1 = sample(gene.summ$Gene, 23),
                  Set2 = sample(gene.summ$Gene, 47),
                  Set3 = sample(gene.summ$Gene, 86),
                  Set4 = sample(gene.summ$Gene, 63),
                  Set5 = sample(gene.summ$Gene, 111),
                  Set6 = sample(gene.summ$Gene, 59))

setTable(results, de.set = gene.summ.de$Gene, gene.sets = gene.sets)

## ----fold-change---------------------------------------------------------
plotFoldChange(data.list = results,
               gene.set = intersect(gene.sets$Set2, gene.summ.de$Gene))

## ----fold-change-counts--------------------------------------------------
counts[c("ENSG00000134802", "ENSG00000131094"), ]

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
    results <- listRegularise(tested)
    
    venn <- geneVenn(results, alpha = 0.05)

    plot.new()
    grid::grid.draw(venn)
    text(0.5, 1, "Comparison of Lists", vfont = c("serif", "bold"), cex = 2)
    text(0.5, 0, "Number of Differentially Expressed Genes",
         vfont = c("serif", "plain"), cex = 1.5)
}

## ----sessionInfo---------------------------------------------------------
sessionInfo()

