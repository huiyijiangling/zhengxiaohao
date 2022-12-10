gdcVolcanoPlot(deg.all, fc = 2, pval = 0.01)

function (deg.all, fc = 2, pval = 0.01) 
{
    geneList <- deg.all
    geneList$threshold <- c()
    geneList$threshold[geneList$logFC > log(fc, 2) & geneList$FDR < 
        pval] <- 1
    geneList$threshold[geneList$logFC >= -log(fc, 2) & geneList$logFC <= 
        log(fc, 2) | geneList$FDR >= pval] <- 2
    geneList$threshold[geneList$logFC < -log(fc, 2) & geneList$FDR < 
        pval] <- 3
    geneList$threshold <- as.factor(geneList$threshold)
    lim <- max(max(geneList$logFC), abs(min(geneList$logFC))) + 
        0.5
    volcano <- ggplot(data = geneList, aes(x = geneList$logFC, 
        y = -log10(geneList$FDR)))
    volcano + geom_point(aes(color = geneList$threshold), alpha = 1, 
        size = 0.8) + xlab("log2(Fold Change)") + ylab("-log10(FDR)") + 
        scale_colour_manual(breaks = geneList$threshold, values = c("red", 
            "black", "green3")) + xlim(c(-lim, lim)) + geom_vline(xintercept = c(-log(fc, 
        2), log(fc, 2)), color = "darkgreen", linetype = 3) + 
        geom_hline(yintercept = -log(pval, 10), color = "darkgreen", 
            linetype = 3) + theme_bw() + theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black"), panel.background = element_blank()) + 
        theme(legend.position = "none") + theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 16))
}
<bytecode: 0x00000184b46b0e48>
<environment: namespace:GDCRNATools>