#' @title ChromQual
#' @param vcf A vcf file 
#' @param chromlist A vector specifying particular chromosomes
#' @param windowSize Specify window size to calculate number of SNPs
#' @param binwidth Specify bindwidth for histogram plot
#' @param FIELD Specify info field to plot
#' @return A histogram of of FIELD
#' @examples ChromQuality(vcf = "General.vcf", chromlist = c("Chr1", "Chr2")), windowSize = 1e+06, scalar = 0.1, ncol = 2,HighLimQuality = 6000,  binwidth1 = 100, binwidth2 =1, DPBINS=10, p1=TRUE, p2=FALSE, p3=TRUE, p4=TRUE, p5=FALSE, p6=TRUE)
#' @export ChromQual

ChromQual <- function(vcf, chromlist=NULL,windowSize=NULL,binwidth=NULL,FIELD=NULL){
vcf <- read.vcfR(file = "freebayes_D2.filtered.vcf.gz")
vcf <- vcfR2tidy(vcf)
SNPset <- vcf
SNPset <- Map(as.data.frame, SNPset)
SNPset <- rbindlist(SNPset, fill = TRUE)
if (!is.null(chromlist)) {
  message("Preparing Data for Quality Control Plotting and removing the following Chromosomes/Contigs: ", 
          paste(unique(SNPset$CHROM)[!unique(SNPset$CHROM) %in% 
                                       chromlist], collapse = ", "))
  SNPset <- SNPset[SNPset$CHROM %in% chromlist, ]
  message("Finishing Chromosome Subset")
}
message("Factoring Chromosome Variable According to Unique Specification")
SNPset$CHROM <- factor(SNPset$CHROM, levels = gtools::mixedsort(unique(SNPset$CHROM)))
message("Selecting Variable Subset")
SNPset <- SNPset %>% dplyr::group_by(CHROM) %>% dplyr::mutate(nSNPs = countSNPs_cpp(POS = POS, windowSize = windowSize))
par(mfrow = c(1, 1))

breaks <- seq(round(min(SNPset$FIELD) - 1, 0), round(max(SNPset$FIELD) + 100, 0), binwidth)
jpeg(file="plot1.jpeg")
hist(x = SNPset$FIELD, breaks = breaks, col = "green", frequency = TRUE, xlab = paste0(FIELD,"Quantities"), main = paste0("Histogram of",FIELD))
     dev.off()
}