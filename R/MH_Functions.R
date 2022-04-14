#' @title ChromQual
#' @param vcf A vcf file 
#' @param chromlist A vector specifying particular chromosomes
#' @param windowSize Specify window size to calculate number of SNPs
#' @param binwidth Specify bindwidth for histogram plot
#' @return A histogram of of Quality FIELD
#' @examples ChromQual(vcf = "General.vcf", chromlist = c("Chr01","Chr02"),windowsize=1e+05,binwidth=10)
#' @export ChromQual

ChromQual <- function(vcf, chromlist=NULL,windowSize=NULL,binwidth=NULL){
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
message("Making breaks Width")
breaks <- seq(round(min(SNPset$QUAL) - 1, 0), round(max(SNPset$QUAL) + 100, 0), binwidth)
jpeg(file="plot1.jpeg")
message("Plotting histogram")
hist(x = SNPset$QUAL, breaks = breaks, col = "green", xlab ="Quality Quantities", main = "Histogram of Quality Quantities")
     dev.off()
z <- hist(x = SNPset$QUAL, breaks = breaks, col = "green", xlab = "Quality Quantities", main ="Histogram of Quality Quantities")
print(z)
}

#' @title ChromDP
#' @param vcf A vcf file 
#' @param chromlist A vector specifying particular chromosomes
#' @param windowSize Specify window size to calculate number of SNPs
#' @param binwidth Specify bindwidth for histogram plot
#' @return A histogram of of Depth FIELD
#' @examples ChromDP(vcf = "General.vcf", chromlist = c("Chr01","Chr02"),windowsize=1e+05,binwidth=10)
#' @export ChromDP

ChromDP <- function(vcf, chromlist=NULL,windowSize=NULL,binwidth=NULL){
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
  message("Making breaks Width")
  breaks <- seq(round(min(SNPset$DP) - 1, 0), round(max(SNPset$DP) + 100, 0), binwidth)
  jpeg(file="plot2.jpeg")
  message("Plotting histogram")
  hist(x = SNPset$DP, breaks = breaks, col = "green", xlab = "Depth Quantities", main = "Histogram of Depth Quantities")
  dev.off()
  z <- hist(x = SNPset$DP, breaks = breaks, col = "green", xlab ="Depth Quantities", main ="Histogram of Depth Quantities")
  print(z)
}



#' @title ChromRO
#' @param vcf A vcf file 
#' @param chromlist A vector specifying particular chromosomes
#' @param windowSize Specify window size to calculate number of SNPs
#' @param binwidth Specify bindwidth for histogram plot
#' @return A histogram of of Depth FIELD
#' @examples ChromRO(vcf = "General.vcf", chromlist = c("Chr01","Chr02"),windowsize=1e+05,binwidth=10)
#' @export ChromRO

ChromRO <- function(vcf, chromlist=NULL,windowSize=NULL,binwidth=NULL){
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
  message("Making breaks Width")
  breaks <- seq(round(min(SNPset$RO) - 1, 0), round(max(SNPset$RO) + 100, 0), binwidth)
  jpeg(file="plot3.jpeg")
  message("Plotting histogram")
  hist(x = SNPset$RO, breaks = breaks, col = "green", xlab = "RO Quantities", main = "Histogram of RO Quantities")
  dev.off()
  z <- hist(x = SNPset$RO, breaks = breaks, col = "green", xlab ="RO Quantities", main ="Histogram of RO Quantities")
  print(z)
}


#' @title ChromMQM
#' @param vcf A vcf file 
#' @param chromlist A vector specifying particular chromosomes
#' @param windowSize Specify window size to calculate number of SNPs
#' @param binwidth Specify bindwidth for histogram plot
#' @return A histogram of of Depth FIELD
#' @examples ChromMQM(vcf = "General.vcf", chromlist = c("Chr01","Chr02"),windowsize=1e+05,binwidth=10)
#' @export ChromMQM

ChromMQM <- function(vcf, chromlist=NULL,windowSize=NULL,binwidth=NULL){
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
  message("Making breaks Width")
  SNPset$MQM <- as.numeric(SNPset$MQM)
  breaks <- seq(round(min(SNPset$MQM) - 1, 0), round(max(SNPset$MQM) + 100, 0), binwidth)
  jpeg(file="plot4.jpeg")
  message("Plotting histogram")
  hist(x = SNPset$MQM, breaks = breaks, col = "green", xlab = "MQM Quantities", main = "Histogram of MQM Quantities")
  dev.off()
  z <- hist(x = SNPset$MQM, breaks = breaks, col = "green", xlab ="MQM Quantities", main ="Histogram of MQM Quantities")
  print(z)
}



