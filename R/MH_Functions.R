# Define global variables otherwise Build check will complain and make a note
globalVariables(c("MQM","AO","CHROM","DP","POS","QUAL","aes","facet_wrap","geom_density","geom_histogram","ggplot","nSNPs","theme_classic","%>%","rbindlist","read.vcfR","vcfR2tidy","countSNPs_cpp"))

#' @importFrom dplyr %>% 
#' @import graphics
#' @import vcfR
#' @import data.table
#' @import dplyr
#' @import QTLseqr

#' @title ChromQual
#' @param vcf A vcf file 
#' @param chromlist A vector specifying particular chromosomes
#' @param windowSize Specify window size to calculate number of SNPs
#' @return A histogram of of Quality FIELD
#' @export ChromQual

ChromQual <- function(vcf, chromlist=NULL,windowSize=NULL){
vcf <- read.vcfR(file = vcf)
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

jpeg(filename="plot1.jpeg")
message("Plotting histogram")
hist(x = SNPset$QUAL, breaks = "sturges", col = "green", xlab ="Quality Quantities", main = "histogram of Quality Quantities")
     dev.off()
z <- hist(x = SNPset$QUAL, breaks = "sturges", col = "green", xlab = "Quality Quantities", main ="histogram of Quality Quantities")
print(z)
}

#' @title ChromDP
#' @param vcf A vcf file 
#' @param chromlist A vector specifying particular chromosomes
#' @param windowSize Specify window size to calculate number of SNPs
#' @return A histogram of of Depth FIELD
#' @export ChromDP

ChromDP <- function(vcf, chromlist=NULL,windowSize=NULL){
  vcf <- read.vcfR(file = vcf)
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
  jpeg(filename="plot2.jpeg")
  message("Plotting histogram")
  hist(x = SNPset$DP, breaks = "sturges", col = "green", xlab = "Depth Quantities", main = "histogram of Depth Quantities")
  dev.off()
  z <- hist(x = SNPset$DP, breaks = "sturges", col = "green", xlab ="Depth Quantities", main ="histogram of Depth Quantities")
  print(z)
}



#' @title ChromRO
#' @param vcf A vcf file 
#' @param chromlist A vector specifying particular chromosomes
#' @param windowSize Specify window size to calculate number of SNPs
#' @return A histogram of of Depth FIELD
#' @export ChromRO

ChromRO <- function(vcf, chromlist=NULL,windowSize=NULL){
  vcf <- read.vcfR(file = vcf)
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
  jpeg(filename="plot3.jpeg")
  message("Plotting histogram")
  hist(x = SNPset$RO, breaks = "sturges", col = "green", xlab = "RO Quantities", main = "histogram of RO Quantities")
  dev.off()
  z <- hist(x = SNPset$RO, breaks = "sturges", col = "green", xlab ="RO Quantities", main ="histogram of RO Quantities")
  print(z)
}


#' @title ChromAO
#' @param vcf A vcf file 
#' @param chromlist A vector specifying particular chromosomes
#' @param windowSize Specify window size to calculate number of SNPs
#' @return A histogram of of Depth FIELD
#' @export ChromAO

ChromAO <- function(vcf, chromlist=NULL,windowSize=NULL){
  vcf <- read.vcfR(file = vcf)
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
  SNPset$AO <- as.numeric(SNPset$AO)
  jpeg(filename="plot4.jpeg")
  message("Plotting histogram")
  hist(x = SNPset$AO, breaks = "sturges", col = "green", xlab = "AO Quantities", main = "histogram of AO Quantities")
  dev.off()
  z <- hist(x = SNPset$AO, breaks = "sturges", col = "green", xlab ="AO Quantities", main ="histogram of AO Quantities")
  print(z)
}


#' @title ChromMQM
#' @param vcf A vcf file 
#' @param chromlist A vector specifying particular chromosomes
#' @param windowSize Specify window size to calculate number of SNPs
#' @return A histogram of of Depth FIELD
#' @export ChromMQM

ChromMQM <- function(vcf, chromlist=NULL,windowSize=NULL){
  vcf <- read.vcfR(file = vcf)
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
  jpeg(filename="plot5.jpeg")
  message("Plotting histogram")
  hist(x = SNPset$MQM, breaks = "sturges", col = "green", xlab = "MQM Quantities", main = "histogram of MQM Quantities")
  dev.off()
  z <- hist(x = SNPset$MQM, breaks = "sturges", col = "green", xlab ="MQM Quantities", main ="histogram of MQM Quantities")
  print(z)
}


#' @title ChromAC
#' @param vcf A vcf file 
#' @param chromlist A vector specifying particular chromosomes
#' @param windowSize Specify window size to calculate number of SNPs
#' @return A histogram of of Depth FIELD
#' @export ChromAC

ChromAC <- function(vcf, chromlist=NULL,windowSize=NULL){
  vcf <- read.vcfR(file = vcf)
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
  SNPset$AC <- as.factor(SNPset$AC)
  jpeg(filename="plot6.jpeg")
  message("Plotting barplot of AC")
  plot(SNPset$AC)
  dev.off()
  totals <- table(SNPset$AC)
  print(totals)
  z <- plot(SNPset$AC)
  print(z)

}

#' @title ChromAN
#' @param vcf A vcf file 
#' @param chromlist A vector specifying particular chromosomes
#' @param windowSize Specify window size to calculate number of SNPs
#' @return A histogram of of Depth FIELD
#' @export ChromAN

ChromAN <- function(vcf, chromlist=NULL,windowSize=NULL){
  vcf <- read.vcfR(file = vcf)
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
  SNPset$AC <- as.factor(SNPset$AN)
  jpeg(filename="plot7.jpeg")
  message("Plotting barplot of AN")
  plot(SNPset$AC)
  dev.off()
  totals <- table(SNPset$AN)
  print(totals)
  z <- plot(SNPset$AN)
  print(z)
  
}

#' @title ChromnSNPs
#' @param vcf A vcf file please
#' @param chromlist A vector specifying particular chromosomes
#' @param windowSize Specify window size to calculate number of SNPs
#' @return A histogram of of Depth FIELD
#' @export ChromnSNPs

ChromnSNPs <- function(vcf, chromlist=NULL,windowSize=NULL){
  vcf <- read.vcfR(file = vcf)
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
  jpeg(filename="plot8.jpeg")
  message("Plotting histogram")
  hist(x = SNPset$nSNPs, breaks = "sturges", col = "green", xlab = "nSNPs Quantities", main = "histogram of nSNPs Quantities")
  dev.off()
  z <- hist(x = SNPset$nSNPs, breaks = "sturges", col = "green", xlab ="nSNPs Quantities", main ="histogram of nSNPs Quantities")
  print(z)
  
}



#' @title FacetChromnSNPs
#' @param vcf A vcf file please
#' @param chromlist A vector specifying particular chromosomes
#' @param windowSize Specify window size to calculate number of SNPs
#' @param ncol An integer representing the number of Chromosomes in your set list
#' @export FacetChromnSNPs

FacetChromnSNPs <- function(vcf, chromlist=NULL,windowSize=NULL,ncol=NULL){
  vcf <- read.vcfR(file = vcf)
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
  SNPset$nSNPs <- as.numeric(SNPset$nSNPs)
  jpeg(filename="plot9.jpeg")
  ggplot(data = SNPset, aes(x = nSNPs)) + geom_histogram(bins = 10, show.legend = TRUE) + facet_wrap(~CHROM, ncol = ncol) + theme_classic()
  dev.off()
  z<-  ggplot(data = SNPset, aes(x = nSNPs)) + geom_histogram(bins = 10, show.legend = TRUE) + facet_wrap(~CHROM, ncol = ncol) + theme_classic()
  print(z)
  
}



#' @title FacetChromQual
#' @param vcf A vcf file please
#' @param chromlist A vector specifying particular chromosomes
#' @param windowSize Specify window size to calculate number of SNPs
#' @param ncol An integer representing the number of Chromosomes in your set list
#' @export FacetChromQual

FacetChromQual <- function(vcf, chromlist=NULL,windowSize=NULL,ncol=NULL){
  vcf <- read.vcfR(file = vcf)
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
  jpeg(filename="plot10.jpeg")
  ggplot(data = SNPset, aes(x = QUAL)) + geom_histogram(bins = 10, show.legend = TRUE) + facet_wrap(~CHROM, ncol = ncol) + theme_classic()
  dev.off()
  z<-  ggplot(data = SNPset, aes(x = QUAL)) + geom_histogram(bins = 10, show.legend = TRUE) + facet_wrap(~CHROM, ncol = ncol) + theme_classic()
  print(z)
  
}


#' @title FacetChromDP
#' @param vcf A vcf file please
#' @param chromlist A vector specifying particular chromosomes
#' @param windowSize Specify window size to calculate number of SNPs
#' @param ncol An integer representing the number of Chromosomes in your set list
#' @export FacetChromDP

FacetChromDP <- function(vcf, chromlist=NULL,windowSize=NULL,ncol=NULL){
  vcf <- read.vcfR(file = vcf)
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
  SNPset$DP <- as.numeric(SNPset$DP)
  jpeg(filename="plot11.jpeg")
  ggplot(data = SNPset, aes(x = DP)) + geom_histogram(bins = 10, show.legend = TRUE) + facet_wrap(~CHROM, ncol = ncol) + theme_classic()
  dev.off()
  z<-  ggplot(data = SNPset, aes(x = DP)) + geom_histogram(bins = 10, show.legend = TRUE) + facet_wrap(~CHROM, ncol = ncol) + theme_classic()
  print(z)
  
}

#' @title FacetChromAO
#' @param vcf A vcf file please please
#' @param chromlist A vector specifying particular chromosomes
#' @param windowSize Specify window size to calculate number of SNPs
#' @param ncol An integer representing the number of Chromosomes in your set list
#' @export FacetChromAO

FacetChromAO <- function(vcf, chromlist=NULL,windowSize=NULL,ncol=NULL){
  vcf <- read.vcfR(file = vcf)
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
  SNPset$AO <- as.numeric(SNPset$AO)
  jpeg(filename="plot9.jpeg")
  ggplot(data = SNPset, aes(x = AO)) + geom_density() + facet_wrap(~CHROM, ncol = ncol) + theme_classic()
  dev.off()
  jpeg(filename="plot12.jpeg")
  ggplot(data = SNPset, aes(x = AO)) + geom_histogram(stat="count") + facet_wrap(~CHROM, ncol = ncol) + theme_classic()
  dev.off()
  z<-  ggplot(data = SNPset, aes(x = AO)) + geom_density() + facet_wrap(~CHROM, ncol = ncol) + theme_classic()
  print(z)
  z1<-ggplot(data = SNPset, aes(x = AO)) + geom_histogram(stat="count") + facet_wrap(~CHROM, ncol = ncol) + theme_classic()
  print(z1)
  
}
#' @title Correlation
#' @description Returns a Correlation Matrix for info fields AC, DP, DPB, QA, RO, AO and most importantly QUAL
#' @param vcffile A vcf file
#' @param chromlist A vector of chromosomes/contigs as specified in VCF file
#' @param p1 A boolean Value either TRUE or FALSE for plotting Correlation matrix and Correlation Tables
#' @param p2 A boolean Value either TRUE or FALSE for plotting Correlation matrix and Correlation Tables
#' @param p3 A boolean Value either TRUE or FALSE for plotting Correlation matrix and Correlation Tables
#' @param p4 A boolean Value either TRUE or FALSE for plotting Correlation matrix and Correlation Tables
#' @param p5 A boolean Value either TRUE or FALSE
#' @export Correlation


Correlation <-
  function (vcffile = NULL, chromlist = NULL,p1 = NULL, p2 = NULL, p3 = NULL, p4 = NULL,p5=TRUE)
  {
    vcf <- read.vcfR(file = vcffile)
    vcf <- vcfR2tidy(vcf)
    message("Extracting unique Chromosome or Contig names reverse compatible to VCF file")
    print(unique(vcf$fix$CHROM))
    SNPset <- vcf
    SNPset <- Map(as.data.frame, SNPset)
    SNPset <- rbindlist(SNPset, fill = TRUE)
    if (!is.null(chromlist)) {
      message("Preparing Data for Quality Control Plotting: ",
              paste(unique(SNPset$CHROM)[!unique(SNPset$CHROM) %in%
                                           chromlist], collapse = ", "))
      SNPset <- SNPset[SNPset$CHROM %in% chromlist,]
      message("Finishing Chromosome Subset")
    }
    
    message("Factoring Chromosome Variable According to Unique Specification")
    SNPset$CHROM <- factor(SNPset$CHROM, levels = gtools::mixedsort(unique(SNPset$CHROM)))
    
    message("Selecting Variable Subset")
    SNPset <- SNPset %>% select(QUAL,AC,DP,DPB,QA,RO,AO)
    SNPset <- as.data.frame(sapply(SNPset, as.numeric))
    message("Mutating SNPS set creating nSNPs variable")
    p1 <- p1
    if (p1 == TRUE){
      t2 <- cor(SNPset)
      round(t2,2)
      print(t2)
    }else if (p1 == FALSE){
      print("Do not plot cor(SNPset)")
    }
    p2 <- p2
    if (p2 == TRUE){
      t3 <- rcorr(as.matrix(SNPset))
      print(t3)
    }else if (p2 == FALSE){
      print("Do not plot rcorr(as.matrix(SNPset))")
    }
    p3 <- p3
    if (p3 == TRUE){
      t4 <- corrplot(cor(SNPset))
      print(t4)
    }else if (p3 == FALSE){
      print("Do not plot corrplot(cor(SNPset))")
    }
    p4 <- p4
    if (p4 == TRUE){
      t5 <- chart.Correlation(SNPset,histogram=TRUE,pch=19)
      print(t5)
    }else if (p4 == FALSE){
      print("Do not plot chart.Correlation(SNPset,histogram=TRUE,pch=19)")
    }
    p5 <- p5
    if (p5 == TRUE){
      col <- colorRampPalette(c("blue","white","red","purple"))(20)
      t2 <- cor(SNPset)
      heatmap(t2,col=col,symm=TRUE)
    }else if (p5 == FALSE){
      print("Do not plot heatmap")
    }
    message("Returning completed Data frame as a SNPSet")
    return(as.data.frame(SNPset))
  }


