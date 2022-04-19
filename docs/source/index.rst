.. VCFstat Package documentation master file, created by
   sphinx-quickstart on Tue Apr 19 13:31:09 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to VCFstat Package's documentation!
===========================================


===============
VCFstat Package
===============

:Author: Michael Hall
:Date:   4/19/2022

VCFstat Package in R
==========

VCFstat is an R package for Downstream analysis of VCF files.

VCFstat Package is still under development and is offered with out any
guarantee.

**For more please read the docs**\ `here <https://github.com/bmansfeld/QTLseqr/raw/master/vignettes/QTLseqr.pdf>`__
---------------------------------------------------------------------------------------------------------------------------------------------



Installation
============

.. raw:: html

You can install VCFstat from github with:

.. code:: r

   # install devtools first to download packages from github
   install.packages("devtools")

   # use devtools to install QTLseqr
   devtools::install_github("PBGLMichaelHall/QTLseqr")

**Note:** QTLseqr is a dependent package so

**If you use QTLseqr in published research, please cite:**

   Mansfeld B.N. and Grumet R, QTLseqr: An R package for bulk segregant
   analysis with next-generation sequencing *The Plant Genome*
   `doi:10.3835/plantgenome2018.01.0006 <https://dl.sciencesocieties.org/publications/tpg/abstracts/11/2/180006>`__

We also recommend citing the paper for the corresponding method you work
with.

QTL-seq method:

   Takagi, H., Abe, A., Yoshida, K., Kosugi, S., Natsume, S., Mitsuoka,
   C., Uemura, A., Utsushi, H., Tamiru, M., Takuno, S., Innan, H., Cano,
   L. M., Kamoun, S. and Terauchi, R. (2013), QTL-seq: rapid mapping of
   quantitative trait loci in rice by whole genome resequencing of DNA
   from two bulked populations. *Plant J*, 74: 174–183.
   `doi:10.1111/tpj.12105 <https://onlinelibrary.wiley.com/doi/full/10.1111/tpj.12105>`__

G prime method:

   Magwene PM, Willis JH, Kelly JK (2011) The Statistics of Bulk
   Segregant Analysis Using Next Generation Sequencing. *PLOS
   Computational Biology* 7(11): e1002255.
   `doi.org/10.1371/journal.pcbi.1002255 <http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002255>`__

Abstract
--------

Examples:
=========

Load/install libraries
======================

.. code:: r 
   
   install.packages(“tinytex”) 
   install.packages(“vcfR”) 
   install.packages(“tidyr”) 
   install.packages(“ggplot2”)
   devtools::install_github(“PBGLMichaelHall/QTLseqr”,force = TRUE)
   devtools::install_github("PBGLMichaelHall/VCFstat",force = TRUE)   
   library(QTLseqr) 
   library(VCFstat)
   library(tinytex) 
   library(vcfR) 
   library(tidyr)
   library(ggplot2)

::

   # Set the Working Directory to where VCF file is stored in computer file system

.. code:: r

   setwd("/home/michael/Desktop/QTLseqr/extdata/")
   chromlist <- c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10")


   VCFstat::ChromQual(vcf = "freebayes_D2.filtered.vcf.gz", chromlist = chromlist,windowSize = 1e+05, binwidth = 50,Maximum=5000)

.. figure:: ../images/1.png
   :alt: 

.. code:: r

   VCFstat::ChromDP(vcf = "freebayes_D2.filtered.vcf.gz", chromlist = chromlist,windowSize = 1e+05, binwidth = 10)

.. figure:: ../images/2.png

  

.. code:: r

   
   VCFstat::ChromRO(vcf = "freebayes_D2.filtered.vcf.gz", chromlist = chromlist,windowSize = 1e+05, binwidth = 10)
.. figure:: ../images/3.png
   :alt: 


.. code:: r

   VCFstat::ChromAO(vcf = "freebayes_D2.filtered.vcf.gz", chromlist = chromlist,windowSize = 1e+05, binwidth = 10)

.. figure:: ../images/4.png


.. code:: r

   VCFstat::ChromAO(vcf = "freebayes_D2.filtered.vcf.gz", chromlist = chromlist,windowSize = 1e+05, binwidth = 10)

.. figure:: ../images/5.png


.. code:: r

   VCFstat::ChromAC(vcf = "freebayes_D2.filtered.vcf.gz", chromlist = chromlist,windowSize = 1e+05, binwidth = 10)

.. figure:: ../images/6.png
   :alt: 

.. code:: r

   VCFstat::ChromAN(vcf = "freebayes_D2.filtered.vcf.gz", chromlist = chromlist,windowSize = 1e+05, binwidth = 10)

.. figure:: ../images/7.png
   :alt: 


.. code:: r

   VCFstat::ChromnSNPs(vcf = "freebayes_D2.filtered.vcf.gz", chromlist = chromlist,windowSize = 1e+05, binwidth = 10)

.. figure:: ../images/8.png
   :alt: 

.. code:: r

   VCFstat::FacetChromnSNPs(vcf = "freebayes_D2.filtered.vcf.gz", chromlist = chromlist,windowSize = 1e+05, ncol=10)

.. figure:: ../images/9.png
  

.. code:: r

   VCFstat::FacetChromQual(vcf = "freebayes_D2.filtered.vcf.gz", chromlist = chromlist,windowSize = 1e+05, ncol=10)

.. figure:: ../images/10.png
   :alt: 


.. code:: r

   VCFstat::FacetChromDP(vcf = "freebayes_D2.filtered.vcf.gz", chromlist = chromlist,windowSize = 1e+05, ncol=10)

.. figure:: ../images/11.png
   :alt: 

.. code:: r
   
   VCFstat::FacetChromAO(vcf = "freebayes_D2.filtered.vcf.gz", chromlist = chromlist,windowSize = 1e+05, ncol=10)

.. figure:: ../images/12.png
   :alt: 

.. figure:: ../images/13.png

.. code:: r

   setwd("/home/michael/Desktop/Variants/Decompressed/")

   chromlist <- c("NC_016131.3","NC_016132.3","NC_016133.3","NC_016134.3","NC_016135.3")


   VCFstat::ChromQual(vcf = "freebayes~bwa~GCF_000005505.3_Brachypodium_distachyon_v3.0~all_samples~filtered-strict.vcf", chromlist = chromlist,windowSize = 1e+05, binwidth = 50,Maximum=5000)

.. figure:: ../images/14.png
   :alt:



.. code:: r

   VCFstat::ChromDP(vcf = "freebayes~bwa~GCF_000005505.3_Brachypodium_distachyon_v3.0~all_samples~filtered-strict.vcf", chromlist = chromlist,windowSize = 1e+05, binwidth = 10)

.. figure:: ../images/15.png
   :alt: 

 

.. code:: r

   VCFstat::ChromRO(vcf = "freebayes~bwa~GCF_000005505.3_Brachypodium_distachyon_v3.0~all_samples~filtered-strict.vcf", chromlist = chromlist,windowSize = 1e+05, binwidth = 10)

.. figure:: ../images/16.png
   :alt: 

  
.. code:: r

   VCFstat::ChromAO(vcf = "freebayes~bwa~GCF_000005505.3_Brachypodium_distachyon_v3.0~all_samples~filtered-strict.vcf", chromlist = chromlist,windowSize = 1e+05, binwidth = 10)

.. figure:: ../images/17.png
   :alt: 

 .. code:: r

   VCFstat::ChromMQM(vcf = "freebayes~bwa~GCF_000005505.3_Brachypodium_distachyon_v3.0~all_samples~filtered-strict.vcf", chromlist = chromlist,windowSize = 1e+05, binwidth = 10)

.. figure:: ../images/18.png
   :alt: 

   .. code:: r

   VCFstat::ChromAC(vcf = "freebayes~bwa~GCF_000005505.3_Brachypodium_distachyon_v3.0~all_samples~filtered-strict.vcf", chromlist = chromlist,windowSize = 1e+05, binwidth = 10)

.. figure:: ../images/19.png
   :alt:

.. code:: r

   VCFstat::ChromAN(vcf = "freebayes~bwa~GCF_000005505.3_Brachypodium_distachyon_v3.0~all_samples~filtered-strict.vcf", chromlist = chromlist,windowSize = 1e+05, binwidth = 10)

.. figure:: ../images/20.png
   :alt: 

.. code:: r

   VCFstat::ChromnSNPs(vcf = "freebayes~bwa~GCF_000005505.3_Brachypodium_distachyon_v3.0~all_samples~filtered-strict.vcf", chromlist = chromlist,windowSize = 1e+05, binwidth = 10)

.. figure:: ../images/21.png
   :alt: 

 

.. code:: r

   VCFstat::FacetChromnSNPs(vcf = "freebayes~bwa~GCF_000005505.3_Brachypodium_distachyon_v3.0~all_samples~filtered-strict.vcf", chromlist = chromlist,windowSize = 1e+05, ncol=10)

.. figure:: ../images/22.png
   :alt:
   
.. code:: r

   VCFstat::FacetChromQual(vcf = "freebayes~bwa~GCF_000005505.3_Brachypodium_distachyon_v3.0~all_samples~filtered-strict.vcf", chromlist = chromlist,windowSize = 1e+05, ncol=10)

.. figure:: ../images/23.png
   :alt: 

.. code:: r

   VCFstat::FacetChromDP(vcf = "freebayes~bwa~GCF_000005505.3_Brachypodium_distachyon_v3.0~all_samples~filtered-strict.vcf", chromlist = chromlist,windowSize = 1e+05, ncol=10)
.. figure:: ../images/24.png
   :alt: 

.. code:: r

   VCFstat::FacetChromAO(vcf = "freebayes~bwa~GCF_000005505.3_Brachypodium_distachyon_v3.0~all_samples~filtered-strict.vcf", chromlist = chromlist,windowSize = 1e+05, ncol=10)

.. figure:: ../images/25.png
   :alt: 

.. figure:: ../images/26.png



