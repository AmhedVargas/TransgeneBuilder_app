
<h1 align="center">Transgene Builder Shiny App</h1>
<p align="center">
Main repository containing pipelines, data, and software described in: <br><b>An online app for comprehensive <i>C. elegans</i> transgene optimization </b>
</p>

## Abstract
We have developed a "one-stop" online app (www.wormbuilder.org/transgenebuilder/) to adapt transgenes for ad-hoc expression in <i>C. elegans</i>. The app combines commonly used codon adaptation routines and optimization steps into a simple interface. These options include codon optimization for high ubiquitous expression or to minimize germline silencing, removal of piRNA homology or restriction enzyme sites, introduction of synthetic or native introns, and addition of 5' UTRs and 3' UTRs. We compare germline and somatic expression of transgenes generated using different algorithms and show that picking codons based on common tissue expression can result in variable expression. Additionally, we test and validate the use of different short 3' UTRs and endogenous introns that are compatible with  PATC-rich DNA sequences which are known to minimize germline silencing. Finally, since synthetic transgenes are becoming increasingly affordable to produce but conflicting non-coding elements (for example, introns and UTRs) can complicate their gene synthesis, we added extra functionalities to the app which analyze and flag transgenes for potential problems in their synthesis. With this tool, we hope to leverage the time required to implement novel molecular tools yet to be seen in <i>C. elegans</i>.

## Repository description
This github repository contains the code used to acquire the data and perform analysis which led to the development of the [transgene builder app](www.wormbuilder.org/transgenebuilder/). For a complete description of its inner workings and procedure of local installation, check this other [github repository](https://github.com/AmhedVargas/GeneOptimizer_2022).

### Folder structure
This repository is structured into two main folders:
1. Scripts
Shell, R, and Markdown documents.
2. Figures
Final figures after analysis. Please note that very few ended up on the final paper

### Data download
We downloaded data for our analysis from different resources, for example:
* Wormbase
* Ahringer's Lab Reg Atlas
* [piRTABASE](http://cosbi6.ee.ncku.edu.tw/piRTarBase/)

The shell scripts aid
We have deposited DNA and RNA libraries under the [GEO](https://www.ncbi.nlm.nih.gov/geo/) accesion number [GSE165210](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165210)

## Additional information
We have build a [shiny web app](https://wormbuilder.dev/piRNAi/) to ease the construction of piRNAi fragments. 
It's code can be found in this repository: https://github.com/AmhedVargas/piRNAi_dev

Originally, the data was processed using BLAST search as a proxy to identify thew number of mis-matches between a piRNAi guide and any other part of the genome. How the data was processed can be found here: https://github.com/AmhedVargas/piRNAi-DB

However, in the latest version of the piRNAi app we opted to use an exhaustive algorithm based on the calculation of the Hamming distance. Our c++ implementation can be found here: https://github.com/AmhedVargas/CelegansHammingAlignments

## Contact us
Please feel free to send any question or comment to [me](mailto:avargas0lcg@gmail.com) or [Christian](mailto:cfjensen@kaust.edu.sa)



# TransgeneBuilder_app
Source code, analysis workflow, and data analyzed for Vargas-Velazquez et al. manuscript
