# ReViewCNV

![GitHub release](https://img.shields.io/github/release/ctglab/ReViewCNV.svg)
[![R-CMD-check](https://github.com/ctglab/ReViewCNV/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ctglab/ReViewCNV/actions/workflows/R-CMD-check.yaml)
![GitHub docker image build](https://github.com/ctglab/ReViewCNV/actions/workflows/docker-image.yml/badge.svg)
![GitHub commit activity](https://img.shields.io/github/commit-activity/m/ctglab/ReViewCNV) 
![GitHub last commit](https://img.shields.io/github/last-commit/ctglab/ReViewCNV)
![GitHub contributors](https://img.shields.io/github/contributors/ctglab/ReViewCNV)
[![GitHub Website](https://img.shields.io/website-up-down-green-red/http/monip.org.svg)](https://ctglab.github.io/)
![GitHub forks](https://img.shields.io/github/forks/ctglab/ReViewCNV)
![GitHub Repo stars](https://img.shields.io/github/stars/ctglab/ReViewCNV)
![GitHub watchers](https://img.shields.io/github/watchers/ctglab/ReViewCNV.svg)


ReViewCNV is a  Shiny App for the visualization of germinal Copy Number Variants (CNVs). It accepts as input [EXCAVATOR2](https://pubmed.ncbi.nlm.nih.gov/27507884/) output files, MIXER output files or a list of CNVs identified by any bioinformatic tool with three columns (without header) indicating for each CNV the genomic coordinates (Chromosome, Start and End) and optionally a fourth column specifying the CNV type (i.e. deletion or duplication). Example input files for each of the three options are available in the EXAMPLE folder.
</p>
<p align="justify">
The App allows the user to compare the  CNVs of interst with population polymorphisms present in different public datasets (AnnotSV, DGV and gnomAD) and to identify the genes present in correspondence of the altered regions (genome-wide). The exons of the genes present in correspondance of the CNVs of interest are highlighted. It is possible to visualize up to three individuals at the same time in synchronized plots, facilitating family studies and the identification of de novo mutations.
</p>
<p align="justify">
You can install the development version of AppPackage from GitHub with:

``` r
# install.packages("pak")
pak::pak("ctglab/ReViewCNV")
```

To run the app just type
``` r
ReViewCNV()
```

  
It is also possible to install the Dockerized version of the app (see below). 
</p>
<p align="justify">
The ShinyApp has been developed in R (v. 4.4.5) using RStudio as IDE and the following R libraries: arrow, bslib, dplyr, htmlwidgets,plotly, shiny, shinyHugePlot and stringr. The base image used for the Dockerfile is rhub/rminimal.
</p>

<br/>


### Building and running the Docker


##### To build the image, download the repository, open the terminal, go to the folder with the downloaded files and run the command below (this will take around 25 minutes):

_sudo docker buildx build  -t ReViewCNV ._

##### To run the Docker use:

_sudo docker run --name ReViewCNV  -p 6868:6868 ReViewCNV Rscript ReViewCNV.R_


<br/>

### Building and running the Singularity/Apptainer

#####  To create the .sif file it is necessary to build the Docker image first (see above). The Docker image should be saved locally as a .tar file using the command below:

_sudo docker save ReViewCNV > ReViewCNV.tar_ 

##### The .tar file  can be converted to a .sif file using one of the two commands below:

_singularity build ReViewCNV.sif docker-archive://ReViewCNV.tar_ (if using singularity)

_apptainer build ReViewCNV.sif docker-archive://ReViewCNV_ (if using apptainer)

##### To run the the app  on apptainer/singularity use one of the two commands below:

_singularity run  ReViewCNV.sif  Rscript ReViewCNV.R_ (if using singularity)

_apptainer run  ReViewCNV.sif  Rscript ReViewCNV.R_ (if using apptainer)



### How to use the app
<p align="justify">
The animated gifs below illustrate how the app works. The first one shows how to upload the input files in the app_CNV_tsv.R. The other gifs show the app functionalities in app_Excavator2.R, they are almost the same in  app_CNV_tsv.R.
</p>

<br/>

### Questions?
Do you have any question or suggestion? Please [let us know](https://github.com/Francesco85P), we look forward to your feedback!


<br/>

#### Uploading the input file in the app_CNV_tsv.R 

![App_CNV_tsv_upload](https://github.com/user-attachments/assets/36d92350-09a8-49ca-ab27-5ce25b7f606f)



<br/><br/>


#### Specifying the genome assembly and interacting with the identified regions and the population polymorphisms plots

![Video_1b](https://github.com/ctglab/ReViewCNV/assets/110105172/c31a6bb3-48df-4e0e-b212-3031c3cdb5ba)


<br/><br/>


#### Showing the genes annotations plot, filtering the identified regions and visualizing a different chromosome



![Video_2b](https://github.com/ctglab/ReViewCNV/assets/110105172/98089ba5-6874-4bdc-b3bb-e03b87772c01)


<br/><br/>

#### Filtering the population polymorphisms, changing public dataset and increasing the number of polymorphisms to visualize


![Video_3_new](https://github.com/ctglab/ReViewCNV/assets/110105172/f30854f3-f0a2-4097-904d-dd2af12d54a1)


<br/><br/>

#### Sharing the y axis range when loading multiple individuals and disabling the x axis sharing
![Video_4b](https://github.com/ctglab/ReViewCNV/assets/110105172/f1bac7ba-adcd-443c-8332-5189edf1b982)



<br/><br/>

#### Downloading an HTML file of a full chromosome or of a chromosomal region

![Video_5b](https://github.com/ctglab/ReViewCNV/assets/110105172/dbffa3da-a028-4970-a4e8-e377d57541f4)

