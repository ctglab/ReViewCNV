# ReViewCNV

![GitHub release](https://img.shields.io/github/release/ctglab/ReViewCNV.svg)
![GitHub docker image build](https://github.com/ctglab/ReViewCNV/actions/workflows/docker-image.yml/badge.svg)
![GitHub commit activity](https://img.shields.io/github/commit-activity/m/ctglab/ReViewCNV) 
![GitHub last commit](https://img.shields.io/github/last-commit/ctglab/ReViewCNV)
![GitHub contributors](https://img.shields.io/github/contributors/ctglab/ReViewCNV)
[![GitHub Website](https://img.shields.io/website-up-down-green-red/http/monip.org.svg)](https://ctglab.github.io/)
![GitHub forks](https://img.shields.io/github/forks/ctglab/ReViewCNV)
![GitHub Repo stars](https://img.shields.io/github/stars/ctglab/ReViewCNV)
![GitHub watchers](https://img.shields.io/github/watchers/ctglab/ReViewCNV.svg)


ReViewCNV is a containerized Shiny App for the visualization of germinal Copy Number Variants (CNVs). It was created to visualize CNVs from exome or gene panels sequencing identified by the algorithm [EXCAVATOR2](https://pubmed.ncbi.nlm.nih.gov/27507884/)

<p align="justify">
ReViewCNV has now increased its compatibility and accepts as input a list of CNVs identified by any bioinformatic tool!
</p>
<p align="justify">
The App allows the user to compare the  CNVs of interst with population polymorphisms present in different public datasets (AnnotSV, DGV and gnomAD) and to identify the genes present in correspondence of the altered regions (genome-wide). For the genes present in correspondance of the CNVs of interest the exons are highlighted. It is possible to visualize up to three individuals at the same time in synchronized plots, facilitating family studies and the identification of de novo mutations.
</p>
<p align="justify">
In the container are present two apps: app_Excavator2.R and app_CNV_tsv.R. The app_Excavator2.R is specific for EXCAVATOR2. The app_CNV_tsv.R accepts as input a tsv file (without header) with three columns indicating for each CNV the genomic coordinates (Chromosome, Start and End) and optionally a fourth column specifying the CNV type (i.e. deletion or duplication). Example input files for the app_Excavator2.R are avialble in the folder TEST_TRIO_app_Excavator2, while for the app_CNV_tsv.R in the folder test_TRIO_app_CNV_tsv.
</p>


<br/>


### Building and running the Docker

<p align="justify">
To build the image, download the repository, open the terminal, go to the folder with the downloaded files and run the command below (this can take a long time):
</p>
_sudo docker buildx build  -t shiny-app-excavator2 ._

To run the app_Excavator2.R use:

_sudo docker run --name shiny-app-excavator2  -p 3838:3838 shiny-app-excavator2 Rscript app_Excavator2.R_

To run the app_CNV_tsv.R use:

_sudo docker run --name shiny-app-excavator2  -p 3838:3838 shiny-app-excavator2 Rscript app_CNV_tsv.R_

<br/>

### Building and running the Singularity/Apptainer

<p align="justify"> ##### To create the .sif file it is necessary to build the Docker image first (see above). The Docker image should be saved locally as a .tar file using the command below:
</p>
_sudo docker save shiny-app-excavator2 > shiny-app-excavator2.tar_ 

##### The .tar file  can be converted to a .sif file using one of the two commands below:

_singularity build shiny-app-excavator2.sif docker-archive://shiny-app-excavator2.tar_ (if using singularity)

_apptainer build shiny-app-excavator2.sif docker-archive://shiny-app-excavator2.tar_ (if using apptainer)

##### To run the the app_Excavator2.R on apptainer/singularity use one of the two commands below:

_singularity run  shiny-app-excavator2.sif  Rscript app_Excavator2.R_ (if using singularity)

_apptainer run  shiny-app-excavator2.sif  Rscript app_Excavator2.R_ (if using apptainer)

##### To run the the app_CNV_tsv.R on apptainer/singularity use one of the two commands below:

_singularity run  shiny-app-excavator2.sif  Rscript app_CNV_tsv.R_ (if using singularity)

_apptainer run  shiny-app-excavator2.sif  Rscript app_CNV_tsv.R_ (if using apptainer)

<br/>

### Example of use gifs:
<p align="justify">
The animated gifs below illustrate how the app works. The first one shows how to upload the input files in the app_CNV_tsv.R, while the following gifs show the app functionalities in the app_Excavator2.R as they are almost the same in the app_CNV_tsv.R.
</p>

<br/>

#### Uploading the input file in the app_CNV_tsv.R 

![App_CNV_tsv_upload](https://github.com/user-attachments/assets/24f2792b-4e30-4719-b60a-b6335c7b877e)


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

