# ReViewCNV

![example workflow](https://github.com/ctglab/ReViewCNV/actions/workflows/docker-image.yml/badge.svg?branch=master)

ReViewCNV is a containerized Shiny App for the visualization of the Copy Number Variants (CNVs) from exome or gene panels sequencing identified by the algorithm  [EXCAVATOR2](https://pubmed.ncbi.nlm.nih.gov/27507884/). The App allows the user to compare the  identified CNVs with population polymorphisms present in different public datasets (AnnotSV, DGV and gnomAD) and to identify the genes present in correspondence of the altered regions (genome-wide). It is possible to visualize the identified regions in up to three individuals (Trio) at the same time. A test dataset is provided in the TEST_TRIO folder.


<br/>


### Building and running the Docker

To build the image, download the repository, open the terminal, go to the folder with the downloaded files and run the command below (this can take a long time):

_sudo docker buildx build  -t shiny-app-excavator2 ._

To run the container run:

_sudo docker run --name shiny-app-excavator2  -p 3838:3838 shiny-app-excavator2_

<br/>

### Building and running the Singularity/Apptainer

To create the .sif file it is necessary to build the Docker image first (see above). The Docker image should be saved locally as a .tar file using the command below:

_sudo docker save shiny-app-excavator2 > shiny-app-excavator2.tar_ 

The .tar file  can be converted to a .sif file using one of the two commands below:

_singularity build shiny-app-excavator2.sif docker-archive://shiny-app-excavator2.tar_ (if using singularity)

_apptainer build shiny-app-excavator2.sif docker-archive://shiny-app-excavator2.tar_ (if using apptainer)

To run the apptainer/singularity use one of the two commands below:

_singularity run    shiny-app-excavator2.sif_ (if using singularity)

_apptainer run    shiny-app-excavator2.sif_ (if using apptainer)

<br/>

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

