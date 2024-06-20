# ReViewCNV

ReViewCNV is a containerized Shiny App for the visualization the Copy Number Variants (CNVs) identified by the algorithm [EXCAVATOR2](https://pubmed.ncbi.nlm.nih.gov/27507884/). The App allows the user to compare the  identified CNVs with population polymorphisms present in different public datasets (AnnotSV, DGV and gnomAD) and to identify the genes present in correspondence of the altered regions (genome-wide).

<br/>


### Building and running the Docker

To build the image, download the repository, open the terminal, go to the folder with the downloaded files and run the command below (this can take some time):

sudo docker build -t shiny-app-excavator2 .

To run the container run the following command:

sudo docker run --name shiny-app-excavator2  -p 3838:3838 shiny-app-excavator2

If working on a Mac with Apple Silicon, run the command below instead:

sudo docker run --name shiny-app-excavator2 --platform linux/x86_64   -p 3838:3838 shiny-app-excavator2 

<br/>

### Building and runnung the Singularity/Apptainer

To build the file .sif is necessary to build the Docker immage first (see above), then the immage should be saved locally as a file .tar using the command below:

sudo docker save shiny-app-excavator2 > shiny-app-excavator2.tar 

The file .tar can be converted to file.sif using one of the two commands below:

singularity build shiny-app-excavator2.sif docker-archive://shiny-app-excavator2.tar (if using singularity)
apptainer build shiny-app-excavator2.sif docker-archive://shiny-app-excavator2.tar (if using apptainer)

To run the apptainer singularity use one of the two commands below:

singularity run    shiny-app-excavator2.sif (if using singularity)
apptainer run    shiny-app-excavator2.sif (if using apptainer)

<br/>

#### Specifying the genome assembly and interacting with the identified regions and the population polymorphisms plots

![Video_1](https://github.com/ctglab/ReViewCNV/assets/110105172/c0f3ad06-f099-4a71-b962-8d7dfb152513)


<br/><br/>


#### Showing the genes annotations plot, filtering the identified regions and to visualizing a different chromosome

![Video_2](https://github.com/ctglab/ReViewCNV/assets/110105172/6c770d9f-b347-4a99-8745-b61ece0def9f)


<br/><br/>

#### Filtering the population polymorphisms, changing the polymorphisms dataset and increasing the number of polymorphisms to visualize

![Video_3](https://github.com/ctglab/ReViewCNV/assets/110105172/27f238ae-884e-4fa8-8fd9-b90be701decd)


<br/><br/>

#### Sharing the y axis range when loading multiple individuals, and disabling the x axis sharing

![Video_4](https://github.com/ctglab/ReViewCNV/assets/110105172/8d4d34bc-4eee-4ab9-a8e7-95a5acaa93f7)


<br/><br/>

#### Downloading an HTML file of a full chromosome or of a chromosomal region

![Video_5](https://github.com/ctglab/ReViewCNV/assets/110105172/91823650-ccd4-4a69-a613-10043d012ff4)
