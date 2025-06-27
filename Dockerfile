#Base R Shiny image
FROM rhub/r-minimal

#Install bash
RUN apk update
RUN apk add --no-cache --update-cache bash

#Install linux libaries
RUN apk update
RUN apk add --no-cache --update-cache R R-dev openssl-dev pandoc \
gcc  tzdata g++ libffi-dev curl 

#Install R packages
RUN installr -d -t gfortran plotly
RUN installr -d \
        -t "zlib-dev" \
        shiny      
RUN installr -d dplyr 
RUN installr -d readr 
RUN installr -d htmlwidgets
RUN installr -d bslib
RUN installr -d \
    -t "make openssl-dev cmake linux-headers apache-arrow-dev" \
    -a "openssl libarrow_dataset libarrow" arrow@14.0.2.1

#Create work directory and copy required files and folders into it
RUN mkdir -p /app
WORKDIR /app
COPY inst/37 /app/37
COPY inst/38 /app/38
COPY inst/genes_annotation_37.rds  /app
COPY inst/genes_annotation_38.rds  /app
COPY inst/Exons_37.rds  /app
COPY inst/Exons_38.rds  /app
COPY inst/www /app/www
COPY inst/hg38_Coordinates.rds /app
COPY inst/hg37_Coordinates.rds /app
COPY tools/app_all_chromosomes.R /app




