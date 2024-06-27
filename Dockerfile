#Base r-minimal image
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
RUN installr -d shinyHugePlot
RUN installr -d bslib


#Create work directory and copy required files and folders into it
RUN mkdir -p /app
WORKDIR /app
COPY Annotations_37.rds /app
COPY Annotations_38.rds /app
COPY genes_annotation_37.rds  /app
COPY genes_annotation_38.rds  /app
COPY www /app/www
COPY app.R /app

#Run app.R
CMD ["Rscript", "/app/app.R"]








