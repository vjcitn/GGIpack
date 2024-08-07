FROM continuumio/miniconda3:latest

# replace this with the name of your shiny app's root folder
# ARG is for build-time, ENV is for runtime (https://stackoverflow.com/a/35562189)
ARG APPROOT=GGIpack
ENV SAPPROOT=$APPROOT

# copy the app to the image
RUN mkdir /root/$APPROOT
COPY $APPROOT /root/$APPROOT

# Install dependencies
WORKDIR /root/$APPROOT

# Some dependencies are available from anaconda
RUN conda env create -f conda.yaml

# Others need custom installation (make sure these installs happen *inside* the conda env we just built!)
# RUN conda run -n ggipack R -e "devtools::install_github('paul-shannon/igvShiny')"
# RUN conda run -n ggipack R -e "install.packages('dashboardthemes', repos='https://cloud.r-project.org/')"

EXPOSE 3838
ENTRYPOINT conda run -n ggipack --no-capture-output R -e "shiny::runApp('/root/$SAPPROOT/R', host='0.0.0.0', port=getOption('shiny.port',3838))"
