FROM pshannon/tmsft3
RUN apt-get install -y libmysqlclient21
RUN apt-get install -y libpq-dev
WORKDIR /home/rstudio
RUN mkdir /home/rstudio/localData
RUN mkdir /home/rstudio/localCode
RUN chown -R rstudio:rstudio /home/rstudio/localData
RUN chown -R rstudio:rstudio /home/rstudio/localCode
RUN chown -R rstudio:rstudio /home/rstudio/data
RUN chown -R rstudio:rstudio /home/rstudio/code
COPY components/installPackagesForTrena2.R /tmp
COPY components/trena_2.0.13.tar.gz /tmp
COPY --chown=rstudio:rstudio compose-v1/data/* /home/rstudio/data/
COPY --chown=rstudio:rstudio compose-v1/code/* /home/rstudio/code/
COPY klf1.Rmd /home/rstudio/code/
RUN R -f /tmp/installPackagesForTrena2.R
RUN R CMD INSTALL /tmp/trena_2.0.13.tar.gz



