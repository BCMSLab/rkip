FROM rocker/verse:latest

RUN apt-get update
RUN apt-get -y upgrade

RUN apt-get -y install libudunits2-dev
RUN apt-get -y install graphviz

RUN Rscript -e "update.packages(ask=FALSE)"

RUN Rscript -e "install.packages('tidygraph', quiet = TRUE)"
RUN Rscript -e "install.packages('ggraph', quiet = TRUE)"
RUN Rscript -e "install.packages('ggimage', quiet = TRUE)"
RUN Rscript -e "install.packages('cowplot', quiet = TRUE)"
RUN Rscript -e "install.packages('BiocManager', quiet = TRUE)"
RUN Rscript -e "BiocManager::install('AnnotationDbi', quiet = TRUE)"
RUN Rscript -e "BiocManager::install('org.Hs.eg.db', quiet = TRUE)"
RUN Rscript -e "BiocManager::install('clusterProfiler', quiet = TRUE)"
RUN Rscript -e "BiocManager::install('impute', quiet = TRUE)"
RUN Rscript -e "BiocManager::install('preprocessCore', quiet = TRUE)"
RUN Rscript -e "BiocManager::install('STRINGdb', quiet = TRUE)"
RUN Rscript -e "install.packages('cRegulome', quiet = TRUE)"
RUN Rscript -e "install.packages('cgdsr', quiet = TRUE)"
RUN Rscript -e "install.packages('WGCNA', quiet = TRUE)"
