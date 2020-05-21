FROM frolvlad/alpine-miniconda3

COPY ./environment.yml /environment.yml
RUN conda env create --file environment.yml