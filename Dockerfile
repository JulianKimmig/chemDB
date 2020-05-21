FROM frolvlad/alpine-miniconda3

COPY ./environment.yml /environment.yml
COPY ./requirements.txt /requirements.txt

RUN conda env create --file environment.yml
RUN conda activate chemDB-env

ADD ./chemDB /usr/src/app

WORKDIR /usr/src/app
EXPOSE 8000

CMD exec gunicorn chemDB.wsgi:application --bind 0.0.0.0:8000 --workers 3