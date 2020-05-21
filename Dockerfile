FROM frolvlad/alpine-miniconda3

COPY ./environment.yml /environment.yml
COPY ./requirements.txt /requirements.txt

RUN conda env create --file environment.yml

ADD ./chemDB /usr/src/app

WORKDIR /usr/src/app
EXPOSE 8000

ENTRYPOINT ["conda", "run", "-n", "chemDB-env", "exec", "gunicorn", "chemDB.wsgi:application", "--bind","0.0.0.0:8000","--workers","3"]