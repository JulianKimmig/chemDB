FROM frolvlad/alpine-miniconda3

RUN apk add default-libmysqlclient-dev

COPY ./environment.yml /environment.yml
COPY ./requirements.txt /requirements.txt

RUN conda env create --file environment.yml


RUN useradd --create-home appuser
ADD ./chemDB /home/appuser/app
WORKDIR /home/appuser/app
USER appuser

EXPOSE 8000

ENTRYPOINT ["conda", "run", "-n", "chemDB-env", "exec", "gunicorn", "chemDB.wsgi:application", "--bind","0.0.0.0:8000","--workers","3"]