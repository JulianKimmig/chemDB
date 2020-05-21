FROM frolvlad/alpine-miniconda3

COPY ./environment.yml /environment.yml
COPY ./requirements.txt /requirements.txt
COPY ./chemDB /app

ADD ./chemDB /usr/src/app
WORKDIR /usr/src/app

# Expose ports
EXPOSE 8000

RUN conda env create --file environment.yml

CMD exec gunicorn chemDB.wsgi:application --bind 0.0.0.0:8001 --workers 3