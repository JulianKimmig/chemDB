#!/bin/bash

# Collect static files
#echo "Collect static files"
#python manage.py collectstatic --noinput

# ony in development
  pip install -r /requirements.txt
  # Apply database migrations
  echo "Apply database migrations"
  python manage.py migrate

    echo "Collect static files"
   python manage.py collectstatic

# Start server
#echo "Starting server"
#python manage.py runserver 0.0.0.0:8000

echo "Starting server with gunicorn"
gunicorn chemDB.wsgi:application --bind 0.0.0.0:8000 --worker-tmp-dir /dev/shm --workers=2 --threads=4 --worker-class=gthread
tail -f /dev/null
