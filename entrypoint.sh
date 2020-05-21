#!/bin/bash
echo "HELLO"
conda run -n $ENVNAME exec gunicorn chemDB.wsgi:application --bind 0.0.0.0:8000 --workers 3
echo "BYE"