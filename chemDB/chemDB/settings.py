"""
Django settings for chemDB project.

Generated by 'django-admin startproject' using Django 2.2.

For more information on this file, see
https://docs.djangoproject.com/en/2.2/topics/settings/

For the full list of settings and their values, see
https://docs.djangoproject.com/en/2.2/ref/settings/
"""

import os

import datetime
from django.utils.timezone import make_aware
from json_dict import JsonDict
import random

# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/2.2/howto/deployment/checklist/

# SECURITY WARNING: keep the secret key used in production secret!

if os.path.exists("/django_settings.json"):
    config = JsonDict("/django_settings.json")
else:
    config = JsonDict("django_settings.json")


SECRET_KEY = config.get('base_settings','SECRET_KEY',default="".join([random.choice("abcdefghijklmnopqrstuvwxyz0123456789!@#$%^&*(-_=+)") for i in range(50)]))
#'l5!r&kcwz5r+#%**eq8j*6p-6@nw&0h)1uf&9c^27sw)9_s^-='

# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = config.get('base_settings','DEBUG',default=False)

ALLOWED_HOSTS = config.get('base_settings','ALLOWED_HOSTS',default=[])
# Application definition

INSTALLED_APPS = [
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',

    #'haystack',
    'polymorphic',
    'rest_framework',
    'guardian',
    'crispy_forms',

    'chemicaldb.apps.ChemicaldbConfig',

    'experiments.apps.ExperimentsConfig',
    'sources.apps.SourcesConfig',

    'chemicaldb_polymers.apps.ChemicaldbPolymersConfig',

    'experiments_nanoparticle.apps.ExperimentsNanoparticleConfig',
]

MIDDLEWARE = [
    'django.middleware.security.SecurityMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
]

ROOT_URLCONF = 'chemDB.urls'

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                'django.template.context_processors.debug',
                'django.template.context_processors.request',
                'django.contrib.auth.context_processors.auth',
                'django.contrib.messages.context_processors.messages',
            ],
        },
    },
]

AUTHENTICATION_BACKENDS = (
    'django.contrib.auth.backends.ModelBackend',  # default
    'guardian.backends.ObjectPermissionBackend',
)

WSGI_APPLICATION = 'chemDB.wsgi.application'

# Database
# https://docs.djangoproject.com/en/2.2/ref/settings/#databases

try:
    DATABASES = config.get('base_settings','DATABASES',default = {
        'default': {
            'ENGINE': 'django.db.backends.mysql',
            'NAME': os.environ['MYSQL_DATABASE'],
            'USER': os.environ['MYSQL_USER'],
            'PASSWORD': os.environ['MYSQL_PASSWORD'],
            'HOST': os.environ.get('MYSQL_HOST', default="127.0.0.1"),
            'PORT': os.environ.get('MYSQL_PORT', default="3306"),
        }
    })
except:
    DATABASES = config.get('base_settings','DATABASES',default = {
        'default': {
            'ENGINE': 'django.db.backends.sqlite3',
            'NAME': os.path.join(BASE_DIR, 'db.sqlite3'),
        }
    })


# Password validation
# https://docs.djangoproject.com/en/2.2/ref/settings/#auth-password-validators

AUTH_PASSWORD_VALIDATORS = [
    {
        'NAME': 'django.contrib.auth.password_validation.UserAttributeSimilarityValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.MinimumLengthValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.CommonPasswordValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.NumericPasswordValidator',
    },
]

# Internationalization
# https://docs.djangoproject.com/en/2.2/topics/i18n/

LANGUAGE_CODE = config.get('base_settings','LANGUAGE_CODE',default='en-us')

TIME_ZONE = config.get('base_settings','TIME_ZONE',default='UTC')

USE_I18N = config.get('base_settings','USE_I18N',default=True)

USE_L10N = config.get('base_settings','USE_L10N',default=True)

USE_TZ = config.get('base_settings','USE_TZ',default=True)

# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/2.2/howto/static-files/

STATIC_URL = config.get('base_settings','STATIC_URL',default='/chemdb_static/')

FORCE_SCRIPT_NAME = config.get('base_settings','FORCE_SCRIPT_NAME',default=None)

STATIC_ROOT = config.get('base_settings','STATIC_ROOT',default='/static')

REST_FRAMEWORK = {
    # Use Django's standard `django.contrib.auth` permissions,
    # or allow read-only access for unauthenticated users.
    'DEFAULT_PERMISSION_CLASSES': [
        'rest_framework.permissions.DjangoModelPermissionsOrAnonReadOnly'
    ]
}

CRISPY_TEMPLATE_PACK = 'bootstrap4'

HAYSTACK_CONNECTIONS = {
    'default': {
        'ENGINE': 'haystack.backends.whoosh_backend.WhooshEngine',
        'PATH': os.path.join(os.path.dirname(__file__), 'whoosh_index'),
    },
}

# import polydeep_base_server.common.settings
# locals().update(
#    polydeep_base_server.common.settings.apply(
#        project_name = "chemDB",
#    )
# )
aware_datetime = make_aware(datetime.datetime.now())

LOGIN_REDIRECT_URL ="/"