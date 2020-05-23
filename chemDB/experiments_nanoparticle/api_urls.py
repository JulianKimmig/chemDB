from django.urls import path, include

from . import api_views

urlpatterns = [
path('', include(api_views.router.urls)),
]