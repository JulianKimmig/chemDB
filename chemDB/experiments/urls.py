from django.urls import path, include

app_name="experiments"

from . import views
urlpatterns = [
    path("",views.main,name="main"),
    path('nanoparticle/',include('experiments_nanoparticle.urls')),
]


