from django.urls import path, include

from . import api_views

app_name="chemicaldb_polymers"
urlpatterns = [
    path('structure_checker', api_views.structure_checker,name="structure_checker"),
]