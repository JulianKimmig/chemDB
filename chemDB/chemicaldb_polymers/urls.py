from django.urls import path, include

app_name="chemicaldb_polymers"

from . import views
urlpatterns = [
    path("",views.main,name="main"),
]


