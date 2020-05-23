from django.urls import path, include
from . import views, api_views
app_name="experiments_nanoparticle"
urlpatterns = [
    path('', views.main_view, name ="main"),
]


