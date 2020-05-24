from django.urls import path, include
from . import views, api_views
app_name="experiments_nanoparticle"
urlpatterns = [
    path('', views.main_view, name ="main"),
    path('batch_particle_creation',views.batch_particle_creation.as_view(),name="batch_particle_creation")
]


