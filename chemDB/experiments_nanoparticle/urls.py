from django.urls import path

from . import views

app_name="experiments_nanoparticle"
urlpatterns = [
    path('', views.main_view, name ="main"),
    path('batch_particle_creation',views.batch_particle_creation.as_view(),name="batch_particle_creation"),
    path('batch_characterizations',views.batch_characterizations,name="get_batch_characterizations"),
    path('characterization/<pk>',views.view_characterizations,name="view_characterization"),
    path('np/<pk>',views.view_particle,name="view_particle"),
]


