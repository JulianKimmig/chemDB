from django.urls import path

from . import views

app_name="experiments_nanoparticle"
urlpatterns = [
    path('', views.main_view, name ="main"),
    path('batch_particle_creation',views.batch_particle_creation.as_view(),name="batch_particle_creation"),
    path('batch_characterizations',views.batch_characterizations,name="get_batch_characterizations"),
    path('characterization/<pk>',views.view_characterizations,name="view_characterization"),

    path('add_preparation_method',views.AddPreparationMethod.as_view(),name="add_preparation_method"),

    path('batch_edit_characterization',views.batch_edit_characterization,name="batch_edit_characterization"),
    path('batch_edit_particles',views.BatchEditParticles.as_view(),name="batch_edit_particles"),



    path('np/<pk>',views.view_particle,name="view_particle"),
]


