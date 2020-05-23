from django.urls import path, include

app_name="experiments"
urlpatterns = [
    path('nanoparticle',include('experiments_nanoparticle.urls')),
]


