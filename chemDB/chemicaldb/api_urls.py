from django.urls import path, include

from . import api_views
app_name="api"

urlpatterns = [
    path('', include(api_views.router.urls)),
    path('search',api_views.search,name="api_search"),
    path("smiles_checker",api_views.smiles_checker,name="smiles_checker"),

    path('polymers/', include('chemicaldb_polymers.api_urls')),

path("api_smiles_to_image",api_views.smiles_to_image,name="api_smiles_to_image"),

path('auth/', include('rest_framework.urls')),
path('experiments/', include('experiments.api_urls')),
]