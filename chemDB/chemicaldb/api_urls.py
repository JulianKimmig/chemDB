from django.urls import path, include

from . import api_views

urlpatterns = [
path('', include(api_views.router.urls)),
path('search',api_views.search,name="api_search"),
path("api_smiles_to_image",api_views.smiles_to_image,name="api_smiles_to_image"),
path('auth/', include('rest_framework.urls')),
path('experiments/', include('experiments.api_urls')),
]