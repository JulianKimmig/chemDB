from django.urls import path, include

#from . import api_views

urlpatterns = [
#path('', include(api_views.router.urls)),
#path('search',api_views.SearchViewSet.as_view(),name="api_search"),
#path('auth/', include('rest_framework.urls')),
path('nanoparticle/', include('experiments_nanoparticle.api_urls')),
]