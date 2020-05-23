from . import views, api_views
from django.contrib.auth import views as auth_views
from django.urls import path, include
from django.contrib.auth.views import LoginView

app_name="chemicaldb"


urlpatterns = [

    path('', views.main_view, name ="main"),
    path('search', views.main_view, name ="search"),

    path('api/', include(api_views.router.urls)),
    path('api/search',api_views.SearchViewSet.as_view(),name="api_search"),
    path('api-auth/', include('rest_framework.urls')),

    path('login/', LoginView.as_view(template_name='chemicaldb/login.html'), name='login'),
    path('logout/', auth_views.auth_logout, name='logout'),

    path('experiments/', include('experiments.urls')),
]


