from django.contrib.auth import views as auth_views
from django.contrib.auth.views import LoginView
from django.urls import path, include

from . import views, api_urls

app_name="chemicaldb"


urlpatterns = [

    path('', views.main_view, name ="main"),
    path('search', views.main_view, name ="search"), # TODO: implement search page

    path('substances', views.substance_main, name ="substance_main"),

    path('structures', views.structures_main, name ="structures_main"),
    path('structures/browse', views.structures_browse, name ="structures_browse"),

    path('substance/<pk>', views.substance_view, name ="substance_view"),

    path('api/', include(api_urls)),


    path('login/', LoginView.as_view(template_name='chemicaldb/login.html'), name='login'),
    path('logout/', auth_views.auth_logout, name='logout'),

    path('experiments/', include('experiments.urls')),
    path('polymers/', include('chemicaldb_polymers.urls')),
]


