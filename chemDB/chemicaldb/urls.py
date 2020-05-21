from django.urls import path, include

from . import views, api_views

app_name="chemicaldb"


urlpatterns = [
    path('', views.main_view, name ="index"),
    path('search', views.main_view, name ="search"),

    path('api/', include(api_views.router.urls)),
    path('api-auth/', include('rest_framework.urls'))
]


