from django.urls import path, include

app_name="chemicaldb_polymers"

from . import views
urlpatterns = [
    path("",views.main,name="main"),
    path("new_polymer",views.NewPolymer.as_view(),name="new_polymer"),
    path("startendgroup",views.StartEndGroup.as_view(),name="new_start_end_group"),
    path("startendgroup/<pk>",views.StartEndGroup.as_view(),name="edit_start_end_group"),
]


