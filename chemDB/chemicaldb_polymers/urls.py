from django.urls import path, include

app_name="chemicaldb_polymers"

from . import views
urlpatterns = [
    path("",views.main,name="main"),
    path("new_polymer",views.NewPolymer.as_view(),name="new_polymer"),
    path("startendgroup",views.StartEndGroupEdit.as_view(),name="new_start_end_group"),
    path("startendgroup/<pk>",views.StartEndGroupEdit.as_view(),name="edit_start_end_group"),

    path("repeating_unit",views.RepeatingUnitEdit.as_view(),name="new_repeating_unit"),
    path("repeating_unit/<pk>",views.RepeatingUnitEdit.as_view(),name="edit_repeating_unit"),
]


