from django.urls import path, include

app_name="chemicaldb_polymers"

from . import views
urlpatterns = [
    path("",views.main,name="main"),
    path("structures",views.structures_main,name="structures_main"),
    path("new_polymer_structure",views.NewPolymerStructure.as_view(),name="new_polymer_structure"),
    path("structure_browse",views.structure_browse,name="structure_browse"),
    path("structure/<pk>",views.structure_view,name="structure"),



    path("new_polymer",views.NewPolymer.as_view(),name="new_polymer"),
    path("startendgroup",views.StartEndGroupEdit.as_view(),name="new_start_end_group"),
    path("startendgroup/<pk>",views.StartEndGroupEdit.as_view(),name="edit_start_end_group"),

    path("repeating_unit",views.RepeatingUnitEdit.as_view(),name="new_repeating_unit"),
    path("repeating_unit/<pk>",views.RepeatingUnitEdit.as_view(),name="edit_repeating_unit"),
]


