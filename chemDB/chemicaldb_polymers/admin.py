from django.contrib import admin

# Register your models here.
from chemicaldb.sub_models import StructureName
from chemicaldb_polymers.models import Polymer, PolymerStructure

class PolymerStructureAdminStructureNameInline(admin.TabularInline):
    model = StructureName

class PolymerStructureAdmin(admin.ModelAdmin):
    inlines = [
        PolymerStructureAdminStructureNameInline,
    ]

admin.site.register(PolymerStructure,PolymerStructureAdmin)

for model in [Polymer,
              ]:
    admin.site.register(model)