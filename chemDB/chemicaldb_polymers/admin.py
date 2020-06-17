from django.contrib import admin

# Register your models here.
from chemicaldb.admin import StructureAdmin
from chemicaldb.sub_models import StructureName
from chemicaldb_polymers.models import Polymer, PolymerStructure, StartEndGroup, RepeatingUnit


class StartEndGroupAdmin(StructureAdmin):
    pass
admin.site.register(StartEndGroup,StartEndGroupAdmin)

class RepeatingUnitAdmin(StructureAdmin):
    pass
admin.site.register(RepeatingUnit,RepeatingUnitAdmin)



class PolymerStructureAdmin(StructureAdmin):
    pass

admin.site.register(PolymerStructure,PolymerStructureAdmin)

for model in [Polymer,
              ]:
    admin.site.register(model)