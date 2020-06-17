from django.contrib import admin

# Register your models here.
import chemicaldb.models as cdb_models

for model in [cdb_models.SimpleSubstance,
              cdb_models.MixedSubstance,
                cdb_models.Substance,
                cdb_models.StructureName
              ]:
    admin.site.register(model)

class StructureAdminAdminStructureNameInline(admin.TabularInline):
    model = cdb_models.StructureName

class StructureAdmin(admin.ModelAdmin):
    inlines = [
        StructureAdminAdminStructureNameInline,
    ]
    readonly_fields = ["structure_image_field","valid","iso_smiles"]

admin.site.register(cdb_models.Structure,StructureAdmin)


for model in []:
    admin.site.register(model)

for model in [cdb_models.ChemdbUser,
              cdb_models.ChemdbInstitute,
              ]:
    admin.site.register(model)