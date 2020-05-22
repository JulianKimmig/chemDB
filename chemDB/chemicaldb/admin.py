from django.contrib import admin

# Register your models here.
import chemicaldb.models as cdb_models

for model in [cdb_models.SimpleSubstance,
              cdb_models.MixedSubstance,
              ]:
    admin.site.register(model)



for model in [cdb_models.Structure,
              cdb_models.StructureName,
              ]:
    admin.site.register(model)
