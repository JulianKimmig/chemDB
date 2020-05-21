from django.contrib import admin

# Register your models here.
from chemDB.chemicaldb_polymers import Polymer, PolymerStructure

for model in [Polymer,
              PolymerStructure,
              ]:
    admin.site.register(model)