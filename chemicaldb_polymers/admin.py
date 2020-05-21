from django.contrib import admin

# Register your models here.
from chemicaldb_polymers.models import Polymer, PolymerStructure

for model in [Polymer,
              PolymerStructure,
              ]:
    admin.site.register(model)