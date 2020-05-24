from django.contrib import admin

# Register your models here.
from experiments.models import ExperimentRawData
from .models import Nanoparticle, NanoparticleCharacterization, NanoparticlePreparationMethod

admin.site.register(Nanoparticle)

class NanoparticleCharacterizationAdminDataInline(admin.TabularInline):
    model = ExperimentRawData

class NanoparticleCharacterizationAdmin(admin.ModelAdmin):
    inlines = [
        NanoparticleCharacterizationAdminDataInline,
    ]

admin.site.register(NanoparticleCharacterization,NanoparticleCharacterizationAdmin)


admin.site.register(NanoparticlePreparationMethod)