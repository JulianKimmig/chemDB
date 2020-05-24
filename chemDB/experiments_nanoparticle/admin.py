from django.contrib import admin

# Register your models here.
from experiments.models import ExperimentRawData, ExperimentData, ExperimentSingleValueData
from .models import Nanoparticle, NanoparticleCharacterization, NanoparticlePreparationMethod, Materials, Additives


class NanoparticleMaterialsInline(admin.TabularInline):
    model = Materials
    fk_name = 'nanoparticle'

class NanoparticleAdditivesInline(admin.TabularInline):
    model = Additives
    fk_name = 'nanoparticle'

class NanoparticleAdmin(admin.ModelAdmin):
    inlines = [
        NanoparticleMaterialsInline,
        NanoparticleAdditivesInline,
    ]

admin.site.register(Nanoparticle,NanoparticleAdmin)

class NanoparticleCharacterizationAdminFileInline(admin.TabularInline):
    model = ExperimentRawData

class NanoparticleCharacterizationAdminDataInline(admin.TabularInline):
    model = ExperimentData

class NanoparticleCharacterizationAdmintSingleValueDataInline(admin.TabularInline):
    model = ExperimentSingleValueData

class NanoparticleCharacterizationAdmin(admin.ModelAdmin):
    inlines = [
        NanoparticleCharacterizationAdminFileInline,
        NanoparticleCharacterizationAdmintSingleValueDataInline,
        NanoparticleCharacterizationAdminDataInline,
    ]

admin.site.register(NanoparticleCharacterization,NanoparticleCharacterizationAdmin)


admin.site.register(NanoparticlePreparationMethod)