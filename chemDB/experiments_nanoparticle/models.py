from django.db import models

# Create your models here.
from chemicaldb.sub_models import Substance, Structure
from experiments.models import GeneraExperiment, ConcentrationUnits


class NanoparticlePreparationMethod(GeneraExperiment):
    pass


class Nanoparticle(Substance):
    preparation_method = models.ForeignKey(NanoparticlePreparationMethod, on_delete=models.SET_NULL, null=True)
    materials = models.ManyToManyField(Substance, through='Materials',related_name="as_nanopartice_material")
    additives = models.ManyToManyField(Substance, through='Additives',related_name="as_nanopartice_additive")


class Materials(models.Model):
    nanoparticle = models.ForeignKey(Nanoparticle, on_delete=models.CASCADE, related_name="np_material_particle")
    material = models.ForeignKey(Substance, on_delete=models.CASCADE,related_name="np_material" )
    relative_content = models.FloatField()

class Additives(models.Model):
    nanoparticle = models.ForeignKey(Nanoparticle, on_delete=models.CASCADE, related_name="np_additives_particle")
    material = models.ForeignKey(Substance, on_delete=models.CASCADE,related_name="np_additive" )
    concentration = models.FloatField()
    concentration_unit = models.CharField(max_length=ConcentrationUnits.max_length, choices=ConcentrationUnits.choices,
                                          default=ConcentrationUnits.MOL_LITER)

