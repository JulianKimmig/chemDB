from django.db import models

# Create your models here.
from chemicaldb.sub_models import Substance, Structure
from experiments.models import GeneraExperiment, ConcentrationUnits, Experiment


class NanoparticleCharacterizationTool(models.TextChoices):
    MALVEN_ZETASIZER_NANO="MZetaNano","Malvern Zetasizer Nano"

NanoparticleCharacterizationTool.max_length = 9

class NanoparticleCharacterization(Experiment):
    tool= models.CharField(max_length=NanoparticleCharacterizationTool.max_length, choices=NanoparticleCharacterizationTool.choices,
                           default=NanoparticleCharacterizationTool.MALVEN_ZETASIZER_NANO)


class NanoparticlePreparationMethod(GeneraExperiment):
    pass


def _validate_prep_method(np):
    if np.preparation_method is None:
        return False
    return True

def _validate_material(np):
    if len(np.materials) <=0:
        return False
    return True

NP_VALIDATION_PREPARATION="np_preparation_method",_validate_prep_method
NP_VALIDATION_MATERIAL="np_material",_validate_material

class Nanoparticle(Substance):
    preparation_method = models.ForeignKey(NanoparticlePreparationMethod, on_delete=models.SET_NULL, null=True)
    materials = models.ManyToManyField(Substance, through='Materials',related_name="as_nanopartice_material")
    additives = models.ManyToManyField(Substance, through='Additives',related_name="as_nanopartice_additive")
    characterizations = models.ManyToManyField(NanoparticleCharacterization)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.validity_checks.update({
            NP_VALIDATION_PREPARATION[0]:NP_VALIDATION_PREPARATION[1]
        })

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
