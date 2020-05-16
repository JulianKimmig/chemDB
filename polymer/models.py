from django.db import models

# Create your models here.
from structure.models import Structure
from substance.models import SimpleSubstance


class PolymerStructure(models.TextChoices):
    LINEAR = 'LINEAR', "linear"
    BRANCHED = 'BRANCHED', "branched"
    NETWORK = 'NETWORK', "network"


class Polymer(SimpleSubstance):
    polymer_structure = models.CharField(max_length=16, choices=PolymerStructure.choices,
                                         default=PolymerStructure.LINEAR)

    class Meta:
        abstract = True


class HomoPolymer(Polymer):
    monomer = models.ForeignKey(Structure, on_delete=models.CASCADE)
    repeating_units = models.PositiveIntegerField()


class MonomerDistribution(models.TextChoices):
    GRADIENT = 'GRADIENT', "gradient"
    STATISTICAL = 'STATISTICAL', "statistical"
    BLOCK = 'BLOCK', "block"


class Copolymer(Polymer):
    monomers = models.ManyToManyField(Structure, through='Copolymer_monomers')
    monomer_distribution = models.CharField(max_length=16, choices=MonomerDistribution.choices,
                                            default=MonomerDistribution.STATISTICAL)
    distribution_factor = models.FloatField(default=1)


class Copolymer_monomers(models.Model):
    copolymer = models.ForeignKey(Copolymer, on_delete=models.CASCADE)
    monomer = models.ForeignKey(Structure, on_delete=models.CASCADE)
    repeating_units = models.PositiveIntegerField()


from .submodels import *
