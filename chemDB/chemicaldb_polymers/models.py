from django.db import models

# Create your models here.
from chemDB.chemicaldb.models import Structure, Substance


class PolymerStructureType(models.TextChoices):
    LINEAR = 'LINEAR', "linear"
    BRANCHED = 'BRANCHED', "branched"
    NETWORK = 'NETWORK', "network"


class PolymerStructure(Structure):
    big_smiles = models.CharField(max_length=200, null=True, blank=True)

class MonomerDistribution(models.TextChoices):
    GRADIENT = 'GRADIENT', "gradient"
    STATISTICAL = 'STATISTICAL', "statistical"
    BLOCK = 'BLOCK', "block"

class Polymer(Substance):
    polymer_structure = models.ForeignKey(PolymerStructure, on_delete=models.CASCADE,related_name="polymer_structure")
    polymer_structure_type = models.CharField(max_length=16, choices=PolymerStructureType.choices,
                                              default=PolymerStructureType.LINEAR)
    monomers = models.ManyToManyField(Structure, through='Monomers')
    monomer_distribution = models.CharField(max_length=16, choices=MonomerDistribution.choices,
                                            default=MonomerDistribution.STATISTICAL)
    distribution_parameter = models.FloatField(default=1)

#    def __str__(self):


class Monomers(models.Model):
    polymer = models.ForeignKey(Polymer, on_delete=models.CASCADE,related_name="p")
    monomer = models.ForeignKey(Structure, on_delete=models.CASCADE,)
    repeating_units = models.PositiveIntegerField()
