from django.db import models

# Create your models here.
from structure.models import Structure


class Substance(models.Model):
    pass



from .submodels import *


class SimpleSubstance(Substance):
    structure=models.ForeignKey(Structure,on_delete=models.CASCADE)

class MixedSubstance(Substance):
    structures = models.ManyToManyField(Structure)

