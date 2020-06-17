from django.db import models

# Create your models here.
# Create your models here.
from django.db.models.signals import pre_save
from django.dispatch import receiver
from polymorphic.models import PolymorphicModel

from .. import Structure, ChemdbUser, ChemDbShareModel, ValidationModel


class Substance(PolymorphicModel,ValidationModel,ChemDbShareModel):
    name = models.CharField(max_length=64)

    def __str__(self):
        if self.code:
            return "{} ({})".format(self.code,self.name)
        else:
            return self.name

from .submodels import *

class SimpleSubstance(Substance):
    structure = models.ForeignKey(Structure, on_delete=models.CASCADE)

class MixedSubstance(Substance):
    structures = models.ManyToManyField(Structure)
