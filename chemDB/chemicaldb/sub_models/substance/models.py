from django.contrib.auth.models import User
from django.db import models

# Create your models here.
from .. import Structure


class Substance(models.Model):
    name = models.CharField(max_length=64)
    code = models.CharField(max_length=16, unique=True, null=True,)
    user = models.ForeignKey(User, on_delete=models.SET_NULL, null=True)

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
