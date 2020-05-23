from django.contrib.auth.models import User
from django.db import models

# Create your models here.
from django.db.models.signals import pre_save
from django.dispatch import receiver

from .. import Structure


class Substance(models.Model):
    name = models.CharField(max_length=64)
    code = models.CharField(max_length=16, unique=True, null=True,)
    user = models.ForeignKey(User, on_delete=models.SET_NULL, null=True)
    valid = models.BooleanField(default=False,editable=False)

    def __str__(self):
        if self.code:
            return "{} ({})".format(self.code,self.name)
        else:
            return self.name

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.validity_checks={}

    def validate(self, save=True):
        valid = True
        for check_name, check in self.validity_checks.items():
            valid = valid and check(self)
            if not valid:
                break
        if self.valid != valid:
            self.valid = valid
            if save:
                self.save()
        return self.valid

from .submodels import *

@receiver(pre_save)
def my_callback(sender, instance: Substance, *args, **kwargs):
    if not issubclass(sender, Substance):
        return
    instance.validate(save=False)

class SimpleSubstance(Substance):
    structure = models.ForeignKey(Structure, on_delete=models.CASCADE)


class MixedSubstance(Substance):
    structures = models.ManyToManyField(Structure)
