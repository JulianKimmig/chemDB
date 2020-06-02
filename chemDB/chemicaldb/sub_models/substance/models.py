from inspect import isclass

from django.db import models

# Create your models here.
from django.db.models.signals import pre_save
from django.dispatch import receiver

from .. import Structure, ChemdbUser


class SubstanceClass(models.Model):
    registered_substance_classes_by_name={}
    registered_substance_classes_by_class={}

    class_name = models.CharField(max_length=100,unique=True)
    def __str__(self):
        return self.class_name

    @staticmethod
    def register_substance_class(substance_class):
        assert issubclass(substance_class,Substance)
        name=substance_class.__module__+"."+substance_class.__name__
        instance,new=SubstanceClass.objects.get_or_create(class_name=name)
        SubstanceClass.registered_substance_classes_by_name[name]=substance_class,instance
        SubstanceClass.registered_substance_classes_by_class[substance_class]=name,instance
        return instance

    @staticmethod
    def get_substance_class(object):
        if not isclass(object):
            object = object.__class__
        if not issubclass(object,Substance):
            return None
        name,instance = SubstanceClass.registered_substance_classes_by_class[object]
        return instance



class Substance(models.Model):
    name = models.CharField(max_length=64)
    code = models.CharField(max_length=24, unique=True, null=True,)
    user = models.ForeignKey(ChemdbUser, on_delete=models.SET_NULL, null=True)
    valid = models.BooleanField(default=False,editable=False)
    substance_class = models.ForeignKey(SubstanceClass,on_delete=models.SET_NULL,null=True,blank=True)


    def __str__(self):
        if self.code:
            return "{} ({})".format(self.code,self.name)
        else:
            return self.name

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.validity_checks={}
        if self.substance_class is None:
            self.substance_class = self.get_substance_class()
            self.save()

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

    @classmethod
    def register_substance_class(cls):
        instance = SubstanceClass.register_substance_class(cls)

    def get_substance_class(self):
        return SubstanceClass.get_substance_class(self)

from .submodels import *

@receiver(pre_save)
def my_callback(sender, instance: Substance, *args, **kwargs):
    if not issubclass(sender, Substance):
        return
    instance.validate(save=False)

class SimpleSubstance(Substance):
    structure = models.ForeignKey(Structure, on_delete=models.CASCADE)
SimpleSubstance.register_substance_class()

class MixedSubstance(Substance):
    structures = models.ManyToManyField(Structure)
MixedSubstance.register_substance_class()