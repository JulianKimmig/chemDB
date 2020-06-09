from inspect import isclass

from crispy_forms.helper import FormHelper
from crispy_forms.layout import Submit
from django import forms
from django.db import models

# Create your models here.
from django.db.models.signals import pre_save
from django.dispatch import receiver
from polymorphic.models import PolymorphicModel

from .. import Structure, ChemdbUser

class Substance(PolymorphicModel):
    name = models.CharField(max_length=64)
    code = models.CharField(max_length=24, unique=True, null=True,)
    user = models.ForeignKey(ChemdbUser, on_delete=models.SET_NULL, null=True)
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


class SubstanceForm(forms.ModelForm):
    class Meta:
        model = Substance
        exclude = []

    def __init__(self,chem_db_user=None,changeable=False, *args, **kwargs):

        super().__init__(*args, **kwargs)
        self.helper = FormHelper()
        if changeable:
            self.helper.form_method = 'post'
            self.helper.add_input(Submit('submit', 'save'))
        if chem_db_user:
            self.fields["user"].initial = chem_db_user.pk
        self.fields["user"].widget = forms.HiddenInput()
#        initial={'owner':chem_db_user.pk,"short_name":}
#   self.fields['raw_data'] = forms.FileInput()

from .submodels import *

@receiver(pre_save)
def my_callback(sender, instance: Substance, *args, **kwargs):
    if not issubclass(sender, Substance):
        return
    instance.validate(save=False)

class SimpleSubstance(Substance):
    structure = models.ForeignKey(Structure, on_delete=models.CASCADE)

class SimpleSubstanceForm(SubstanceForm):
    class Meta:
        model = SimpleSubstance
        exclude = []

class MixedSubstance(Substance):
    structures = models.ManyToManyField(Structure)
