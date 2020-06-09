from django.db import models
from django.db.models.signals import pre_save
from django.dispatch import receiver


class ValidationModel(models.Model):
    class Meta:
        abstract = True

    valid = models.BooleanField(default=True, editable=False)

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

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.validity_checks = {}

@receiver(pre_save)
def my_callback(sender, instance: ValidationModel, *args, **kwargs):
    if not issubclass(sender, ValidationModel):
        return
    instance.validate(save=False)
# Create your models here.
from .sub_models import *


