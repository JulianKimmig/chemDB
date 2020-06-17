from django.db import models
from django.db.models.signals import pre_save
from django.dispatch import receiver


class ValidationError(Exception):
    pass


class ValidationModel(models.Model):
    class Meta:
        abstract = True

    valid = models.BooleanField(default=True, editable=False)

    def validate(self):
        valid = True

        for check_name, check in self.validity_checks.items():
            exception = None
            next_valid = check(self)
            valid = valid and next_valid
        return self.valid

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.validity_checks = {}


@receiver(pre_save)
def _validation_model_pre_save(sender, instance: ValidationModel, *args, **kwargs):
    if not issubclass(sender, ValidationModel):
        return
    try:
        instance.validate()
    except (ValidationError, ValueError):
        instance.valid = False



# Create your models here.
from .sub_models import *
