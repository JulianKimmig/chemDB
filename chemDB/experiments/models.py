from django.utils import timezone
from django.contrib.auth.models import User
from django.db import models


class ConcentrationUnits(models.TextChoices):
    MOL_LITER = 'M_L', "mol/L"
    GRAM_LITER = 'G_L', "g/L"

ConcentrationUnits.max_length = 3


class GeneraExperiment(models.Model):
    name = models.CharField(max_length=64)

    def __str__(self):
        return self.name

    class Meta:
        abstract = True


def upload_to_experiment(instance, filename):
    return 'experiments/{username}/{classname}/{date}_{filename}'.format(
        username=instance.user.user.username, classname=instance.__class__.__name__, date=timezone.now(),
        filename=filename)


class Experiment(models.Model):
    name = models.CharField(max_length=64)
    short_name = models.CharField(max_length=16)
    owner = models.ForeignKey(User, null=True, on_delete=models.CASCADE)
    run_date = models.DateTimeField(default=timezone.now, blank=True)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    raw_data = models.FileField(upload_to=upload_to_experiment)


class ExperimentDataXY(models.Model):
    experiment = models.ForeignKey(Experiment, on_delete=models.CASCADE)
    x = models.FloatField()
    y = models.FloatField()

    class Meta:
        abstract = True
        unique_together = ('experiment', 'x',)


from .sub_models import *
