from django.db import models

from experiments.models import Experiment, ExperimentDataXY


class SECTypeAdditives(models.Model):
    SEC_type = models.ForeignKey('SECType', on_delete=models.CASCADE)
    structure = models.ForeignKey('chemicaldb.Structure', on_delete=models.CASCADE,related_name="SEC_additive")
    concentration = models.FloatField()


class SECType(models.Model):
    solvent = models.ForeignKey('chemicaldb.Structure', on_delete=models.SET_NULL, null=True,related_name="SEC_solvent")
    additives = models.ManyToManyField('chemicaldb.Structure', through=SECTypeAdditives)


class SEC(Experiment):
    type = SECType


class SECRawData(ExperimentDataXY):
    experiment = models.ForeignKey(SEC, on_delete=models.CASCADE)
    x = models.FloatField()
    y = models.FloatField()

class SECCalibration(SEC):
    pass
