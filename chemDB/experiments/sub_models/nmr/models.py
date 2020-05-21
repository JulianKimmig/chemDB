from django.db import models

from chemDB.experiments import Experiment, ExperimentDataXY


class Nmr(Experiment):
    class Meta:
        abstract = True


class D1NmrChoices(models.TextChoices):
    H1 = "H1", "1HNmr"
    C13 = "C13", "13CNMR"

class D1Nmr(Nmr):
    type = models.CharField(max_length=16, choices=D1NmrChoices.choices, default=D1NmrChoices.H1)


class D1NmrData(ExperimentDataXY):
    experiment = models.ForeignKey(D1Nmr, on_delete=models.CASCADE)
    x = models.FloatField()
    y = models.FloatField()

