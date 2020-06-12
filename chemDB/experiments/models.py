from django.db import models
from django.utils import timezone

from chemicaldb.sub_models import ChemdbUser, ChemDbShareModel


class ConcentrationUnits(models.TextChoices):
    MOL_LITER = 'M_L', "mol/L"
    GRAM_LITER = 'G_L', "g/L"


ConcentrationUnits.max_length = 3


class GeneraExperiment(ChemDbShareModel):
    name = models.CharField(max_length=64)
    procedure = models.TextField(default="", blank=True)

    def __str__(self):
        return self.name

    class Meta:
        abstract = True


def upload_to_experiment(instance, filename):
    return 'data/experiments/{classname}/{date}_{username}_{filename}'.format(
        username=instance.experiment.owner.auth_user.username, classname=instance.experiment.__class__.__name__,
        date=timezone.now(),
        filename=filename)


class Experiment(models.Model):
    name = models.CharField(max_length=64)
    short_name = models.CharField(max_length=16)
    owner = models.ForeignKey(ChemdbUser, null=True, on_delete=models.SET_NULL)
    run_date = models.DateTimeField(default=timezone.now, blank=True)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    batch_experiment = models.ForeignKey("self", default=None, on_delete=models.CASCADE, null=True,
                                         related_name="sub_experiments")
    batch_experiment_index = models.PositiveIntegerField(null=True, default=None)

    def get_data(self):
        return self.data

    class Meta:
        unique_together = ('short_name', 'owner',)

    def __str__(self):
        return "{} ({})".format(self.short_name, self.owner)


class ExperimentRawData(models.Model):
    name = models.CharField(max_length=32)
    raw_data = models.FileField(upload_to=upload_to_experiment)
    experiment = models.ForeignKey(Experiment, related_name="data", on_delete=models.CASCADE)

    class Meta:
        unique_together = ('name', 'experiment')

    def __str__(self):
        return "{} ({})".format(self.name, self.experiment)


class DataTypesChoices(models.TextChoices):
    STRING = ("str", "string")
    FLOAT = ("float", "decimal")
    INT = ("int", "integer")
    BOOL = ("bool", "boolean")


class ExperimentData(models.Model):
    experiment = models.ForeignKey(Experiment, on_delete=models.CASCADE)
    name = models.CharField(max_length=32)
    type = models.CharField(max_length=5, choices=DataTypesChoices.choices, default=DataTypesChoices.FLOAT)

    class Meta:
        unique_together = ('experiment', 'name',)


class ExperimentSingleValueData(ExperimentData):
    value = models.CharField(max_length=32)


class ExperimentDataPointXY(models.Model):
    experiment_data = models.ForeignKey(ExperimentData, on_delete=models.CASCADE)
    x = models.CharField(max_length=32)
    y = models.CharField(max_length=32)

    class Meta:
        unique_together = ('experiment_data', 'x',)
