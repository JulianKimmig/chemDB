from django.db import models

# Create your models here.
import chemDB
print(chemDB)
from chemDB.experiments.models import Experiment


class Source(models.Model):
    pass

class ExperimentsSource(Source):
    experiment = models.ForeignKey(Experiment,on_delete=models.CASCADE)