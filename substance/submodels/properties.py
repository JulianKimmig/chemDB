from django.db import models

from substance.models import Substance


class SubstanceProperty(models.Model):
    substance = models.ForeignKey(Substance, on_delete=models.CASCADE)
    class Meta:
        abstract = True

class MeltingPoint(SubstanceProperty):
    value = models.FloatField()

class BoilingPoint(SubstanceProperty):
    value = models.FloatField()

