from django.db import models

from chemDB.sources import Source
from ..models import Substance


class SubstanceProperty(models.Model):
    substance = models.ForeignKey(Substance, on_delete=models.CASCADE)
    value = models.Field()
    source = models.ForeignKey(Source,on_delete=models.SET_NULL,null=True)
    false_flag = models.BooleanField(default=False)

    class Meta:
        abstract = True

class MeltingPoint(SubstanceProperty):
    value = models.FloatField()


class BoilingPoint(SubstanceProperty):
    value = models.FloatField()
