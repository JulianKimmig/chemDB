from django.db import models

from structure.models import Structure

class StructureName(models.Model):
    structure = models.ForeignKey(Structure, on_delete=models.CASCADE)
    name = models.CharField(max_length=200)

