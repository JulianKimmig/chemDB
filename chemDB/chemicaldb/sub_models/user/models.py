from django.contrib.auth.models import User
from django.db import models
from django.db.models.signals import post_save
from django.dispatch import receiver


class ChemdbInstitute(models.Model):
    name=models.CharField(max_length=64)
    short=models.CharField(max_length=4,unique=True)

class ChemdbUser(models.Model):
    auth_user = models.ForeignKey(User, on_delete=models.CASCADE)
    institute = models.ForeignKey(ChemdbInstitute,default=None,null=True,on_delete=models.SET_NULL)
    short=models.CharField(max_length=4,unique=True)

    class Meta:
        unique_together=("short","institute")