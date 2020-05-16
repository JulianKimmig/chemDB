from django.contrib.auth.models import User
from django.db import models

from structure.models import Identifier
from substance.models import Substance


class Internal_Id(Identifier):
    user = models.ForeignKey(User,on_delete=models.SET_NULL,null=True)
    id = models.CharField(max_length=20,primary_key=True)
    substance = models.ForeignKey(Substance,on_delete=models.CASCADE)