from django.contrib.auth.models import User
from django.db import models


class ChemdbInstitute(models.Model):
    name=models.CharField(max_length=64)
    short=models.CharField(max_length=4,unique=True)

    def __str__(self):
        return self.name

class ChemdbUser(models.Model):
    auth_user = models.ForeignKey(User, on_delete=models.CASCADE)
    institute = models.ForeignKey(ChemdbInstitute,default=None,null=True,on_delete=models.SET_NULL,blank=True)
    short=models.CharField(max_length=4,unique=True)

    class Meta:
        unique_together=("short","institute")

    def __str__(self):
        if self.auth_user.first_name and self.auth_user.last_name:
            name = "{} {}".format(self.auth_user.first_name, self.auth_user.last_name)
            if self.institute:
                name+=" ({})".format(self.institute.short)
            return name

        return self.auth_user.username