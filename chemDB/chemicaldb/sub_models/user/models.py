from django.contrib.auth.models import User
from django.db import models


class ChemdbInstitute(models.Model):
    name = models.CharField(max_length=64)
    short = models.CharField(max_length=4, unique=True)
    parent_institute = models.ForeignKey("self",on_delete=models.SET_NULL,null=True,related_name="child_institute",blank=True)

    def __str__(self):
        return "{} ({})".format(self.name,self.get_short())

    def get_short(self):
        if self.parent_institute:
            return self.parent_institute.get_short() + "_" + self.short
        return self.short

class ChemdbUser(models.Model):
    auth_user = models.OneToOneField(User, on_delete=models.CASCADE, primary_key=True)
    institute = models.ForeignKey(ChemdbInstitute, default=None, null=True, on_delete=models.SET_NULL, blank=True)
    short = models.CharField(max_length=4, unique=True)

    class Meta:
        unique_together = ("short", "institute")

    def __str__(self):
        if self.auth_user.first_name and self.auth_user.last_name:
            name = "{} {}".format(self.auth_user.first_name, self.auth_user.last_name)
            if self.institute:
                name += " ({})".format(self.institute.get_short())
            return name

        return self.auth_user.username

    def get_prefix(self):
        s = ""
        if self.institute:
            s += self.institute.get_short() + "_"
        s += self.short
        return s

    def prefix_string(self, string):
        prefix=self.get_prefix()
        if string.startswith(prefix):
            return string
        return "{}_{}".format(prefix, string)
