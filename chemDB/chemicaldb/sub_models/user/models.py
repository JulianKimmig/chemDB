from django.contrib.auth.models import User
from django.db import models
from django.db.models.signals import pre_save
from django.dispatch import receiver


class ChemdbInstitute(models.Model):
    name = models.CharField(max_length=64)
    short = models.CharField(max_length=4, unique=True)
    parent_institute = models.ForeignKey("self", on_delete=models.SET_NULL, null=True, related_name="child_institutes",
                                         blank=True)

    def __str__(self):
        return "{} ({})".format(self.name, self.get_short())

    def get_all_child(self, include_self=False, stop_at=None):
        children = []
        if include_self:
            children.append(self)
        for c in self.child_institutes.all():
            if stop_at:
                if stop_at(c):
                    children.append(c)
                    return children

            if not include_self:
                children.append(c)

            children.extend(c.get_all_child(include_self=include_self))
        return children

    def get_short(self):
        if self.parent_institute:
            return self.parent_institute.get_short() + "_" + self.short
        return self.short

    def is_parent(self, other_institute):
        for c in self.get_all_child(stop_at=lambda c: c.pk == other_institute.pk):
            if c.pk == other_institute.pk:
                return True
        return False


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
        prefix = self.get_prefix()
        if string.startswith(prefix):
            return string
        return "{}_{}".format(prefix, string)


class ChemDbShareModel(models.Model):
    PUBLIC_TO_USER = True
    PUBLIC = True
    owner = models.ForeignKey(ChemdbUser, on_delete=models.SET_NULL, null=True,related_name="%(app_label)s_%(class)s_owner")
    admin = models.ManyToManyField(ChemdbUser,related_name="%(app_label)s_%(class)s_admin")

    can_edit_user = models.ManyToManyField(ChemdbUser,related_name="%(app_label)s_%(class)s_can_edit")
    can_edit_institute = models.ManyToManyField(ChemdbInstitute,related_name="%(app_label)s_%(class)s_edit_institute")
    can_view_user = models.ManyToManyField(ChemdbUser,related_name="%(app_label)s_%(class)s_can_view")
    can_view_institute = models.ManyToManyField(ChemdbInstitute,related_name="%(app_label)s_%(class)s_view_institute")

    public_to_user = models.BooleanField()
    public=models.BooleanField()

    code = models.CharField(max_length=64)
    code_prefix = models.CharField(max_length=64)

    class Meta:
        abstract = True
        unique_together = ("code", "code_prefix")

    def check_can_view(self, chemdb_user: ChemdbUser):
        if self.owner == chemdb_user:
            return True
        if self.admin.filter(pk=chemdb_user.pk).exists():
            return True
        if self.can_view_user.filter(pk=chemdb_user.pk).exists():
            return True
        if self.can_view_institute.filter(pk=chemdb_user.institute.pk).exists():
            return True

        for i in self.can_view_institute.all():
            if i.is_parent(chemdb_user.institute):
                return True

        return False

    def check_can_edit(self, chemdb_user: ChemdbUser):
        if self.owner == chemdb_user:
            return True
        if self.admin.filter(pk=chemdb_user.pk).exists():
            return True
        if self.can_edit_user.filter(pk=chemdb_user.pk).exists():
            return True
        if self.can_edit_institute.filter(pk=chemdb_user.institute.pk).exists():
            return True

        for i in self.can_edit_institute.all():
            if i.is_parent(chemdb_user.institute):
                return True

        return False

    def __str__(self):
        return "{}_{}".format(self.code_prefix, self.code)


@receiver(pre_save)
def code_checker(sender, instance: ChemDbShareModel, *args, **kwargs):
    if not issubclass(sender, ChemDbShareModel):
        return
    if instance.owner:
        instance.code_prefix = instance.owner.get_prefix()
    if instance.public_to_user is None:
        instance.public_to_user=instance.PUBLIC_TO_USER
    if instance.public is None:
        instance.public=instance.PUBLIC