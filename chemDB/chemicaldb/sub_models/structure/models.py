import rdkit
from django.db import models
from django.db.models.signals import pre_save
from django.dispatch import receiver

from rdkit import Chem
import rdkit.Chem.Descriptors as rdk_descriptors


class Structure(models.Model):
    class Meta:
        permissions = (
            ('add structure', 'Add Structure'),
        )


    iso_smiles = models.CharField(max_length=200, null=True, blank=True)
    standard_inchi = models.TextField(null=True, blank=True)
    inchi_key = models.CharField(max_length=27, null=True, blank=True)
    molar_mass = models.FloatField(null=True, blank=True, )
    valid = models.BooleanField(default=True, )

    # external references
    cas_number = models.CharField(max_length=12, unique=True, null=True, blank=True)
    chemspider_id = models.PositiveIntegerField(unique=True, null=True, blank=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.mol = None

    def get_mol(self):
        if self.mol is None:
            if self.standard_inchi:
                self.mol = rdkit.Chem.MolFromInchi(self.standard_inchi)
            elif self.iso_smiles:
                self.mol = rdkit.Chem.MolFromSmiles(self.iso_smiles)
        return self.mol

    def _validate_structure_identifiers(self, mol):
        if mol is None:
            return False
        valid = True

        smiles = rdkit.Chem.MolToSmiles(mol, isomericSmiles=True)
        if self.iso_smiles:
            if self.iso_smiles != smiles:
                valid = False
        else:
            self.iso_smiles = smiles

        inchi = rdkit.Chem.MolToInchi(mol)
        if self.standard_inchi:
            if self.standard_inchi != inchi:
                valid = False
        else:
            self.standard_inchi = inchi

        ikey = rdkit.Chem.MolToInchiKey(mol)
        if self.inchi_key:
            if self.inchi_key != ikey:
                valid = False
        else:
            self.inchi_key = ikey
        return valid

    def _validate_properties(self, mol):
        if mol is None:
            return False
        valid = True

        mass = rdk_descriptors.MolWt(mol)
        if self.molar_mass:
            if self.molar_mass != mass:
                valid = False
        else:
            self.molar_mass = mass

        return valid

    def validate(self, save=True):
        mol = self.get_mol()
        valid = True
        if mol is None:
            valid = False
        else:
            valid = valid and self._validate_structure_identifiers(mol)
            valid = valid and self._validate_properties(mol)

        if self.valid != valid:
            self.valid = valid
            if save:
                self.save()
        return self.valid

    def __str__(self):
        name = self.names.first()

        if self.cas_number:
            if name:
                return "{} ({})".format(self.cas_number,name)
            return self.cas_number
        if self.iso_smiles:
            if name:
                return "{} ({})".format(self.iso_smiles,name)
            return self.iso_smiles
        return super().__str__()


@receiver(pre_save)
def my_callback(sender, instance: Structure, *args, **kwargs):
    if not issubclass(sender, Structure):
        return
    instance.validate(save=False)



class StructureName(models.Model):
    structure = models.ForeignKey(Structure, on_delete=models.CASCADE,related_name="names")
    name = models.CharField(max_length=200)

    def __str__(self):
        return self.name



from .submodels import *
