import rdkit
from django.db import models

# Create your models here.
# structure specific
from rdkit import Chem


class Structure(models.Model):
    iso_smiles = models.CharField(max_length=200)
    standard_inchi = models.TextField()
    inchi_key = models.CharField(max_length=27)
    molar_mass = models.FloatField()
    valid = models.BooleanField(default=True)

    cas_number = models.CharField(max_length=12, primary_key=True)
    chemspider_id = models.PositiveIntegerField(primary_key=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def get_mol(self):
        if self.mol is None:
            self.mol = rdkit.Chem.MolFromInchi(self.standard_inchi)
        return self.mol

    def validate(self,save=True):
        mol = self.get_mol()
        valid=True
        if self.iso_smiles != rdkit.Chem.MolToSmiles(mol):
            valid = False
        if self.inchi_key != rdkit.Chem.MolToInchiKey(mol):
            valid = False
        if self.valid != valid:
            self.valid = valid
            if save:
                self.save()
        return self.valid

from .submodels import *
