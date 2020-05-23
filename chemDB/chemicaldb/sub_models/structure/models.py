import rdkit
from django.db import models
from django.db.models.signals import pre_save
from django.dispatch import receiver

from rdkit import Chem
import rdkit.Chem.Descriptors as rdk_descriptors


def _validate_structure_identifiers(structure):
    mol = structure.get_mol()
    if mol is None:
        return False
    valid = True

    smiles = rdkit.Chem.MolToSmiles(mol, isomericSmiles=True)
    if structure.iso_smiles:
        if structure.iso_smiles != smiles:
            valid = False
    else:
        structure.iso_smiles = smiles

    inchi = rdkit.Chem.MolToInchi(mol)
    if structure.standard_inchi:
        if structure.standard_inchi != inchi:
            valid = False
    else:
        structure.standard_inchi = inchi

    ikey = rdkit.Chem.MolToInchiKey(mol)
    if structure.inchi_key:
        if structure.inchi_key != ikey:
            valid = False
    else:
        structure.inchi_key = ikey
    return valid


def _validate_structure_mass(structure):
    mol = structure.get_mol()
    if mol is None:
        return False
    valid = True

    mass = rdk_descriptors.MolWt(mol)
    if structure.molar_mass:
        if structure.molar_mass != mass:
            valid = False
    else:
        structure.molar_mass = mass

    return valid


def _validate_structure_mol(structure):
    return structure.get_mol() is not None


VALIDATION_HAS_MOL = "structure_has_mol", _validate_structure_mol
VALIDATION_IDENTIFIER = "structure_identifiers", _validate_structure_identifiers
VALIDATION_MASS = "structure_mass", _validate_structure_mass


class Structure(models.Model):
    class Meta:
        permissions = (
            ('add structure', 'Add Structure'),
        )

    iso_smiles = models.CharField(max_length=200, null=True, blank=True)
    standard_inchi = models.TextField(null=True, blank=True)
    inchi_key = models.CharField(max_length=27, null=True, blank=True)
    molar_mass = models.FloatField(null=True, blank=True, )
    valid = models.BooleanField(default=True, editable=False)

    # external references
    cas_number = models.CharField(max_length=12, unique=True, null=True, blank=True)
    chemspider_id = models.PositiveIntegerField(unique=True, null=True, blank=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.mol = None
        self.validity_checks = {
            VALIDATION_HAS_MOL[0]: VALIDATION_HAS_MOL[1],
            VALIDATION_IDENTIFIER[0]: VALIDATION_IDENTIFIER[1],
            VALIDATION_MASS[0]: VALIDATION_MASS[1],
        }

    def get_mol(self):
        if self.mol is None:
            if self.standard_inchi:
                self.mol = rdkit.Chem.MolFromInchi(self.standard_inchi)
            elif self.iso_smiles:
                self.mol = rdkit.Chem.MolFromSmiles(self.iso_smiles)
        return self.mol

    def validate(self, save=True):
        valid = True
        for check_name, check in self.validity_checks.items():
            valid = valid and check(self)
            if not valid:
                break
        if self.valid != valid:
            self.valid = valid
            if save:
                self.save()
        return self.valid

    def __str__(self):
        name = self.names.first()

        if self.cas_number:
            if name:
                return "{} ({})".format(self.cas_number, name)
            return self.cas_number
        if self.iso_smiles:
            if name:
                return "{} ({})".format(self.iso_smiles, name)
            return self.iso_smiles
        return super().__str__()


@receiver(pre_save)
def my_callback(sender, instance: Structure, *args, **kwargs):
    if not issubclass(sender, Structure):
        return
    instance.validate(save=False)


class StructureName(models.Model):
    structure = models.ForeignKey(Structure, on_delete=models.CASCADE, related_name="names")
    name = models.CharField(max_length=200)

    def __str__(self):
        return self.name


from .submodels import *
