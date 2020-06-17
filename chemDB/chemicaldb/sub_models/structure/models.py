import base64

import rdkit
import rdkit.Chem.Descriptors
import rdkit.Chem.Draw
from django.db import models
from django.utils.safestring import mark_safe
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

from .. import ChemDbShareModel
from ...models import ValidationModel, ValidationError
from ...services import structure_image, mol_to_image_mol


def _validate_structure_identifiers(structure):
    mol = structure.get_mol()
    if mol is None:
        raise ValidationError("generated mol is None")
    valid = True

    smiles = rdkit.Chem.MolToSmiles(mol, isomericSmiles=True)
    if structure.smiles:
        gen_smiles=rdkit.Chem.MolToSmiles(rdkit.Chem.MolFromSmiles(structure.smiles), isomericSmiles=True)
        if gen_smiles != smiles:
            valid = False
            reason = "smiles doesn't match ({}, {})".format(smiles,gen_smiles)
    else:
        structure.smiles = smiles

    inchi = rdkit.Chem.MolToInchi(mol)
    if structure.standard_inchi:
        if structure.standard_inchi != inchi:
            valid = False
            reason = "inchi doesn't match"
    else:
        structure.standard_inchi = inchi

    ikey = rdkit.Chem.MolToInchiKey(mol)
    if structure.inchi_key:
        if structure.inchi_key != ikey:
            reason = "inchikey doesn't match"
            valid = False
    else:
        structure.inchi_key = ikey
    if not valid:
        raise ValidationError(reason)
    return valid


def _validate_structure_mass(structure):
    mol = structure.get_mol()
    if mol is None:
        return False
    valid = True

    mass = rdkit.Chem.Descriptors.MolWt(mol)
    if structure.molar_mass:
        if structure.molar_mass != mass:
            valid = False
            raise ValidationError("mass does not match")
    else:
        structure.molar_mass = mass

    return valid


def _validate_structure_mol(structure):
    if structure.get_mol() is not None:
        return True
    else:
        raise ValidationError("generated mol is None")


def _validate_structre_iso_identifier(structure):
    if structure.smiles and structure.iso_smiles:
        smile_mol = rdkit.Chem.MolFromSmiles(structure.smiles)
        new_smiles = rdkit.Chem.MolToSmiles(smile_mol, isomericSmiles=True)
        structure.smiles = new_smiles
    if structure.standard_inchi:
        inchi_mol = rdkit.Chem.MolFromInchi(structure.standard_inchi)
        new_inchi = rdkit.Chem.MolToInchi(inchi_mol)
        structure.standard_inchi = new_inchi
    return True


VALIDATION_ISO_IDENTIFIER = "iso_identifier", _validate_structre_iso_identifier
VALIDATION_HAS_MOL = "structure_has_mol", _validate_structure_mol
VALIDATION_IDENTIFIER = "structure_identifiers", _validate_structure_identifiers
VALIDATION_MASS = "structure_mass", _validate_structure_mass


class Structure(ValidationModel, ChemDbShareModel):
    PUBLIC_TO_USER = True
    PUBLIC = True

    SMILES_INFO=""

    class Meta:
        permissions = (
            ('add structure', 'Add Structure'),
        )

    smiles = models.CharField(max_length=200, null=True, blank=True)
    standard_inchi = models.TextField(null=True, blank=True)
    inchi_key = models.CharField(max_length=27, null=True, blank=True)
    molar_mass = models.FloatField(null=True, blank=True, )
    iso_smiles = models.BooleanField(default=True, editable=False)

    # external references
    cas_number = models.CharField(max_length=12, unique=True, null=True, blank=True)
    chemspider_id = models.PositiveIntegerField(unique=True, null=True, blank=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.mol = None
        self.validity_checks.update({
            VALIDATION_ISO_IDENTIFIER[0]: VALIDATION_ISO_IDENTIFIER[1],
            VALIDATION_HAS_MOL[0]: VALIDATION_HAS_MOL[1],
            VALIDATION_IDENTIFIER[0]: VALIDATION_IDENTIFIER[1],
            VALIDATION_MASS[0]: VALIDATION_MASS[1],
        })

    def get_mol(self):
        if self.mol is None:
            self.mol = self.create_mol()
        if self.mol is None:
            raise ValidationError("cannot generate mol for molecule")
        return self.mol

    def create_mol(self, smiles=None, inchi=None):
        mol = None
        if smiles is not None:
            mol = rdkit.Chem.MolFromSmiles(smiles)
        if mol is not None:
            return mol
        if inchi is not None:
            mol = rdkit.Chem.MolFromInchi(self.smiles)
        if mol is not None:
            return mol

        if self.smiles:
            mol = rdkit.Chem.MolFromSmiles(self.smiles)
        elif self.standard_inchi:
            mol = rdkit.Chem.MolFromInchi(self.standard_inchi)
        return mol

    def structure_image_field(self,size=200):
        mol = self.create_mol()
        if mol is None:
            raise ValidationError("cannot get mol for molecule")
        return mark_safe(mol_to_image_mol(mol, "svg",size=size))

    def structure_image_field_small(self):
        return self.structure_image_field(size=100)

    def __str__(self):
        name = self.names.first()

        if self.cas_number:
            if name:
                return "{} ({})".format(self.cas_number, name)
            return self.cas_number
        if self.smiles:
            if name:
                return "{} ({})".format(self.smiles, name)
            return self.smiles

        return str(name)

    def get_all_names(self):
        return [n for n in self.names.values_list("name",flat=True)]

    all_names=property(get_all_names)

class StructureName(models.Model):
    structure = models.ForeignKey(Structure, on_delete=models.CASCADE, related_name="names")
    name = models.CharField(max_length=200)

    def __str__(self):
        return self.name

    class Meta:
        unique_together = ("structure", "name")

