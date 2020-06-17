import rdkit
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Submit, Layout, Row, Column
from django import forms
from django.db import models

# Create your models here.
from django.forms import widgets, TextInput
from django.urls import reverse
from rdkit.Chem import rdDepictor, MolFromSmiles
from rdkit.Chem.Draw import rdMolDraw2D

from chemDB import settings
from chemicaldb.models import Structure, Substance, mark_safe, ValidationModel, VALIDATION_IDENTIFIER, ValidationError, \
    SubstanceProperty
from chemicaldb.services import mol_to_image_mol


def _validate_repeating_units(structure):
    mol = structure.get_mol()
    if mol is None:
        raise ValidationError("generated mol is None")
    assert rdkit.Chem.Descriptors.NumRadicalElectrons(mol) >= 2, "repeating units need at least two radicals"
    valid = True
    return valid


def _validate_start_end_groups(structure):
    mol = structure.get_mol()
    if mol is None:
        raise ValidationError("generated mol is None")
    assert rdkit.Chem.Descriptors.NumRadicalElectrons(mol) >= 1, "Start and end groups need at least one radical"
    valid = True
    return valid


VALIDATION_REPEATING_UNIT = "repeating_unit_validation", _validate_repeating_units
VALIDATION_START_END_GROUP = "start_end_group_validation", _validate_start_end_groups


class PolymerStructureType(models.TextChoices):
    LINEAR = 'LINEAR', "linear"
    BRANCHED = 'BRANCHED', "branched"
    NETWORK = 'NETWORK', "network"


class StartEndGroup(Structure):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.validity_checks.update({
            VALIDATION_START_END_GROUP[0]: VALIDATION_START_END_GROUP[1]
        })
        del self.validity_checks[VALIDATION_IDENTIFIER[0]]

    def structure_image_field(self,size=200):
        mol = self.create_mol(with_chain=True)
        if mol is None:
            raise ValidationError("cannot get mol for molecule")
        return mark_safe(mol_to_image_mol(mol, "svg",size=size))

    def create_mol(self, smiles=None, with_chain=False):
        if smiles is None:
            smiles = self.smiles
        if with_chain:
            smiles += "[*]"
        return super().create_mol(smiles=smiles)


class RepeatingUnit(Structure):
    big_smiles = models.CharField(max_length=200, null=True, blank=True)
    SMILES_INFO = "smiles should contain at least two radicals in the form of [*]"

    def structure_image_field(self,size=200):
        mol = self.create_mol(with_chain=True)
        if mol is None:
            raise ValidationError("cannot get mol for molecule")
        return mark_safe(mol_to_image_mol(mol, "svg",size=size))

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.validity_checks.update({
            VALIDATION_REPEATING_UNIT[0]: VALIDATION_REPEATING_UNIT[1]
        })
        del self.validity_checks[VALIDATION_IDENTIFIER[0]]

    def create_mol(self, smiles=None, n=4, with_chain=False):
        if smiles is None:
            smiles = self.smiles
        smiles = smiles * n
        if with_chain:
            smiles = "[*]" + smiles + "[*]"
        return super().create_mol(smiles=smiles)

class MonomerDistribution(models.TextChoices):
    GRADIENT = 'GRADIENT', "gradient"
    STATISTICAL = 'STATISTICAL', "statistical"
    BLOCK = 'BLOCK', "block"
    MIXED = 'MIXED', "mixed"


MonomerDistribution.max_length = 12


class PolymerStructure(RepeatingUnit):
    repeating_units = models.ManyToManyField(RepeatingUnit, related_name="repeating_units")

    terminal_start_group = models.ForeignKey(StartEndGroup, on_delete=models.SET_NULL, null=True,
                                             related_name="as_polymer_startgroup", blank=True)
    terminal_end_group = models.ForeignKey(StartEndGroup, on_delete=models.SET_NULL, null=True,
                                           related_name="as_polymer_endgroup", blank=True)

    monomer_distribution = models.CharField(max_length=MonomerDistribution.max_length,
                                            choices=MonomerDistribution.choices,
                                            default=MonomerDistribution.STATISTICAL)

    polymer_structure_type = models.CharField(max_length=16, choices=PolymerStructureType.choices,
                                              default=PolymerStructureType.LINEAR)


    def __init__(self, *args, **kwargs):
        self._meta.get_field('iso_smiles').default = False
        super().__init__(*args, **kwargs)
        del self.validity_checks[VALIDATION_REPEATING_UNIT[0]]

    def create_mol(self, repeating_units=None, terminal_start_group=None, terminal_end_group=None, with_chain=False):
        smiles=""
        if repeating_units is None:
            repeating_units = list(self.repeating_units.all())
        if terminal_start_group is None:
            terminal_start_group = self.terminal_start_group
        if terminal_end_group is None:
            terminal_end_group = self.terminal_end_group

        if terminal_start_group:
            smiles+=terminal_start_group.smiles
        elif with_chain:
            smiles+="[*]"

        for ru in repeating_units:
            smiles+=ru.smiles

        if terminal_end_group:
            smiles+=terminal_end_group.smiles
        elif with_chain:
            smiles+="[*]"

        return super().create_mol(smiles=smiles,n=1)


class StartEndGroupForm(forms.ModelForm):
    class Meta:
        model = Structure
        exclude = []


def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)


class Polymer(Substance):
    polymer_structure = models.ForeignKey(PolymerStructure, on_delete=models.CASCADE, related_name="polymer_structure")

class PolymerPropertie(SubstanceProperty):
    substance = models.ForeignKey(Polymer, on_delete=models.CASCADE)

    class Meta:
        abstract = True

class GlasTransition(PolymerPropertie):
    value = models.PositiveIntegerField()

class NumbeAverageMolarMass(PolymerPropertie):
    value = models.PositiveIntegerField()

class MassAverageMolarMass(PolymerPropertie):
    value = models.PositiveIntegerField()