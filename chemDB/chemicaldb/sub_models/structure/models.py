import rdkit
import rdkit.Chem.Descriptors as rdk_descriptors
import rdkit.Chem.Draw
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Submit
from django import forms
from django.db import models
from django.utils.safestring import mark_safe
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

from .. import ChemDbShareModel
from ...models import ValidationModel


def _validate_structure_identifiers(structure):
    mol = structure.get_mol()
    if mol is None:
        return False
    valid = True

    smiles = rdkit.Chem.MolToSmiles(mol, isomericSmiles=True)
    if structure.smiles:
        if rdkit.Chem.MolToSmiles(rdkit.Chem.MolFromSmiles(structure.smiles), isomericSmiles=True) != smiles:
            valid = False
    else:
        structure.smiles = smiles

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
            if self.smiles:
                self.mol = rdkit.Chem.MolFromSmiles(self.smiles)
            elif self.standard_inchi:
                self.mol = rdkit.Chem.MolFromInchi(self.standard_inchi)
        return self.mol

    def structure_image(self):
        mol = self.get_mol()
        # mc = Chem.Mol(mol.ToBinary())
        if not mol.GetNumConformers():
            rdDepictor.Compute2DCoords(mol)
        drawer = rdMolDraw2D.MolDraw2DSVG(200, 200)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        return mark_safe(svg)

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
        return super().__str__()


class StructureName(models.Model):
    structure = models.ForeignKey(Structure, on_delete=models.CASCADE, related_name="names")
    name = models.CharField(max_length=200)

    def __str__(self):
        return self.name

    class Meta:
        unique_together = ("structure", "name")


class StructureForm(forms.ModelForm):
    class Meta:
        model = Structure
        exclude = []

    def __init__(self, chem_db_user=None, changeable=False, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.helper = FormHelper()
        if changeable:
            self.helper.form_method = 'post'
            self.helper.add_input(Submit('submit', 'save'))


class StructureSmilesForm(forms.ModelForm):
    names = forms.CharField(widget=forms.Textarea, required=False)

    class Meta:
        model = Structure
        fields = ['smiles']

    def __init__(self, chem_db_user=None, changeable=False, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.helper = FormHelper()
        if self.instance:
            self.fields["names"].initial = "\n".join([n.name for n in self.instance.names.all()])

        if changeable:
            self.helper.form_method = 'post'
            self.helper.add_input(Submit('submit', 'save'))
