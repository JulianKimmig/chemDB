from crispy_forms.helper import FormHelper
from crispy_forms.layout import Submit
from django import forms
from django.db import models

# Create your models here.
from rdkit.Chem import rdDepictor, MolFromSmiles
from rdkit.Chem.Draw import rdMolDraw2D

from chemicaldb.models import Structure, Substance, mark_safe, ValidationModel


class PolymerStructureType(models.TextChoices):
    LINEAR = 'LINEAR', "linear"
    BRANCHED = 'BRANCHED', "branched"
    NETWORK = 'NETWORK', "network"


class PolymerStructure(Structure):
    big_smiles = models.CharField(max_length=200, null=True, blank=True)

    def __init__(self, *args, **kwargs):
        self._meta.get_field('iso_smiles').default = False
        super().__init__(*args, **kwargs)

    def structure_image(self):
        mol = MolFromSmiles(self.smiles * 4)
        # mc = Chem.Mol(mol.ToBinary())
        if not mol.GetNumConformers():
            rdDepictor.Compute2DCoords(mol)
        drawer = rdMolDraw2D.MolDraw2DSVG(200, 200)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        return mark_safe(svg)


class MonomerDistribution(models.TextChoices):
    GRADIENT = 'GRADIENT', "gradient"
    STATISTICAL = 'STATISTICAL', "statistical"
    BLOCK = 'BLOCK', "block"
    MIXED = 'MIXED', "mixed"


MonomerDistribution.max_length = 12

class StartEndGroup(ValidationModel):
    structure = models.ForeignKey

    class Meta:
        abstract = True


class StartGroup(StartEndGroup):
    pass


class EndGroup(StartEndGroup):
    pass

class Polymer(Substance):
    polymer_structure = models.ForeignKey(PolymerStructure, on_delete=models.CASCADE, related_name="polymer_structure")
    polymer_structure_type = models.CharField(max_length=16, choices=PolymerStructureType.choices,
                                              default=PolymerStructureType.LINEAR)
    monomers = models.ManyToManyField(Structure, through='Monomers')
    monomer_distribution = models.CharField(max_length=MonomerDistribution.max_length,
                                            choices=MonomerDistribution.choices,
                                            default=MonomerDistribution.STATISTICAL)
    distribution_parameter = models.FloatField(default=1)
    start_group = models.ForeignKey(StartGroup, on_delete=models.SET_NULL, null=True,
                                    related_name="as_polymer_startgroup")
    end_group = models.ForeignKey(EndGroup, on_delete=models.SET_NULL, null=True, related_name="as_polymer_endgroup")


class StartEndGroupForm(forms.ModelForm):
    class Meta:
        model = Structure
        exclude = []


def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)


class Monomers(models.Model):
    polymer = models.ForeignKey(Polymer, on_delete=models.CASCADE, related_name="p")
    monomer = models.ForeignKey(Structure, on_delete=models.CASCADE, )
    repeating_units = models.CharField(max_length=16, choices=MonomerDistribution.choices,
                                       default=MonomerDistribution.STATISTICAL)


class NewPolymerForm(forms.ModelForm):
    class Meta:
        model = Polymer
        exclude = []

    def __init__(self, chem_db_user=None, changeable=False, *args, **kwargs):

        super().__init__(*args, **kwargs)
        self.helper = FormHelper()
        if changeable:
            self.helper.form_method = 'post'
            self.helper.add_input(Submit('submit', 'save'))
        if chem_db_user:
            self.fields["user"].initial = chem_db_user.pk
        self.fields["user"].widget = forms.HiddenInput()
#        initial={'owner':chem_db_user.pk,"short_name":}
#   self.fields['raw_data'] = forms.FileInput()
