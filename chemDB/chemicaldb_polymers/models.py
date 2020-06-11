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
        try:
            mol = MolFromSmiles(self.smiles * 4)
            # mc = Chem.Mol(mol.ToBinary())
            if not mol.GetNumConformers():
                rdDepictor.Compute2DCoords(mol)
            drawer = rdMolDraw2D.MolDraw2DSVG(200, 200)
            drawer.DrawMolecule(mol)
            drawer.FinishDrawing()
            svg = drawer.GetDrawingText()
            return mark_safe(svg)
        except:
            return ""


class MonomerDistribution(models.TextChoices):
    GRADIENT = 'GRADIENT', "gradient"
    STATISTICAL = 'STATISTICAL', "statistical"
    BLOCK = 'BLOCK', "block"
    MIXED = 'MIXED', "mixed"


MonomerDistribution.max_length = 12

class StartEndGroup(ValidationModel):
    structure = models.ForeignKey(Structure,on_delete=models.CASCADE)

    def __str__(self):
        return str(self.structure)

class RepeatingUnit(ValidationModel):
    structure = models.ForeignKey(Structure,on_delete=models.CASCADE)

    def __str__(self):
        return str(self.structure)


class Polymer(Substance):
    polymer_structure = models.ForeignKey(PolymerStructure, on_delete=models.CASCADE, related_name="polymer_structure")
    polymer_structure_type = models.CharField(max_length=16, choices=PolymerStructureType.choices,
                                              default=PolymerStructureType.LINEAR)
    monomers = models.ManyToManyField(RepeatingUnit, through='Monomers')
    monomer_distribution = models.CharField(max_length=MonomerDistribution.max_length,
                                            choices=MonomerDistribution.choices,
                                            default=MonomerDistribution.STATISTICAL)
    distribution_parameter = models.FloatField(default=1)
    start_group = models.ForeignKey(StartEndGroup, on_delete=models.SET_NULL, null=True,
                                    related_name="as_polymer_startgroup",blank=True)
    end_group = models.ForeignKey(StartEndGroup, on_delete=models.SET_NULL, null=True, related_name="as_polymer_endgroup",blank=True)


class StartEndGroupForm(forms.ModelForm):
    class Meta:
        model = Structure
        exclude = []


def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)


class Monomers(models.Model):
    polymer = models.ForeignKey(Polymer, on_delete=models.CASCADE, related_name="p")
    monomer = models.ForeignKey(RepeatingUnit, on_delete=models.CASCADE, )
    repeating_units = models.CharField(max_length=16, choices=MonomerDistribution.choices,
                                       default=MonomerDistribution.STATISTICAL)



class RelatedFieldWidgetCanAdd(widgets.SelectMultiple):

    def __init__(self, related_model, related_url=None, *args, **kw):

        super(RelatedFieldWidgetCanAdd, self).__init__(*args, **kw)

        if not related_url:
            rel_to = related_model
            info = (rel_to._meta.app_label, rel_to._meta.object_name.lower())
            related_url = 'admin:%s_%s_add' % info

        # Be careful that here "reverse" is not allowed
        self.related_url = related_url

    def render(self, name, value, *args, **kwargs):
        self.related_url = reverse(self.related_url)
        output = [super(RelatedFieldWidgetCanAdd, self).render(name, value, *args, **kwargs)]
        output.append('<a href="%s" class="add-another" id="add_id_%s" onclick="return showAddAnotherPopup(this);"> ' % \
                      (self.related_url, name))
        output.append('<img src="%sadmin/img/icon_addlink.gif" width="10" height="10" alt="%s"/></a>' % (settings.STATIC_URL, 'Add Another'))
        return mark_safe(''.join(output))

class NewPolymerForm(forms.ModelForm):

    code_prefix = forms.CharField(initial="hello",
                                  disabled=True,
                                  label="Prefix",
                                  required=False,

                             )
    class Meta:
        model = Polymer
        fields = ['name','code','polymer_structure_type','monomer_distribution','distribution_parameter','start_group','end_group']

    def __init__(self, chem_db_user, changeable=False, *args, **kwargs):

        super().__init__(*args, **kwargs)
        self.helper = FormHelper()
        if changeable:
            self.helper.form_method = 'post'
            self.helper.add_input(Submit('submit', 'save'))

        self.fields['code_prefix'].initial=chem_db_user.get_prefix()
        self.helper.layout = Layout(
            'name',
            Row(
                Column(
                    'code_prefix',
                    css_class='form-group col-md-2 mb-0'),
                Column(
                    'code', css_class='form-group col-md-10 mb-0'),
                css_class='form-row'
            ),
        'polymer_structure_type','monomer_distribution','distribution_parameter','start_group','end_group'
        )
#        initial={'owner':chem_db_user.pk,"short_name":}
#   self.fields['raw_data'] = forms.FileInput()
