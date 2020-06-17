from crispy_forms.bootstrap import PrependedText
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Submit, Layout
from django import forms

from chemicaldb.sub_models import Structure, Substance, SimpleSubstance, ChemDbShareModel


class ChemDbShareModelForm(forms.ModelForm):
    owner = forms.IntegerField()

    class Meta:
        model = ChemDbShareModel
        fields = ["code", "public", "public_to_user", "public_to_institute"]

    def __init__(self, chem_db_user, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.helper = FormHelper()
        self.helper.form_method = 'post'
        self.fields["owner"].initial = chem_db_user.pk
        self.fields["owner"].widget = forms.HiddenInput()
        self.helper.layout = Layout(*list(self.fields.keys()))
        if "code" in self.fields:
            self.helper['code'].wrap(PrependedText, chem_db_user.get_prefix())

        self.helper.add_input(Submit('submit', 'Save'))

class SubstanceForm(ChemDbShareModelForm):
    owner = forms.IntegerField()

    class Meta:
        model = Substance
        fields = ["name"]+ChemDbShareModelForm.Meta.fields



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

    def __init__(self, chem_db_user=None, changeable=False, submodel=Structure, *args, **kwargs):
        super().__init__(*args, **kwargs)

        if not issubclass(submodel, Structure):
            raise TypeError("{} not of sublcass Structure".format(submodel))
        self._meta.model = submodel

        self.helper = FormHelper()
        if self.instance:
            self.fields["names"].initial = "\n".join([n.name for n in self.instance.names.all()])
            self.fields["smiles"].help_text = submodel.SMILES_INFO
        if changeable:
            self.helper.form_method = 'post'
            self.helper.add_input(Submit('submit', 'save'))



#        initial={'owner':chem_db_user.pk,"short_name":}
#   self.fields['raw_data'] = forms.FileInput()


class SimpleSubstanceForm(SubstanceForm):
    class Meta:
        model = SimpleSubstance
        exclude = []
