from crispy_forms.helper import FormHelper
from crispy_forms.layout import Submit
from django import forms

from chemicaldb.sub_models import Structure, Substance, SimpleSubstance


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

class SubstanceForm(forms.ModelForm):
    class Meta:
        model = Substance
        exclude = []

    def __init__(self,chem_db_user=None,changeable=False, *args, **kwargs):

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


class SimpleSubstanceForm(SubstanceForm):
    class Meta:
        model = SimpleSubstance
        exclude = []
