from django import forms

from chemicaldb.forms import SubstanceForm
from chemicaldb_polymers.models import Polymer


class PolymerForm(SubstanceForm):
    class Meta:
        model = Polymer
        fields = ["polymer_structure"]+SubstanceForm.Meta.fields

        widgets = {
            'polymer_structure': forms.TextInput(attrs={'class':"search_model_autocomplete_input_dummy"}),
        }

#    def __init__(self, chem_db_user, *args, **kwargs):
 #       super().__init__(chem_db_user, *args, **kwargs)