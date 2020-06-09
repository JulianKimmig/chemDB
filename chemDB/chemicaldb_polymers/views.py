from django.forms import modelform_factory
from django.shortcuts import render

# Create your views here.
from django.views import View

from chemicaldb.sub_models import Structure, SubstanceForm, SimpleSubstance, SimpleSubstanceForm, StructureForm, \
    StructureSmilesForm
from chemicaldb_polymers.models import NewPolymerForm


def main(request):
    return render(request,"chemicaldb_polymers/main.html")


class NewPolymer(View):
    def get(self,request):
        chem_db_user = request.user.chemdbuser
        context = {"polymer_form": NewPolymerForm(chem_db_user=chem_db_user,changeable=True)}
        return render(request, "chemicaldb_polymers/new_polymer.html", context)


class StartEndGroup(View):
    def get(self,request,pk=None):
        chem_db_user = request.user.chemdbuser

        context={"form":StructureSmilesForm(changeable=True)}
        return render(request,"chemicaldb/structure.html",context=context)