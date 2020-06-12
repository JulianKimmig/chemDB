from django.forms import modelform_factory
from django.shortcuts import render, redirect

# Create your views here.
from django.views import View

from chemicaldb.forms import StructureSmilesForm
from chemicaldb.sub_models import Structure, SimpleSubstance, StructureName
from chemicaldb_polymers.models import NewPolymerForm, StartEndGroup, RepeatingUnit


def main(request):
    return render(request, "chemicaldb_polymers/main.html")


class NewPolymer(View):
    def get(self, request):
        chem_db_user = request.user.chemdbuser
        context = {"polymer_form": NewPolymerForm(chem_db_user=chem_db_user, changeable=True)}
        return render(request, "chemicaldb_polymers/new_polymer.html", context)


class RepeatingUnitEdit(View):
    def get(self, request, pk=None):
        chem_db_user = request.user.chemdbuser
        instance = None
        structure = None
        if pk is not None:
            instance = RepeatingUnit.objects.get(pk=pk)
            structure = instance.structure

        form = StructureSmilesForm(instance=structure, changeable=True)
        context = {"form": form}
        return render(request, "chemicaldb/structure_edit.html", context=context)

    def post(self, request, pk=None):
        chem_db_user = request.user.chemdbuser
        instance = None
        structure = None
        if pk is not None:
            instance = RepeatingUnit.objects.get(pk=pk)
            structure = instance.structure

        if structure is None:
            try:
                structure = Structure.objects.get(smiles=request.POST.get("smiles"))
            except:
                pass
        form = StructureSmilesForm(data=request.POST, instance=structure, changeable=True)
        if form.is_valid():
            names = [n.strip() for n in form.cleaned_data.get("names").strip().split('\n')]
            structure = form.save(commit=False)
            structure.iso_smiles = False
            structure.save()

            names = [StructureName(structure=structure, name=n) for n in names if len(n) > 0]
            StructureName.objects.bulk_create(names, ignore_conflicts=True)

            if instance is None:
                instance = RepeatingUnit.objects.create(structure=structure)

            return redirect("chemicaldb:chemicaldb_polymers:edit_repeating_unit", pk=instance.pk)

        context = {"form": form}
        return render(request, "chemicaldb/structure_edit.html", context=context)


class StartEndGroupEdit(View):
    def get(self, request, pk=None):
        chem_db_user = request.user.chemdbuser
        instance = None
        structure = None
        if pk is not None:
            instance = StartEndGroup.objects.get(pk=pk)
            structure = instance.structure

        form = StructureSmilesForm(instance=structure, changeable=True)
        context = {"form": form}
        return render(request, "chemicaldb/structure_edit.html", context=context)

    def post(self, request, pk=None):
        chem_db_user = request.user.chemdbuser
        instance = None
        structure = None
        if pk is not None:
            instance = StartEndGroup.objects.get(pk=pk)
            print(instance)
            structure = instance.structure

        if structure is None:
            try:
                structure = Structure.objects.get(smiles=request.POST.get("smiles"))
            except:
                pass
        form = StructureSmilesForm(data=request.POST, instance=structure, changeable=True)
        names = form.data.get("names").strip().split('\n')
        if form.is_valid():
            names = [n.strip() for n in form.cleaned_data.get("names").strip().split('\n')]
            structure = form.save(commit=False)
            structure.iso_smiles = False
            structure.save()

            names = [StructureName(structure=structure, name=n) for n in names if len(n) > 0]
            StructureName.objects.bulk_create(names, ignore_conflicts=True)

            if instance is None:
                instance = StartEndGroup.objects.create(structure=structure)

            return redirect("chemicaldb:chemicaldb_polymers:edit_start_end_group", pk=instance.pk)

        context = {"form": form}
        return render(request, "chemicaldb/structure_edit.html", context=context)
