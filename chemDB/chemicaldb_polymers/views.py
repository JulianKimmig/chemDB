from django.forms import modelform_factory
from django.shortcuts import render, redirect

# Create your views here.
from django.views import View

from chemicaldb.forms import StructureSmilesForm
from chemicaldb.services import mol_to_image_mol
from chemicaldb.sub_models import Structure, SimpleSubstance, StructureName, mark_safe
from chemicaldb_polymers.forms import PolymerForm
from chemicaldb_polymers.models import StartEndGroup, RepeatingUnit, PolymerStructure


def main(request):
    return render(request, "chemicaldb_polymers/main.html")

def structures_main(request):
    return render(request, "chemicaldb_polymers/structures_main.html")

class NewPolymerStructure(View):
    def get(self, request):
        chem_db_user = request.user.chemdbuser
        context = {}
        return render(request, "chemicaldb_polymers/new_polymer_structure.html", context)

    def post(self,request):
        chem_db_user = request.user.chemdbuser

        repeating_units_pks = request.POST.getlist("repeating_unit")
        start_group = request.POST.get("start_group")
        end_group = request.POST.get("end_group")
        names = [n.strip() for n in request.POST.get("names").strip().split('\n')]
        valid=True
        error=""
        if start_group:
            start_group = StartEndGroup.objects.get(pk=start_group)
            if not start_group.check_can_view(chem_db_user):
                valid=False
                start_group=None
                error+="invalid start group\n"


        if end_group:
            end_group = StartEndGroup.objects.get(pk=end_group)
            if not end_group.check_can_view(chem_db_user):
                valid=False
                end_group=None
                error+="invalid end group\n"

        repeating_units = []
        for repeating_unit in repeating_units_pks:
            repeating_unit = RepeatingUnit.objects.get(pk=repeating_unit)
            if not repeating_unit.check_can_view(chem_db_user):
                valid=False
                error+="invalid repeating unit\n"
            else:
                repeating_units.append(repeating_unit)


        if len(repeating_units_pks)==0:
            valid=False
            error+="no valid repeating units\n"

        if valid:
            structure = PolymerStructure.create_new_chemdbshare(chemdb_user=chem_db_user)
            names = [StructureName(structure=structure, name=n) for n in names if len(n) > 0]
            StructureName.objects.bulk_create(names, ignore_conflicts=True)
            for repeating_unit in repeating_units:
                structure.repeating_units.add(repeating_unit)
            if end_group:
                structure.terminal_end_group=end_group
            if start_group:
                structure.terminal_start_group=start_group
            structure.save()
        else:
            structure = PolymerStructure()
        mol = structure.create_mol(repeating_units=repeating_units, terminal_start_group=start_group,
                                        terminal_end_group=end_group, with_chain=True)

        context={"error":error,
                 "start_group":start_group,
                 "end_group":end_group,
                 "repeating_units":repeating_units,
                 "img": mark_safe(mol_to_image_mol(mol, "svg")),
                 "names":'\n'.join((str(name) for name in names))
                 }
        return render(request, "chemicaldb_polymers/new_polymer_structure.html", context)

def structure_browse(request):
    chem_db_user = request.user.chemdbuser
    structures=[]
    for structure in PolymerStructure.objects.all():
        if structure.check_can_view(chem_db_user):
            structures.append(structure)
    context={"structures":structures}

    return render(request, "chemicaldb_polymers/browse_structures.html",context=context)


def structure_view(request,pk):
    chem_db_user = request.user.chemdbuser
    structure= PolymerStructure.objects.get(pk=pk)
    assert structure.check_can_view(chem_db_user)
    context={"structure":structure}
    return render(request,"chemicaldb_polymers/structure_view.html",context=context)

class NewPolymer(View):
    def get(self, request):
        chem_db_user = request.user.chemdbuser

        try:
            intital_structure = request.GET.get("structure")
            intital_structure=PolymerStructure.objects.get(pk=intital_structure)
            if not intital_structure.check_can_view(chem_db_user):
                intital_structure=None
        except:
            intital_structure=None

        form = PolymerForm(chem_db_user=chem_db_user)
        context = {"polymer_form": form,
                   'intital_structure':intital_structure
                   }
        return render(request, "chemicaldb_polymers/new_polymer.html", context)

    def post(self, request):
        chem_db_user = request.user.chemdbuser
        form = PolymerForm(data=request.POST,chem_db_user=chem_db_user)

        try:
            intital_structure = request.POST.get("polymer_structure")
            intital_structure=PolymerStructure.objects.get(pk=intital_structure)
            if not intital_structure.check_can_view(chem_db_user):
                intital_structure=None
        except:
            intital_structure=None

        if form.is_valid():
            pass
        context = {"polymer_form": form,
                   'intital_structure':intital_structure
                   }
        return render(request, "chemicaldb_polymers/new_polymer.html", context)

class RepeatingUnitEdit(View):
    def get(self, request, pk=None):
        chem_db_user = request.user.chemdbuser
        instance = None
        if pk is not None:
            instance = RepeatingUnit.objects.get(pk=pk)

        form = StructureSmilesForm(instance=instance, changeable=True,submodel=RepeatingUnit)
        context = {"form": form,
                   "structure_model":"{}.{}".format(RepeatingUnit._meta.app_label,RepeatingUnit.__name__),
                   }
        return render(request, "chemicaldb/structure_edit.html", context=context)

    def post(self, request, pk=None):
        chem_db_user = request.user.chemdbuser
        instance = None

        if pk is not None:
            instance = RepeatingUnit.objects.get(pk=pk)

        form = StructureSmilesForm(data=request.POST, instance=instance, changeable=True)
        if form.is_valid():
            names = [n.strip() for n in form.cleaned_data.get("names").strip().split('\n')]
            instance = form.save(commit=False)
            instance.iso_smiles = False
            instance.save()

            names = [StructureName(structure=instance, name=n) for n in names if len(n) > 0]
            StructureName.objects.bulk_create(names, ignore_conflicts=True)


            return redirect("chemicaldb:chemicaldb_polymers:edit_repeating_unit", pk=instance.pk)

        context = {"form": form,
                   "structure_model":"{}.{}".format(RepeatingUnit._meta.app_label,RepeatingUnit.__name__),
                   }
        return render(request, "chemicaldb/structure_edit.html", context=context)


class StartEndGroupEdit(View):
    def get(self, request, pk=None):
        chem_db_user = request.user.chemdbuser
        instance = None
        if pk is not None:
            instance = StartEndGroup.objects.get(pk=pk)

        form = StructureSmilesForm(instance=instance, changeable=True,submodel=StartEndGroup)
        context = {"form": form,
                   "structure_model":"{}.{}".format(StartEndGroup._meta.app_label,StartEndGroup.__name__),
                   }
        return render(request, "chemicaldb/structure_edit.html", context=context)

    def post(self, request, pk=None):
        chem_db_user = request.user.chemdbuser
        instance = None
        if pk is not None:
            instance = StartEndGroup.objects.get(pk=pk)

        form = StructureSmilesForm(data=request.POST, instance=instance, changeable=True,submodel=StartEndGroup)
        names = form.data.get("names").strip().split('\n')
        if form.is_valid():
            names = [n.strip() for n in form.cleaned_data.get("names").strip().split('\n')]
            instance = form.save(commit=False)
            instance.iso_smiles = False
            instance.save()

            names = [StructureName(structure=instance, name=n) for n in names if len(n) > 0]
            StructureName.objects.bulk_create(names, ignore_conflicts=True)


            return redirect("chemicaldb:chemicaldb_polymers:edit_start_end_group", pk=instance.pk)

        context = {"form": form,
                   "structure_model":"{}.{}".format(StartEndGroup._meta.app_label,StartEndGroup.__name__),
                   }
        return render(request, "chemicaldb/structure_edit.html", context=context)


