from datetime import datetime

from crispy_forms.layout import Hidden
from django.core import serializers
from django.db.models import Count
from django.forms import modelformset_factory
from django.shortcuts import render, redirect

# Create your views here.
from django.utils.safestring import mark_safe
from django.views import View

from chemicaldb.sub_models import Substance
from experiments.models import ExperimentRawData, ExperimentData, DataTypesChoices, ExperimentSingleValueData, \
    ConcentrationUnits
from experiments_nanoparticle.characterization_tools import NanoparticleCharacterizationTool
from experiments_nanoparticle.models import NanoparticleBatchCharacterizationForm, Nanoparticle, NanoparticleForm, \
    NanoparticleCharacterization, NanoparticleCharacterizationForm, NanoparticlePreparationMethodForm, \
    NanoparticlePreparationMethod, Materials, Additives


def main_view(request):
    return render(request, "experiments_nanoparticle/main.html")


class batch_particle_creation(View):
    def get(self, request):
        chem_db_user = request.user.chemdbuser
        characterization_form = NanoparticleBatchCharacterizationForm(chem_db_user=chem_db_user)
        characterization_form.helper.add_input(Hidden("next_step", 1))
        context = {'characterization_form': characterization_form}
        return render(request, "experiments_nanoparticle/batch_create_1.html", context)

    def post(self, request):
        step = int(request.POST.get("next_step", "1"))
        if step == 1:
            return self.post1(request)
        if step == 2:
            return self.post2(request)
        if step == 3:
            return self.post3(request)
        return self.get(request)

    def post1(self, request):
        chem_db_user = request.user.chemdbuser
        characterization_form = NanoparticleBatchCharacterizationForm(chem_db_user=chem_db_user, data=request.POST,
                                                                      files=request.FILES, )
        if characterization_form.is_valid():
            form_data = characterization_form.cleaned_data

            tool = getattr(NanoparticleCharacterizationTool, form_data['tool'])
            df = tool.read_batch_data(data=form_data['exported_data'].read(), name=form_data['name'],
                                      short_name=form_data['short_name'])

            batch_experiment = characterization_form.save(commit=True)
            batch_experiment.owner = chem_db_user

            ExperimentRawData.objects.create(name="exported_data", raw_data=form_data['exported_data'],
                                             experiment=batch_experiment)
            ExperimentRawData.objects.create(name="raw_data", raw_data=form_data['raw_data'],
                                             experiment=batch_experiment)

            earliest = df['run_date'].min()
            nanoparticles = []
            characterizations = []
            step = ""
            try:
                step = "NanoparticleCharacterization"
                for row, data in df.iterrows():
                    #   if earliest is None:
                    #       earliest = data['run_date']
                    #   if data['run_date'] < earliest:
                    #       earliest = data['run_date']

                    character = NanoparticleCharacterization.objects.create(
                        name="{} in {}".format(data['sample_name'], batch_experiment.short_name),
                        short_name="{}_np-characterization".format(data['sample_name']),
                        owner=chem_db_user,
                        tool=batch_experiment.tool,
                        run_date=data['run_date'],
                        batch_experiment=batch_experiment,
                        batch_experiment_index=row,
                    )

                    if "mean_count_rate" in data:
                        mcr = ExperimentSingleValueData.objects.create(experiment=character,
                                                                       name="mean_count_rate",
                                                                       type=DataTypesChoices.INT,
                                                                       value=int(data["mean_count_rate"]))

                    characterizations.append(character)
                    step = "Nanoparticle"
                    nanoparticle = Nanoparticle.objects.create(
                        name=data['sample_name'],
                        z_average=data['z_average'],
                        mean_diameter_by_volume=data['mean_diameter_by_volume'],
                        mean_diameter_by_number=data['mean_diameter_by_number'],
                        mean_diameter_by_intensity=data['mean_diameter_by_intensity'],
                        pdi=data['pdi'],
                        )
                    nanoparticle.characterizations.add(character)

                    nanoparticles.append(nanoparticle)

                batch_experiment.run_date = earliest
                batch_experiment.save()

                request.session['batch_particle_creation_1'] = {'batch_experiment': batch_experiment.pk,
                                                                "nanoparticles": [np.pk for np in nanoparticles],
                                                                "characterizations": [c.pk for c in characterizations],
                                                                }
                np_formset = []
                for npi in nanoparticles:
                    np_formset.append(NanoparticleForm(instance=npi, readonly="__all__", remove=[
                        "preparation_method",
                        "materials",
                        "additives",
                        "characterizations",
                    ]))

                np_char_formset = []
                for character in characterizations:
                    np_char_formset.append(
                        NanoparticleCharacterizationForm(instance=character, readonly="__all__", hide=[
                            "preparation_method",
                            "materials",
                            "additives",
                            "characterizations",
                        ]))

                context = {'np_formset': np_formset,
                           'np_char_formset': np_char_formset,
                           'df': mark_safe(df.to_html),
                           'np_character_batch': batch_experiment
                           }
                return render(request, "experiments_nanoparticle/batch_create_2.html", context)
            except Exception as e:
                batch_experiment.delete()
                for np in nanoparticles:
                    np.delete()
                characterization_form.add_error(None, step + ": " + str(e))

        characterization_form.helper.add_input(Hidden("next_step", 1))
        context = {'characterization_form': characterization_form}
        return render(request, "experiments_nanoparticle/batch_create_1.html", context)

    def post2(self, request):
        chem_db_user = request.user.chemdbuser
        batch_characterization = NanoparticleCharacterization.objects.get(owner=chem_db_user,
                                                                          pk=request.POST["np_character_batch"])
        np_session_data = request.session['batch_particle_creation_1']
        assert batch_characterization.pk == np_session_data['batch_experiment']

        return redirect("chemicaldb:experiments:experiments_nanoparticle:view_characterization",
                        pk=batch_characterization.pk)

    def post3(self, request):
        chem_db_user = request.user.chemdbuser
        batch_characterization = NanoparticleCharacterization.objects.get(owner=chem_db_user,
                                                                          pk=request.POST["np_character_batch"])
        np_session_data = request.session['batch_particle_creation_1']
        assert batch_characterization.pk == np_session_data['batch_experiment']

        experiment_np_sets = [subex.nanoparticlecharacterization.characterized_nanoparticle.all() for subex in
                              batch_characterization.sub_experiments.all()]
        experiment_np = []
        for np_set in experiment_np_sets:
            experiment_np.extend(np_set)
        assert set([np.pk for np in experiment_np]) == set(np_session_data["nanoparticles"])
        for np in experiment_np:
            np.delete()
        batch_characterization.delete()
        return self.get(request)


def batch_characterizations(request):
    batches = NanoparticleCharacterization.objects.annotate(sub_experiment_count=Count('sub_experiments')).filter(
        sub_experiment_count__gte=1)
    context = {'batches': batches}
    return render(request, "experiments_nanoparticle/list_batch_characterizations.html", context)


def view_characterizations(request, pk):
    characterization = NanoparticleCharacterization.objects.get(pk=pk)
    context = {'characterization': characterization}
    return render(request, "experiments_nanoparticle/view_characterization.html", context)


def view_particle(request, pk):
    np = Nanoparticle.objects.get(pk=pk)
    context = {'np': np}
    return render(request, "experiments_nanoparticle/view_particle.html", context)


def batch_edit_characterization(request):
    chem_db_user = request.user.chemdbuser
    pks = [int(pk) for pk in request.GET.get("pk").split(",")]
    characterizations = NanoparticleCharacterization.objects.filter(pk__in=pks, owner=chem_db_user)
    context = {"characterizations": characterizations}
    return render(request, "experiments_nanoparticle/batch_edit_characterization.html", context)


class BatchEditParticles(View):
    def get(self, request):
        chem_db_user = request.user.chemdbuser
        pks = [int(pk) for pk in request.GET.get("pk").split(",")]
        particles = Nanoparticle.objects.filter(pk__in=pks, user=chem_db_user)
        context = {"particles": particles,
                   "formulations":NanoparticlePreparationMethod.objects.all(),
                   "concentration_units":ConcentrationUnits.choices,
                   }
        return render(request, "experiments_nanoparticle/batch_edit_particles.html", context)

    def post(self, request):
        post = dict(request.POST)
        print(post)

        chem_db_user = request.user.chemdbuser
        pks = [int(pk) for pk in request.GET.get("pk").split(",")]
        particles = Nanoparticle.objects.filter(pk__in=pks, user=chem_db_user)

        pk_particle = [int(pk) for pk in request.POST.getlist('pk_particle', [])]
        preparation_method_pk = [None if pk == "" else int(pk) for pk in request.POST.getlist('preparation_method_pk', [])]
        length=len(particles)
        equal_particle_length = all(len(l) == length for l in [pk_particle,preparation_method_pk]) and all(p.pk in pk_particle for p in particles)

        material_np=[int(pk) for pk in request.POST.getlist('material_np', [])]
        material_pk=[int(pk) for pk in request.POST.getlist('material_pk', [])]
        material_relative_content=[float(c) for c in request.POST.getlist('material_relative_content', [])]
        mat_len=len(material_np)
        material_check = all(np in pk_particle for np in material_np) and all(len(l) == mat_len for l in [material_pk,material_relative_content])

        additive_np=[int(pk) for pk in request.POST.getlist('additive_np', [])]
        additive_pk=[int(pk) for pk in request.POST.getlist('additive_pk', [])]
        additive_concentration=[float(c) for c in request.POST.getlist('additive_concentration', [])]
        additive_unit=request.POST.getlist('additive_unit', [])
        add_len=len(additive_np)
        additive_check = all(np in pk_particle for np in additive_np) and all(len(l) == add_len for l in [additive_pk,additive_concentration,additive_unit])

        print(equal_particle_length , material_check , additive_check)
        if equal_particle_length and material_check and additive_check:
            prep_method_by_pk={}
            particles_by_pk ={}
            material_by_pk={}

            for p in particles:
                particles_by_pk[p.pk]=p

            for pk in preparation_method_pk:
                if pk not in prep_method_by_pk:
                    if pk is None:
                        method = None
                    else:
                        method = NanoparticlePreparationMethod.objects.get(pk=pk)
                    prep_method_by_pk[pk] = method

            for p in material_pk:
                material_by_pk[p]=Substance.objects.get(pk=p)
            for p in additive_pk:
                material_by_pk[p]=Substance.objects.get(pk=p)

            for i in range(len(pk_particle)):
                p=particles_by_pk[pk_particle[i]]
                method = prep_method_by_pk[preparation_method_pk[i]]
                update=False
                if p.preparation_method != method:
                    print(p.preparation_method,method)
                    p.preparation_method = method
                    update = True

                if update:
                    p.save()

            for i in range(add_len):
                np=particles_by_pk[additive_np[i]]
                mat=material_by_pk[additive_pk[i]]
                conc=additive_concentration[i]
                unit = additive_unit[i]
                additive, created = Additives.objects.get_or_create(nanoparticle=np, material=mat,concentration=conc,concentration_unit=unit)

            for i in range(mat_len):
                np=particles_by_pk[material_np[i]]
                mat=material_by_pk[material_pk[i]]
                cont=material_relative_content[i]
                material, created = Materials.objects.get_or_create(nanoparticle=np, material=mat,relative_content=cont)

        return self.get(request)

class AddPreparationMethod(View):
    def get(self, request):
        chem_db_user = request.user.chemdbuser
        context = {"preparation_form": NanoparticlePreparationMethodForm(chem_db_user=chem_db_user,changeable=True)}
        return render(request, "experiments_nanoparticle/add_prepation_method.html", context)

    def post(self, request):
        chem_db_user = request.user.chemdbuser
        prep_form = NanoparticlePreparationMethodForm(chem_db_user=chem_db_user, data=request.POST,changeable=True)
        if prep_form.is_valid():
            prep_form.save()


class ViewPreparationMethod(View):
    def get(self, request,pk):
        chem_db_user = request.user.chemdbuser
        instance=NanoparticlePreparationMethod.objects.get(pk=pk)
        context = {"preparation_form": NanoparticlePreparationMethodForm(instance=instance,changeable=instance.owner==chem_db_user)}
        return render(request, "experiments_nanoparticle/add_prepation_method.html", context)

    def post(self, request,pk):
        chem_db_user = request.user.chemdbuser
        instance=NanoparticlePreparationMethod.objects.get(pk=pk)
        assert chem_db_user == instance.owner
        assert chem_db_user.pk == int(request.POST.get('owner'))
        prep_form = NanoparticlePreparationMethodForm(instance=instance,data=request.POST,changeable=instance.owner==chem_db_user)
        if prep_form.is_valid():
            prep_form.save()
        else:
            print(prep_form.errors)

        context = {"preparation_form": prep_form}
        return render(request, "experiments_nanoparticle/add_prepation_method.html", context)