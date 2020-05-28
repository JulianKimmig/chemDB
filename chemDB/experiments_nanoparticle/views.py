from datetime import datetime

from crispy_forms.layout import Hidden
from django.core import serializers
from django.db.models import Count
from django.forms import modelformset_factory
from django.shortcuts import render, redirect

# Create your views here.
from django.utils.safestring import mark_safe
from django.views import View

from experiments.models import ExperimentRawData, ExperimentData, DataTypesChoices, ExperimentSingleValueData
from experiments_nanoparticle.characterization_tools import NanoparticleCharacterizationTool
from experiments_nanoparticle.models import NanoparticleBatchCharacterizationForm, Nanoparticle, NanoparticleForm, \
    NanoparticleCharacterization, NanoparticleCharacterizationForm


def main_view(request):
    return render(request, "experiments_nanoparticle/main.html")


class batch_particle_creation(View):
    def get(self,request):
        chem_db_user=request.user.chemdbuser
        characterization_form=NanoparticleBatchCharacterizationForm(chem_db_user=chem_db_user)
        characterization_form.helper.add_input(Hidden("next_step", 1))
        context = {'characterization_form':characterization_form}
        return render(request, "experiments_nanoparticle/batch_create_1.html",context)

    def post(self,request):
        step = int(request.POST.get("next_step","1"))
        if step == 1:
            return self.post1(request)
        if step == 2:
            return self.post2(request)
        if step == 3:
            return self.post3(request)
        return self.get(request)

    def post1(self, request):
        chem_db_user=request.user.chemdbuser
        characterization_form=NanoparticleBatchCharacterizationForm(chem_db_user=chem_db_user,data=request.POST,files=request.FILES,)
        if characterization_form.is_valid():
            form_data=characterization_form.cleaned_data

            tool = getattr(NanoparticleCharacterizationTool,form_data['tool'])
            df = tool.read_batch_data(data=form_data['exported_data'].read(),name=form_data['name'],short_name=form_data['short_name'])

            batch_experiment = characterization_form.save(commit=True)
            batch_experiment.owner=chem_db_user

            ExperimentRawData.objects.create(name="exported_data",raw_data=form_data['exported_data'],experiment=batch_experiment)
            ExperimentRawData.objects.create(name="raw_data",raw_data=form_data['raw_data'],experiment=batch_experiment)

            earliest = df['run_date'].min()
            nanoparticles=[]
            characterizations=[]
            step=""
            try:
                step="NanoparticleCharacterization"
                for row,data in df.iterrows():
                    #   if earliest is None:
                 #       earliest = data['run_date']
                 #   if data['run_date'] < earliest:
                 #       earliest = data['run_date']

                    character = NanoparticleCharacterization.objects.create(
                        name="{} in {}".format(data['sample_name'], batch_experiment.short_name),
                        short_name="{}_np-characterization".format(data['sample_name']),
                        owner = chem_db_user,
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
                    step="Nanoparticle"
                    nanoparticle = Nanoparticle.objects.create(
                        name=data['sample_name'],
                        code=chem_db_user.prefix_string("{}_{}".format(batch_experiment.short_name, row)),
                        z_average = data['z_average'],
                        mean_diameter_by_volume = data['mean_diameter_by_volume'],
                        mean_diameter_by_number = data['mean_diameter_by_number'],
                        mean_diameter_by_intensity = data['mean_diameter_by_intensity'],
                        pdi = data['pdi'],
                    )
                    nanoparticle.characterizations.add(character)

                    nanoparticles.append(nanoparticle)

                batch_experiment.run_date = earliest
                batch_experiment.save()

                request.session['batch_particle_creation_1'] = {'batch_experiment':batch_experiment.pk,
                                                                "nanoparticles":[np.pk for np in nanoparticles],
                                                                "characterizations":[c.pk for c in characterizations],
                                                                }
                np_formset=[]
                for npi in nanoparticles:
                    np_formset.append(NanoparticleForm(instance=npi,readonly="__all__",remove=[
                        "preparation_method",
                        "materials",
                        "additives",
                        "characterizations",
                    ]))

                np_char_formset=[]
                for character in characterizations:
                    np_char_formset.append(NanoparticleCharacterizationForm(instance=character,readonly="__all__",hide=[
                        "preparation_method",
                        "materials",
                        "additives",
                        "characterizations",
                    ]))

                context={'np_formset':np_formset,
                         'np_char_formset':np_char_formset,
                         'df':mark_safe(df.to_html),
                         'np_character_batch':batch_experiment
                         }
                return render(request, "experiments_nanoparticle/batch_create_2.html",context)
            except Exception as e:
                batch_experiment.delete()
                for np in nanoparticles:
                    np.delete()
                characterization_form.add_error(None,step+": "+str(e))

        characterization_form.helper.add_input(Hidden("next_step", 1))
        context = {'characterization_form':characterization_form}
        return render(request, "experiments_nanoparticle/batch_create_1.html",context)

    def post2(self, request):
        chem_db_user=request.user.chemdbuser
        batch_characterization = NanoparticleCharacterization.objects.get(owner=chem_db_user,pk=request.POST["np_character_batch"])
        np_session_data = request.session['batch_particle_creation_1']
        assert batch_characterization.pk == np_session_data['batch_experiment']

        return redirect("chemicaldb:experiments:experiments_nanoparticle:view_characterization",pk=batch_characterization.pk)

    def post3(self, request):
        chem_db_user=request.user.chemdbuser
        batch_characterization = NanoparticleCharacterization.objects.get(owner=chem_db_user,pk=request.POST["np_character_batch"])
        np_session_data = request.session['batch_particle_creation_1']
        assert batch_characterization.pk == np_session_data['batch_experiment']

        experiment_np_sets = [subex.nanoparticlecharacterization.characterized_nanoparticle.all() for subex in batch_characterization.sub_experiments.all()]
        experiment_np=[]
        for np_set in experiment_np_sets:
            experiment_np.extend(np_set)
        assert set([np.pk for np in experiment_np]) == set(np_session_data["nanoparticles"])
        for np in experiment_np:
            np.delete()
        batch_characterization.delete()
        return self.get(request)

def batch_characterizations(request):
   batches =  NanoparticleCharacterization.objects.annotate(sub_experiment_count=Count('sub_experiments')).filter(sub_experiment_count__gte=1)
   context={'batches':batches}
   return render(request, "experiments_nanoparticle/list_batch_characterizations.html",context)

def view_characterizations(request,pk):
    characterization = NanoparticleCharacterization.objects.get(pk=pk)
    context={'characterization':characterization}
    return render(request, "experiments_nanoparticle/view_characterization.html",context)

def view_particle(request,pk):
    np = Nanoparticle.objects.get(pk=pk)
    context={'np':np}
    return render(request, "experiments_nanoparticle/view_particle.html",context)


def batch_edit_characterization(request):
    chem_db_user=request.user.chemdbuser
    pks = [int(pk) for pk in request.GET.get("pk").split(",")]
    characterizations = NanoparticleCharacterization.objects.filter(pk__in=pks,owner=chem_db_user)
    context={"characterizations":characterizations}
    return render(request, "experiments_nanoparticle/batch_edit_characterization.html",context)