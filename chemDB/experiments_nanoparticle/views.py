from crispy_forms.layout import Hidden
from django.core import serializers
from django.forms import modelformset_factory
from django.shortcuts import render

# Create your views here.
from django.views import View

from experiments.models import ExperimentRawData
from experiments_nanoparticle.characterization_tools import NanoparticleCharacterizationTool
from experiments_nanoparticle.models import NanoparticleBatchCharacterizationForm, Nanoparticle, NanoparticleForm, \
    NanoparticleCharacterization


def main_view(request):
    return render(request, "experiments_nanoparticle/main.html")


class batch_particle_creation(View):
    def get(self,request):
        characterization_form=NanoparticleBatchCharacterizationForm(initial={'owner':request.user.chemdb_user.id})
        characterization_form.helper.add_input(Hidden("next_step", 1))
        context = {'characterization_form':characterization_form}
        return render(request, "experiments_nanoparticle/batch_create_1.html",context)

    def post(self,request):
        step = int(request.POST.get("next_step","1"))
        if step == 1:
            return self.post1(request)
        return self.get(request)

    def post1(self, request):
        characterization_form=NanoparticleBatchCharacterizationForm(request.POST,request.FILES,)
        if characterization_form.is_valid():
            batch_experiment = characterization_form.save(commit=True)
            batch_experiment.owner=request.user.chemdb_user

            form_data=characterization_form.cleaned_data
            tool = getattr(NanoparticleCharacterizationTool,form_data['tool'])
            _batch_experiment,_sub_experiments,_np,df = tool.read_batch_data(data=form_data['exported_data'].read(),name=form_data['name'],short_name=form_data['short_name'])

            ExperimentRawData.objects.create(name="exported_data",raw_data=form_data['exported_data'],experiment=batch_experiment)
            ExperimentRawData.objects.create(name="raw_data",raw_data=form_data['raw_data'],experiment=batch_experiment)

            earliest=None
            nanoparticles=[]
            characterizations = []
            step=""
            try:
                for row,data in df.iterrows():
                    if earliest is None:
                        earliest = data['run_date']
                    if data['run_date'] < earliest:
                        earliest = data['run_date']
                    step="NanoparticleCharacterization"
                    character = NanoparticleCharacterization.objects.create(
                        name="{} in {}".format(data['sample_name'], batch_experiment.short_name),
                        short_name="{}_np-characterization".format(data['sample_name']),
                        owner = request.user.chemdb_user,
                        tool=batch_experiment.tool,
                        run_date=data['run_date'],
                        batch_experiment=batch_experiment,
                        batch_experiment_index=row,
                    )
                    characterizations.append(character)
                    step="Nanoparticle"
                    nanoparticle = Nanoparticle.objects.create(
                        name=data['sample_name'],
                        code="{}_{}".format(batch_experiment.short_name,row),
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

                #request.session['batch_particle_creation_1'] = characterization_form.cleaned_data

                np_formset=[]
                for npi in nanoparticles:
                    np_formset.append(NanoparticleForm(instance=npi,readonly="__all__",hide=[
                        "preparation_method",
                        "materials",
                        "additives",
                        "characterizations",
                    ]))
                context={'np_formset':np_formset}
                return render(request, "experiments_nanoparticle/batch_create_2.html",context)
            except Exception as e:
                batch_experiment.delete()
                for np in nanoparticles:
                    np.delete()
                characterization_form.add_error(None,step+": "+str(e))

        characterization_form.helper.add_input(Hidden("next_step", 1))
        context = {'characterization_form':characterization_form}
        return render(request, "experiments_nanoparticle/batch_create_1.html",context)