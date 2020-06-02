from crispy_forms.helper import FormHelper
from crispy_forms.layout import Submit
from django import forms
from django.db import models

# Create your models here.
from django.templatetags.static import static
from django.utils.safestring import mark_safe

from chemicaldb.sub_models import Substance, Structure
from experiments.models import GeneraExperiment, ConcentrationUnits, Experiment
from experiments_nanoparticle.characterization_tools import NanoparticleCharacterizationTool


class NanoparticleCharacterization(Experiment):
    tool = models.CharField(max_length=NanoparticleCharacterizationTool.max_length,
                            choices=NanoparticleCharacterizationTool.choices,
                            default=str(NanoparticleCharacterizationTool.MZetaNano1))

    def deep_characterized_nanoparticle(self,ignore=[]):
        particle=[p for p in self.characterized_nanoparticle.all()]
        ignore.append(self.pk)
        for sub_experiment in self.sub_experiments.all():
            if sub_experiment.pk not in ignore:
                particle.extend(sub_experiment.nanoparticlecharacterization.deep_characterized_nanoparticle(ignore=ignore))
        return particle

class NanoparticleCharacterizationForm(forms.ModelForm):
    class Meta:
        model = NanoparticleCharacterization
        exclude = []

    def __init__(self, readonly=None, hide=None, *args, **kwargs):
        super().__init__(*args, **kwargs)

        if hide is None:
            hide = []
        if readonly is None:
            readonly = []

        for field_name in hide:
            if field_name in self.fields:
                self.fields[field_name].widget = forms.HiddenInput()

        if readonly == "__all__":
            for field_name, field in self.fields.items():
                field.widget.attrs['readonly'] = True
                field.widget.attrs['disabled'] = True
        else:
            for field_name in readonly:
                if field_name in self.fields:
                    self.fields[field_name].widget.attrs['readonly'] = True

class NanoparticleBatchCharacterizationForm(forms.ModelForm):
    class Meta:
        model = NanoparticleCharacterization
        exclude = ["run_date", "batch_experiment", "batch_experiment_index"]

    raw_data = forms.FileField(help_text="the experimented raw data, will not be interpreted (until know)")
    exported_data = forms.FileField(widget=forms.FileInput(attrs={'accept': '.txt,.csv'}),
                                    help_text=mark_safe(
                                        "the experimented raw data, will not be interpreted (until know).</br>"
                                        "For '{}' get the export template <a href='{}' download>here</a>.".format(
                                            NanoparticleCharacterizationTool.MZetaNano1.long_name,
                                            static(
                                                'experiments_nanoparticle/characterization_templates/chemdb_zetasizer_np_export_1.zip')
                                        )),
                                    required=True,
                                    )

    def __init__(self,chem_db_user, *args, **kwargs):

        super().__init__(*args, **kwargs)
        self.helper = FormHelper()
        self.helper.form_method = 'post'
        self.helper.add_input(Submit('submit', 'Continue'))
        self.fields["owner"].initial = chem_db_user.pk
        self.fields["owner"].widget = forms.HiddenInput()
        self.fields["short_name"].widget.attrs['placeholder']=chem_db_user.prefix_string("np_data")
        self.fields["short_name"].help_text ="if does not start with '{}' the prefix will be added automatically.So you don't have to enter it".format(chem_db_user.get_prefix())
#        initial={'owner':chem_db_user.pk,"short_name":}
    #   self.fields['raw_data'] = forms.FileInput()


class NanoparticlePreparationMethod(GeneraExperiment):
   pass

class NanoparticlePreparationMethodForm(forms.ModelForm):
    class Meta:
        model = NanoparticlePreparationMethod
        exclude = []

    def __init__(self,chem_db_user=None,changeable=False, *args, **kwargs):

        super().__init__(*args, **kwargs)
        self.helper = FormHelper()
        if changeable:
            self.helper.form_method = 'post'
            self.helper.add_input(Submit('submit', 'save'))
        if chem_db_user:
            self.fields["owner"].initial = chem_db_user.pk
        self.fields["owner"].widget = forms.HiddenInput()
#        initial={'owner':chem_db_user.pk,"short_name":}
#   self.fields['raw_data'] = forms.FileInput()

def _validate_prep_method(np):
    if np.preparation_method is None:
        return False
    return True


def _validate_material(np):
    if len(np.materials) <= 0:
        return False
    return True


NP_VALIDATION_PREPARATION = "np_preparation_method", _validate_prep_method
NP_VALIDATION_MATERIAL = "np_material", _validate_material



class Nanoparticle(Substance):
    preparation_method = models.ForeignKey(NanoparticlePreparationMethod, on_delete=models.SET_NULL, null=True)
    materials = models.ManyToManyField(Substance, through='Materials', related_name="as_nanopartice_material")
    additives = models.ManyToManyField(Substance, through='Additives', related_name="as_nanopartice_additive")
    characterizations = models.ManyToManyField(NanoparticleCharacterization,related_name="characterized_nanoparticle")

    z_average = models.PositiveIntegerField(null=True)
    mean_diameter_by_volume = models.PositiveIntegerField(null=True)
    mean_diameter_by_number = models.PositiveIntegerField(null=True)
    mean_diameter_by_intensity = models.PositiveIntegerField(null=True)
    pdi = models.FloatField(null=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.validity_checks.update({
            NP_VALIDATION_PREPARATION[0]: NP_VALIDATION_PREPARATION[1]
        })

Nanoparticle.register_substance_class()


class NanoparticleForm(forms.ModelForm):
    class Meta:
        model = Nanoparticle
        exclude = ["user"]

    def __init__(self, readonly=None, hide=None,remove=None, *args, **kwargs):
        super().__init__(*args, **kwargs)

        if hide is None:
            hide = []
        if readonly is None:
            readonly = []
        if remove is None:
            remove = []

        for field_name in remove:
            if field_name in self.fields:
                del self.fields[field_name]

        for field_name in hide:
            if field_name in self.fields:
                self.fields[field_name].widget = forms.HiddenInput()

        if readonly == "__all__":
            for field_name, field in self.fields.items():
                field.widget.attrs['readonly'] = True
                field.widget.attrs['disabled'] = True
        else:
            for field_name in readonly:
                if field_name in self.fields:
                    self.fields[field_name].widget.attrs['readonly'] = True


class Materials(models.Model):
    nanoparticle = models.ForeignKey(Nanoparticle, on_delete=models.CASCADE, related_name="np_material_particle")
    material = models.ForeignKey(Substance, on_delete=models.CASCADE, related_name="np_material")
    relative_content = models.FloatField()

    class Meta:
        unique_together=("nanoparticle","material")

class Additives(models.Model):
    nanoparticle = models.ForeignKey(Nanoparticle, on_delete=models.CASCADE, related_name="np_additives_particle")
    material = models.ForeignKey(Substance, on_delete=models.CASCADE, related_name="np_additive")
    concentration = models.FloatField()
    concentration_unit = models.CharField(max_length=ConcentrationUnits.max_length, choices=ConcentrationUnits.choices,
                                          default=ConcentrationUnits.MOL_LITER)
