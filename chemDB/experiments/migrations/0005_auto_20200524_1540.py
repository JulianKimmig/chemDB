# Generated by Django 3.0.6 on 2020-05-24 13:40

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('experiments_nanoparticle', '0003_delete_nanoparticlecharacterizationproperty'),
        ('experiments', '0004_auto_20200524_1531'),
    ]

    operations = [
        migrations.CreateModel(
            name='ExperimentSingleValueData',
            fields=[
                ('experimentdata_ptr', models.OneToOneField(auto_created=True, on_delete=django.db.models.deletion.CASCADE, parent_link=True, primary_key=True, serialize=False, to='experiments.ExperimentData')),
                ('value', models.CharField(max_length=32)),
            ],
            bases=('experiments.experimentdata',),
        ),
        migrations.DeleteModel(
            name='ExperimentDataPointX',
        ),
    ]
