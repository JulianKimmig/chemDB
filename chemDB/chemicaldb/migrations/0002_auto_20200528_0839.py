# Generated by Django 3.0.6 on 2020-05-28 08:39

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    initial = True

    dependencies = [
        ('sources', '0001_initial'),
        ('chemicaldb', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='meltingpoint',
            name='source',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.SET_NULL, to='sources.Source'),
        ),
        migrations.AddField(
            model_name='meltingpoint',
            name='substance',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='chemicaldb.Substance'),
        ),
        migrations.AddField(
            model_name='chemdbuser',
            name='institute',
            field=models.ForeignKey(blank=True, default=None, null=True, on_delete=django.db.models.deletion.SET_NULL, to='chemicaldb.ChemdbInstitute'),
        ),
        migrations.AddField(
            model_name='boilingpoint',
            name='source',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.SET_NULL, to='sources.Source'),
        ),
        migrations.AddField(
            model_name='boilingpoint',
            name='substance',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='chemicaldb.Substance'),
        ),
        migrations.AddField(
            model_name='simplesubstance',
            name='structure',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='chemicaldb.Structure'),
        ),
        migrations.AddField(
            model_name='mixedsubstance',
            name='structures',
            field=models.ManyToManyField(to='chemicaldb.Structure'),
        ),
        migrations.AlterUniqueTogether(
            name='chemdbuser',
            unique_together={('short', 'institute')},
        ),
    ]