# Generated by Django 3.0.6 on 2020-05-30 03:40

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('chemicaldb', '0002_auto_20200528_0839'),
    ]

    operations = [
        migrations.AddField(
            model_name='chemdbinstitute',
            name='parent_institute',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='child_institute', to='chemicaldb.ChemdbInstitute'),
        ),
        migrations.AlterField(
            model_name='substance',
            name='code',
            field=models.CharField(max_length=24, null=True, unique=True),
        ),
    ]
