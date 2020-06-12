# Generated by Django 3.0.6 on 2020-06-09 13:01

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('chemicaldb', '0010_auto_20200609_1501'),
        ('chemicaldb_polymers', '0003_auto_20200602_0610'),
    ]

    operations = [
        migrations.CreateModel(
            name='EndGroup',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='StartGroup',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.AddField(
            model_name='polymer',
            name='end_group',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='as_polymer_endgroup', to='chemicaldb.Structure'),
        ),
        migrations.AddField(
            model_name='polymer',
            name='start_group',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='as_polymer_startgroup', to='chemicaldb.Structure'),
        ),
        migrations.AlterField(
            model_name='monomers',
            name='repeating_units',
            field=models.CharField(choices=[('GRADIENT', 'gradient'), ('STATISTICAL', 'statistical'), ('BLOCK', 'block'), ('MIXED', 'mixed')], default='STATISTICAL', max_length=16),
        ),
        migrations.AlterField(
            model_name='polymer',
            name='monomer_distribution',
            field=models.CharField(choices=[('GRADIENT', 'gradient'), ('STATISTICAL', 'statistical'), ('BLOCK', 'block'), ('MIXED', 'mixed')], default='STATISTICAL', max_length=12),
        ),
    ]