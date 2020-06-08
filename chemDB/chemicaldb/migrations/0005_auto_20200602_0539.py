# Generated by Django 3.0.6 on 2020-06-02 03:39

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('contenttypes', '0002_remove_content_type_name'),
        ('chemicaldb', '0004_auto_20200601_0909'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='mixedsubstance',
            options={'base_manager_name': 'objects'},
        ),
        migrations.AlterModelOptions(
            name='simplesubstance',
            options={'base_manager_name': 'objects'},
        ),
        migrations.AlterModelOptions(
            name='substance',
            options={'base_manager_name': 'objects'},
        ),
        migrations.AddField(
            model_name='substance',
            name='polymorphic_ctype',
            field=models.ForeignKey(editable=False, null=True, on_delete=django.db.models.deletion.CASCADE, related_name='polymorphic_chemicaldb.substance_set+', to='contenttypes.ContentType'),
        ),
        migrations.AlterField(
            model_name='substance',
            name='substance_class',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, to='chemicaldb.SubstanceClass'),
        ),
        migrations.AlterField(
            model_name='substanceclass',
            name='class_name',
            field=models.CharField(max_length=100, unique=True),
        ),
    ]