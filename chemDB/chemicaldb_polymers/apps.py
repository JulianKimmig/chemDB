from django.apps import AppConfig


class ChemicaldbPolymersConfig(AppConfig):
    name = 'chemicaldb_polymers'

    def ready(self):
        from django.contrib.contenttypes.models import ContentType
        from chemicaldb_polymers.models import PolymerStructure

        new_ct = ContentType.objects.get_for_model(PolymerStructure)
        PolymerStructure.objects.filter(polymorphic_ctype__isnull=True).update(polymorphic_ctype=new_ct)