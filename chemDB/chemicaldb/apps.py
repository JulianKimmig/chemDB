from django.apps import AppConfig


class ChemicaldbConfig(AppConfig):
    name = 'chemicaldb'

    def ready(self):
        from chemicaldb.api_views import register_search_model
        from chemicaldb.sub_models import Substance, Structure, StructureName
        register_search_model(Substance,qs=["name__icontains","code__icontains"],fields=["name","code"],name="Substance")
        register_search_model(Structure,qs=["code__icontains"],fields=["code"],name="Structure")