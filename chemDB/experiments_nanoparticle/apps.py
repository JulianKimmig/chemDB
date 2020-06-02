from django.apps import AppConfig


class ExperimentsNanoparticleConfig(AppConfig):
    name = 'experiments_nanoparticle'

    def ready(self):
        from chemicaldb.api_views import register_search_model
        from .models import Nanoparticle, NanoparticlePreparationMethod
        register_search_model(Nanoparticle,qs=["name__icontains","code__icontains"],fields=["name","code"])

        register_search_model(NanoparticlePreparationMethod,qs=["name__icontains"],fields=["name"])

