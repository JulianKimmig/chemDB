from django.apps import AppConfig
from django.urls import reverse


class ChemicaldbPolymersConfig(AppConfig):
    name = 'chemicaldb_polymers'

    def ready(self):
        from chemicaldb.views import register_structure_plugin
        from django.templatetags.static import static
        from chemicaldb.api_views import register_search_model
        from chemicaldb_polymers.models import RepeatingUnit, StartEndGroup,PolymerStructure,Polymer

        register_structure_plugin(name=self.name,title="Polymers",
                                  entry_link=reverse('chemicaldb:chemicaldb_polymers:structures_main'),
                                  image=static('chemicaldb_polymers/PMAA_polymer.svg'))

        register_search_model(RepeatingUnit,qs=["code__icontains"],fields=["code","all_names"],name="RepeatingUnit")
        register_search_model(StartEndGroup,qs=["code__icontains"],fields=["code","all_names"],name="StartEndGroup")
        register_search_model(PolymerStructure,qs=["code__icontains"],fields=["code","all_names"],name="PolymerStructure")
        register_search_model(Polymer,qs=["code__icontains","name__icontains"],fields=["code","name"],name="Polymer")