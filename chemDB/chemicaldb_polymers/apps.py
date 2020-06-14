from django.apps import AppConfig
from django.urls import reverse


class ChemicaldbPolymersConfig(AppConfig):
    name = 'chemicaldb_polymers'

    def ready(self):
        from chemicaldb.views import register_structure_plugin
        from django.templatetags.static import static
        register_structure_plugin(name=self.name,title="Polymers",
                                  entry_link=reverse('chemicaldb:chemicaldb_polymers:structures_main'),
                                  image=static('chemicaldb_polymers/PMAA_polymer.svg'))