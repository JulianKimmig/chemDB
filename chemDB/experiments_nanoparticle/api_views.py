from rest_framework import routers, serializers, viewsets, generics
from rest_framework.views import APIView

from rest_framework.response import Response

from experiments_nanoparticle.models import Nanoparticle


class StructureSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = Nanoparticle
        fields = ['iso_smiles', 'standard_inchi', 'inchi_key', 'molar_mass', 'valid', 'cas_number', 'chemspider_id']


class StructureViewSet(viewsets.ModelViewSet):
    queryset = Nanoparticle.objects.all()
    serializer_class = StructureSerializer






class SearchModel:
    all_searches = []
    models_by_weight = {}

    def __init__(self, weight=10):
        SearchModel.all_searches.append(self)
        if weight not in SearchModel.models_by_weight:
            SearchModel.models_by_weight[weight] = []
        SearchModel.models_by_weight[weight].append(self)

router = routers.DefaultRouter()
router.register(r'structure', StructureViewSet)