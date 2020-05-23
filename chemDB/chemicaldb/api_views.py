from rest_framework import routers, serializers, viewsets
from rest_framework.response import Response
from rest_framework.views import APIView

from chemicaldb.models import Structure


class StructureSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = Structure
        fields = ['iso_smiles', 'standard_inchi', 'inchi_key', 'molar_mass', 'valid', 'cas_number', 'chemspider_id']


class StructureViewSet(viewsets.ModelViewSet):
    queryset = Structure.objects.all()
    serializer_class = StructureSerializer






class SearchModel:
    all_searches = []
    models_by_weight = {}

    def __init__(self, weight=10):
        SearchModel.all_searches.append(self)
        if weight not in SearchModel.models_by_weight:
            SearchModel.models_by_weight[weight] = []
        SearchModel.models_by_weight[weight].append(self)


class SearchViewSet(APIView):

    serializer_class = StructureSerializer

    def get(self, request, format=None, **kwargs):
        structures = Structure.objects.all()
        structures_serializer = StructureSerializer(structures)

        return Response({
            'structures': structures_serializer.data
        })
    #def list(self, request):
    #    queryset = self.get_queryset()
    #    serializer = self.get_serializer()(queryset, many=True)
    #    return Response(serializer.data)


router = routers.DefaultRouter()
router.register(r'structure', StructureViewSet)

