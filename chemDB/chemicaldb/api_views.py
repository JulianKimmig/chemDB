from rest_framework import routers, serializers, viewsets

from chemicaldb.models import Structure


class StructureSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = Structure
        fields = ['iso_smiles', 'standard_inchi', 'inchi_key', 'molar_mass', 'valid', 'cas_number', 'chemspider_id']


class StructureViewSet(viewsets.ModelViewSet):
    queryset = Structure.objects.all()
    serializer_class = StructureSerializer


router = routers.DefaultRouter()
router.register(r'structure', StructureViewSet)
print(router)


class SearchModel:
    all_searches = []
    models_by_weight = {}

    def __init__(self, weight=10):
        SearchModel.all_searches.append(self)
        if weight not in SearchModel.models_by_weight:
            SearchModel.models_by_weight[weight] = []
        SearchModel.models_by_weight[weight].append(self)


class SearchViewSet(viewsets.ListAPIView):
    def list(self, request):
        queryset = list(itertools.chain(Tweet.objects.all(), Article.objects.all()))
        serializer = TimelineSerializer(queryset, many=True)
        return Response(serializer.data)
