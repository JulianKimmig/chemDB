from functools import partial

from django.db.models import Q
from django.http import JsonResponse
from rest_framework import routers, serializers, viewsets, permissions
from rest_framework.decorators import permission_classes
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


SEARCHABLE_MODELS = {}


def _search(query, model, queries, search_fields=None):
    _queries = []
    for q in queries:
        if isinstance(q, tuple):
            _queries.append(q)
        else:
            _queries.append((q, None))
    queries = _queries

    if search_fields is not None:
        _queries = []
        for q in queries:
            if q[0].split("__")[0] in search_fields:
                _queries.append(q)
        queries = _queries

    qs = []
    for q in queries:
        if q[1]:
            try:
                qs.append(Q(**{q[0]: q[1](query)}))
            except:
                pass
        else:
            qs.append(Q(**{q[0]: query}))

    if len(qs) == 0:
        return[]

    obj_query = qs.pop()
    for q in qs:
        obj_query |= q

    return model.objects.filter(obj_query).distinct().all()


def register_search_model(model, qs=[], fields=[], name=None):
    model_class = model.__module__ + "." + model.__name__
    qs.append(("pk", int))
    qs = list(set(qs))
    if name is None:
        name = model.__name__
    SEARCHABLE_MODELS[name] = (
        (model, partial(_search, model=model, queries=qs), fields, model_class)
    )


def search(request):
    query = request.GET.get("q","")
    models = request.GET.get("model")
    search_fields = request.GET.get("sf")

    if search_fields is not None:
        search_fields = search_fields.split(",")
    if models is not None:
        models = models.split(",")
    else:
        models = list(SEARCHABLE_MODELS.keys())
    results = []
    if query is not None:
        for model_name, (model, model_search, model_fields, model_class) in SEARCHABLE_MODELS.items():
            if model_name in models:
                model_results = model_search(query, search_fields=search_fields)

                model_data = [
                    {**{'__name__': model_name, '__class__': model_class, "__str__": str(res), "pk": res.pk},
                     **{v: getattr(res, v) for v in model_fields}
                     } for res in model_results
                ]
                results.extend(model_data)

    return JsonResponse({"results": results})


class SearchViewSet(APIView):
    serializer_class = StructureSerializer

    @permission_classes((permissions.AllowAny,))
    def get(self, request, format=None, **kwargs):
        structures = Structure.objects.all()
        structures_serializer = StructureSerializer(structures)

        return Response({
            'structures': structures_serializer.data
        })
    # def list(self, request):
    #    queryset = self.get_queryset()
    #    serializer = self.get_serializer()(queryset, many=True)
    #    return Response(serializer.data)


router = routers.DefaultRouter()
router.register(r'structure', StructureViewSet)
