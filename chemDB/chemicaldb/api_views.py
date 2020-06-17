from functools import partial

import rdkit
import rdkit.Chem
import rdkit.Chem.Descriptors
import rdkit.Chem.PandasTools
from django.apps import apps
from django.db.models import Q
from django.db.models.base import ModelBase
from django.http import JsonResponse, HttpResponse
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rest_framework import routers, serializers, viewsets, permissions
from rest_framework.decorators import permission_classes
from rest_framework.response import Response
from rest_framework.views import APIView

from chemicaldb.models import Structure, StructureName, ValidationError
from chemicaldb.services import structure_image, mol_to_image_mol, structure_check


class StructureSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = Structure
        fields = ['smiles', 'standard_inchi', 'inchi_key', 'molar_mass', 'valid', 'cas_number', 'chemspider_id']


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

    return model.objects.filter(obj_query)


def register_search_model(model, qs=[], fields=[], name=None,return_map=None):
    model_class = model.__module__ + "." + model.__name__
    qs.append(("pk", int))
    qs = list(set(qs))
    if name is None:
        name = model.__name__
    SEARCHABLE_MODELS[name] = (
        (model, partial(_search, model=model, queries=qs), fields, model_class,return_map)
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
        for model_name, (model, model_search, model_fields, model_class,return_map) in SEARCHABLE_MODELS.items():
            if model_name in models:
                model_results = model_search(query, search_fields=search_fields)
                if issubclass(model,Structure):
                    named = StructureName.objects.filter(name__icontains=query).values_list('structure', flat=True)
                    if model == Structure:
                        model_results = model_results | model.objects.filter(pk__in=named)
                    else:
                        model_results = model_results | model.objects.filter(structure_ptr_id__in=named)

                model_results=model_results.distinct().all()

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

def smiles_checker(request):
    try:
        smiles = request.GET.get("smiles","").strip()
        assert smiles is not None, "no smiles found"
        image_format = request.GET.get("image_format")
        structure_model = request.GET.get("model","chemicaldb.Structure").split(".")
        model = apps.get_model(".".join(structure_model[:-1]), structure_model[-1])
        validate = request.GET.get("validate",True)
        temp_instance:Structure=model(smiles=smiles)
        mol=temp_instance.create_mol()
        return JsonResponse(structure_check(temp_instance,mol=mol,image_format=image_format,validate=validate))
    except ValidationError as e:
        return JsonResponse({"success":False,
                                      "reason":str(e)})
def smiles_to_image(request):
    smiles = request.GET.get("smiles")
    format = request.GET.get("type")
    validate = request.GET.get("validate")

    mol = rdkit.Chem.MolFromSmiles(smiles)

    if not mol.GetNumConformers():
        rdDepictor.Compute2DCoords(mol)

    drawer = rdMolDraw2D.MolDraw2DSVG(200,200)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    return HttpResponse(svg)