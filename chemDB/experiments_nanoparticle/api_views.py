from rest_framework import routers, serializers, viewsets

from experiments_nanoparticle.models import Nanoparticle


class NanoparticleSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = Nanoparticle
        fields = ['preparation_method', 'materials', 'additives']


class NanoparticleViewSet(viewsets.ModelViewSet):
    queryset = Nanoparticle.objects.all()
    serializer_class = NanoparticleSerializer

router = routers.DefaultRouter()
router.register(r'nanopoarticle', NanoparticleViewSet)