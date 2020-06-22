from haystack import indexes

from chemicaldb.models import Structure


class StructureIndex(indexes.SearchIndex, indexes.Indexable):
    text = indexes.CharField(document=True, use_template=True)

    #author = indexes.CharField(model_attr='owner')

   # pub_date = indexes.DateTimeField(model_attr='pub_date')

    def get_model(self):
        return Structure

    def index_queryset(self, using=None):
        """Used when the entire index for model is updated."""
        return self.get_model().objects.all()