from django.db import models

# Create your models here.
from rdkit.Chem import rdDepictor, MolFromSmiles
from rdkit.Chem.Draw import rdMolDraw2D

from chemicaldb.models import Structure, Substance, mark_safe


class PolymerStructureType(models.TextChoices):
    LINEAR = 'LINEAR', "linear"
    BRANCHED = 'BRANCHED', "branched"
    NETWORK = 'NETWORK', "network"



class PolymerStructure(Structure):
    big_smiles = models.CharField(max_length=200, null=True, blank=True)

    def __init__(self, *args, **kwargs):
        self._meta.get_field('iso_smiles').default = False
        super().__init__(*args, **kwargs)

    def structure_image(self):
        mol = MolFromSmiles(self.smiles*4)
        #mc = Chem.Mol(mol.ToBinary())
        if not mol.GetNumConformers():
            rdDepictor.Compute2DCoords(mol)
        drawer = rdMolDraw2D.MolDraw2DSVG(200,200)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        return mark_safe(svg)


class MonomerDistribution(models.TextChoices):
    GRADIENT = 'GRADIENT', "gradient"
    STATISTICAL = 'STATISTICAL', "statistical"
    BLOCK = 'BLOCK', "block"
    MIXED = 'MIXED', "mixed"

MonomerDistribution.max_length = 12

class Polymer(Substance):
    polymer_structure = models.ForeignKey(PolymerStructure, on_delete=models.CASCADE,related_name="polymer_structure")
    polymer_structure_type = models.CharField(max_length=16, choices=PolymerStructureType.choices,
                                              default=PolymerStructureType.LINEAR)
    monomers = models.ManyToManyField(Structure, through='Monomers')
    monomer_distribution = models.CharField(max_length=MonomerDistribution.max_length, choices=MonomerDistribution.choices,
                                            default=MonomerDistribution.STATISTICAL)
    distribution_parameter = models.FloatField(default=1)
    starting_group=models.ForeignKey(Structure, on_delete=models.SET_NULL,null=True)
    end_group=models.ForeignKey(Structure, on_delete=models.SET_NULL,null=True)




class Monomers(models.Model):
    polymer = models.ForeignKey(Polymer, on_delete=models.CASCADE,related_name="p")
    monomer = models.ForeignKey(Structure, on_delete=models.CASCADE,)
    repeating_units = models.CharField(max_length=16, choices=MonomerDistribution.choices,
                                       default=MonomerDistribution.STATISTICAL)
