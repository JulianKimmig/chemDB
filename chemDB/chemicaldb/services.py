import base64

import rdkit
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

from chemicaldb.models import ValidationError


def mol_to_image_mol(mol, image_format,size=300):
    if image_format == 'svg':
        drawer = rdMolDraw2D.MolDraw2DSVG(size, size)
        mode = 'w'
    elif image_format == 'png':
        drawer = rdMolDraw2D.MolDraw2DCairo(size, size)
        mode = 'wb'
    else:
        raise NotImplemented("unnknown format: {}".format(image_format))
    if not mol.GetNumConformers():
        rdDepictor.Compute2DCoords(mol)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    if image_format == "svg":
        return drawer.GetDrawingText()
    elif image_format == "png":
        return base64.b64encode(drawer.GetDrawingText()).decode('utf-8')
    raise NotImplemented("unnknown format: {}".format(image_format))


def structure_image(structure, image_format, *args, **kwargs):
    mol = structure.get_mol(*args, **kwargs)
    if mol is None:
        raise ValidationError("cannot get mol for molecule")
    return mol_to_image_mol(mol, image_format)


def structure_check(structure, mol, validate=True, image_format=None, raise_error=False):
    ans = {}

    structure.mol = mol
    if image_format:
        ans['img'] = mol_to_image_mol(mol, image_format)
    if validate:
        structure.validate()

    ans["success"] = True
    return ans