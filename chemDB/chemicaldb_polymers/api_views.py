from django.http import JsonResponse

from chemicaldb.models import ValidationError
from chemicaldb.services import structure_check
from chemicaldb_polymers.models import StartEndGroup, PolymerStructure, RepeatingUnit


def structure_checker(request):
    try:
        chem_db_user = request.user.chemdbuser

        repeating_units_pks = request.GET.getlist("repeating_unit")

        if len(repeating_units_pks)==0:
            raise ValidationError("No repeating units")

        start_group = request.GET.get("start_group")

        temp_structure = PolymerStructure()

        if start_group:
            start_group = StartEndGroup.objects.get(pk=start_group)
            assert start_group.check_can_view(chem_db_user), "invalid start group"
            temp_structure.terminal_start_group = start_group

        end_group = request.GET.get("end_group")
        if end_group:
            end_group = StartEndGroup.objects.get(pk=end_group)
            assert end_group.check_can_view(chem_db_user), "invalid end group"
            temp_structure.terminal_end_group = end_group

        repeating_units = []
        for repeating_unit in repeating_units_pks:
            repeating_unit = RepeatingUnit.objects.get(pk=repeating_unit)
            assert repeating_unit.check_can_view(chem_db_user), "invalid repeating unit"
            repeating_units.append(repeating_unit)

        image_format = request.GET.get("image_format")
        validate = request.GET.get("validate", True)

        mol = temp_structure.create_mol(repeating_units=repeating_units, terminal_start_group=start_group,
                                        terminal_end_group=end_group, with_chain=True)
        return JsonResponse(
            structure_check(temp_structure, mol=mol, image_format=image_format, validate=validate, raise_error=True)
        )
    except ValidationError as e:
        return JsonResponse({"success": False,
                             "reason": str(e)})
