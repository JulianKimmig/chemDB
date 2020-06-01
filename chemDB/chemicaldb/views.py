from django.shortcuts import render


# Create your views here.

def main_view(request):
    return render(request, "chemicaldb/main.html")


def substance_view(request,pk):
    return None