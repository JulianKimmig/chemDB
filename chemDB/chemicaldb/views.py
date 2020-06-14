from django.shortcuts import render


# Create your views here.

def main_view(request):
    return render(request, "chemicaldb/main.html")

STRUCTURE_PLUGINS={}
def register_structure_plugin(name,title,entry_link,image=None):
    STRUCTURE_PLUGINS[name]={'name':name,'title':title,'link':entry_link,'image':image}

def structures_main(request):
    return render(request, "chemicaldb/structures_main.html",context={'structure_plugins':STRUCTURE_PLUGINS})


def substance_view(request,pk):
    return None



def structures_browse(request):
    return render(request, "chemicaldb/structures_browse.html")


def substance_main(request):
    return render(request, "chemicaldb/substance_main.html")