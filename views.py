from django.shortcuts import render
from django.http import HttpResponse
from django.template import Context, loader, RequestContext, Template

def index(request):
    context = {}
    if request.user.is_authenticated:
        context["user_authenticated"]=True
        context["username"]=request.user.username
    return render(request, "translate_DNA/index.html", context)

# This function activates the cgi script.
def results(request):
    if request.method == 'POST':
        # Process data a bit
        data = request.POST
        
	# Read file in chunks if it exists.
        userinput = data['userinput']
        resolvecharacter = data['resolvecharacter'] 
        if (resolvecharacter == ''): resolvecharacter = 'X'

        flagselection = int(data['flagselection'])
        readingframe = int(data['readingframe'])
       
        # I dunno why its a string of a bool...
        if 'highlight' in data:
                highlight = data['highlight']
        else:
                highlight = 'False'

        if "runtranslate" in data:
                button = "run"
        else:
                button = "none"  # what... ?
        
        # Run actual calulation (by passing data)
        from . import translate_DNA
        output_t = translate_DNA.run(userinput, resolvecharacter, flagselection, highlight, readingframe, button)
        template = Template(output_t)
        context = RequestContext(request)
        return HttpResponse(template.render(context))
    else:
        return HttpResponse("Please use the form to submit data.")
