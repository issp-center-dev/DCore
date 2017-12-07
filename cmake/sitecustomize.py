def application_pytriqs_import(name,*args,**kwargs):
    if name.startswith('pytriqs.applications.dcore'):
        name = name[len('pytriqs.applications')+1:]
    return builtin_import(name,*args,**kwargs)

import __builtin__
__builtin__.__import__, builtin_import = application_pytriqs_import, __builtin__.__import__

