from __future__ import print_function
import pytriqs
import os
path = pytriqs.__path__[0]
print(str('/'.join(os.path.abspath(path).split('/')[0:-1])), end='')
