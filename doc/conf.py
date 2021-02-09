# -*- coding: utf-8 -*-

import sys

extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.mathjax',
              'sphinx.ext.intersphinx',
              'matplotlib.sphinxext.plot_directive',
              'sphinx.ext.doctest',
              'sphinx.ext.todo',
              'sphinx.ext.viewcode',
              'sphinx.ext.autosummary',
]

source_suffix = '.rst'
todo_include_todos = True

project = u'DCore'
copyright = u'2017, The University of Tokyo'

html_theme = 'wild'
import wild_sphinx_theme
html_theme_path = [wild_sphinx_theme.get_theme_dir()]

html_show_sphinx = False
html_context = {'header_title': 'DCore',
                'header_subtitle': 'integrated DMFT software for CORrelated Electrons',
                'header_links': [['Install', 'install'],
                                 ['Documentation', 'documentation'],
                                 #['Presentatation', 'presentation'],
                                 ['Issues', 'issues'],
                                 ['About DCore', 'about']]}
#html_static_path = ['@CMAKE_SOURCE_DIR@/doc/_static']
html_static_path = ['_static']
#html_sidebars = {'index': ['sideb.html', 'searchbox.html']}
html_sidebars = {'**': ['globaltoc.html', 'relations.html', 'searchbox.html']}

#html_title = "DMFT software"
#html_logo = 'logo_dcore1.png'

# no 'module' in the header
html_domain_indices = False

# no 'index' in the header
html_use_index = False

htmlhelp_basename = 'DCoredoc'

rst_epilog = '.. |BRANCH| replace:: {}'.format('@GIT_BRANCH_NAME@')

#intersphinx_mapping = {'python': ('http://docs.python.org/2.7', None), 'triqslibs': ('http://triqs.ipht.cnrs.fr/1.x', None), 'triqscthyb': ('https://triqs.ipht.cnrs.fr/1.x/applications/cthyb/', None), 'triqsdfttools': ('https://triqs.ipht.cnrs.fr/1.x/applications/dft_tools/', None)}

# overwrite css of html_theme
def setup(app):
    app.add_css_file('dcore.css')
