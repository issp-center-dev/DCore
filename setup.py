"""
Setup script for irbasis_x
"""
import io
import os
import re

from distutils.dep_util import newer
from distutils import log
from setuptools import setup, find_packages
from setuptools.extension import Extension
from setuptools.command.build_ext import build_ext
from setuptools.command.sdist import sdist


_HEREPATH = os.path.abspath(os.path.dirname(__file__))
_VERSION_RE = re.compile(r"(?m)^__version__\s*=\s*['\"]([^'\"]*)['\"]")
_DOCLINK_RE = re.compile(r"(?m)^\s*\[\s*([^\]\n\r]+)\s*\]:\s*(doc/[./\w]+)\s*$")
_PYXFILE_RE = re.compile(r"(?i)\.pyx$")

def fullpath(path):
    """Return the full path to a file"""
    if path[0] == '/':
        raise ValueError("Do not supply absolute paths")
    return os.path.join(_HEREPATH, *path.split("/"))


def readfile(path):
    """Return contents of file with path relative to script directory"""
    return io.open(fullpath(path), 'r').read()


def extract_version(path):
    """Extract value of __version__ variable by parsing python script"""
    return _VERSION_RE.search(readfile(path)).group(1)


def rebase_links(text, base_url):
    """Rebase links to doc/ directory to ensure they work online."""
    result, nsub = _DOCLINK_RE.subn(r"[\1]: %s/\2" % base_url, text)
    return result


#VERSION = extract_version('src/dcore/__init__.py')
VERSION = '3.0.0b1'
REPO_URL = "https://github.com/issp-center-dev/DCore.git"
#DOCTREE_URL = "%s/tree/v%s" % (REPO_URL, VERSION)
#LONG_DESCRIPTION = rebase_links(readfile('README.md'), DOCTREE_URL)
#LONG_DESCRIPTION = rebase_links(readfile('README.md'), DOCTREE_URL)
LONG_DESCRIPTION = ""

setup(
    name='dcore',
    version=VERSION,

    description='DMFT software for CORrelated Electrons',
    #long_description=LONG_DESCRIPTION,
    #long_description_content_type='text/markdown',
    keywords=' '.join([
        'condensed-matter',
        'dmft',
        ]),
    classifiers=[
        'Development Status :: 5 - Stable',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Physics',
        'Programming Language :: Python :: 3',
        ],

    url=REPO_URL,
    author=', '.join([
        'H. Shinaoka',
        'J. Otsuki',
        'K. Yoshimi',
        'M. Kawamura',
        'N. Takemori',
        'Y. Motoyama'
        ]),
    author_email='h.shinaoka@gmail.com',

    python_requires='>=3, <4',
    install_requires=[
        'numpy',
        'scipy',
        ],
    extras_require={
        'dev': ['pytest', 'sphinx', 'matplotlib', 'wild_sphinx_theme'],
        },

    setup_requires=[
        'numpy',
        'scipy',
        ],

    package_dir={'': 'src'},
    packages=find_packages(where='src'),

    zip_safe=True,      # reconsider when adding data files
    )
