"""
Setup script for dcorelib
"""
from setuptools import setup, find_packages
import os, re, io

REPO_URL = "https://github.com/issp-center-dev/DCore.git"
LONG_DESCRIPTION = ""

def readfile(*parts):
    """Return contents of file with path relative to script directory"""
    herepath = os.path.abspath(os.path.dirname(__file__))
    fullpath = os.path.join(herepath, *parts)
    with io.open(fullpath, 'r') as f:
        return f.read()

def extract_version(*parts):
    """Extract value of __version__ variable by parsing python script"""
    initfile = readfile(*parts)
    version_re = re.compile(r"(?m)^__version__\s*=\s*['\"]([^'\"]*)['\"]")
    match = version_re.search(initfile)
    return match.group(1)

VERSION = extract_version('src', 'dcorelib', '__init__.py')

setup(
    name='dcorelib',
    version=VERSION,

    description='Backend of DMFT software for CORrelated Electrons',
    keywords=' '.join([
        'condensed-matter',
        'dmft',
        ]),
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Physics',
        'Programming Language :: Python :: 3',
        ],

    url=REPO_URL,
    author=', '.join([
        'H. Shinaoka',
        ]),
    author_email='h.shinaoka@gmail.com',

    python_requires='>=3.6, <4',
    install_requires=[
        'numpy',
        'scipy',
        'h5py',
        'irbasis3'
        ],
    extras_require={
        'dev': ['pytest', 'matplotlib', 'versioneer'],
        },

    setup_requires=[
        'numpy',
        'scipy',
        'h5py',
        'irbasis3'
        ],

    package_dir={'': 'src'},
    packages=find_packages(where='src'),

    zip_safe=False,
    )
