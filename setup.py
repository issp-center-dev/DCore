"""
Setup script for DCore
"""
from setuptools import setup, find_packages
import versioneer
import re, os, io

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

setup(
    name='dcore',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),

    description='DMFT software for CORrelated Electrons',
    #long_description=LONG_DESCRIPTION,
    #long_description_content_type='text/markdown',
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
        'J. Otsuki',
        'K. Yoshimi',
        'M. Kawamura',
        'N. Takemori',
        'Y. Motoyama'
        ]),
    author_email='h.shinaoka@gmail.com',

    python_requires='>=3.6, <4',
    install_requires=[
        'numpy',
        'scipy',
        # h5py 2.10.0 has a bug.
        # Import h5py imports mpi4py automatically.
        'h5py!=2.10.0',
        'dcorelib',
        ],
    extras_require={
        'dev': ['pytest', 'sphinx', 'matplotlib', 'wild_sphinx_theme', 'versioneer'],
        },

    setup_requires=[
        'numpy',
        'scipy',
        ],

    package_dir={'': 'src'},
    packages=find_packages(where='src'),
    entry_points={
        "console_scripts": [
            "dcore_pre = dcore.dcore_pre:run",
            "dcore = dcore.dcore:run",
            "dcore_post = dcore.dcore_post:run",
            "dcore_check = dcore.dcore_check:run",
            "dcore_bse = dcore.dcore_bse:run",
            "dcore_gk = dcore.dcore_gk:run",
            "dcore_mpicheck = dcore.dcore_mpicheck:run",
        ]
    },

    zip_safe=False,
    )
