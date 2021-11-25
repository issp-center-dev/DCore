"""
Setup script for dcore_backend
"""
from setuptools import setup, find_packages
import versioneer

REPO_URL = "https://github.com/issp-center-dev/DCore.git"
LONG_DESCRIPTION = ""

setup(
    name='dcore_backend',
    version='0.1',

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
