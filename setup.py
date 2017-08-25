"""
Setup file for the drift speed package.
"""
from numpy.distutils.core import Extension, setup, find_packages
import os

setup(
    name='calc_drift',
    use_scm_version=True,
    description='calculates drift velocity',
    long_description=open(os.path.join(
        os.path.dirname(__file__), 'Readme.md')).read(),
    url='http://www.til-birnstiel.de',
    author='Til Birnstiel',
    author_email='birnstiel@me.com',
    packages=find_packages(),
    license='GPLv3',
    include_package_data=True,
    setup_requires=['setuptools_scm'],
    install_requires=['numpy'],
    zip_safe=False,
    ext_modules=[
        Extension(
            name='calc_drift.routines',
            sources=[
                'calc_drift/routines.f90',
                'calc_drift/constants.f90',
                'calc_drift/switches.f90',
                'calc_drift/variables.f90',
                ],
            include_dirs='calc_drift/'
            )
        ],
    )
