"""
Initializaition file for package `disktemperature`.
"""
from setuptools_scm import get_version as __get_version
__version__ = __get_version(root='..', relative_to=__file__)

from . import routines

__all__ = ['routines']
