import pathlib
from setuptools import setup, find_packages

HERE = pathlib.Path(__file__).parent

VERSION = '0.1.0'
PACKAGE_NAME = 'adlib'
AUTHOR = 'Sevy Harris'
AUTHOR_EMAIL = 'sevy.harris@gmail.com'
URL = 'https://github.com/sevyharris/adlib'

LICENSE = 'MIT License'
DESCRIPTION = 'Manager for computing dft adsorption energy in ASE and Quantum Espresso'
LONG_DESCRIPTION = (HERE / "README.md").read_text()
LONG_DESC_TYPE = "text/markdown"

INSTALL_REQUIRES = []

setup(
    name=PACKAGE_NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    long_description_content_type=LONG_DESC_TYPE,
    author=AUTHOR,
    license=LICENSE,
    author_email=AUTHOR_EMAIL,
    url=URL,
    install_requires=INSTALL_REQUIRES,
    packages=find_packages()
)
