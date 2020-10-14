from distutils.core import setup, Extension

# define the extension module
KLevenshtein = Extension('KLevenshtein', sources=['KLevenshtein.c'])

# run the setup
setup(ext_modules=[KLevenshtein])

