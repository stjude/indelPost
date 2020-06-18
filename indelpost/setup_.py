#from distutils.core import setup, Extension
#from distutils.core import Extension
from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext


from pysam import get_include as pysam_get_include


extensions = [
    Extension(
        "variant",
        ["variant.pyx"],
        include_dirs=pysam_get_include(),
    ),
    Extension(
        "utilities",
        ["utilities.pyx"],
        include_dirs=pysam_get_include(),
    ),
    
    Extension(
        "pileup",
        ["pileup.pyx"],
        include_dirs=pysam_get_include(),
    ),
]



setup(
    cmdclass = {"build_ext": build_ext},
    ext_modules=cythonize(extensions),
)

# setup(ext_modules=cythonize("./variant_cy.pyx"))
# setup(ext_modules=cythonize("./utilities_cy.pyx"))
