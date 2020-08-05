#from distutils.core import setup, Extension
#from distutils.core import Extension
from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext


from pysam import get_include as pysam_get_include

extensions = [
    Extension(
        "indelpost.variant",
        ["indelpost/variant.pyx"],
        include_dirs=pysam_get_include(),
    ),
    Extension(
        "indelpost.utilities",
        ["indelpost/utilities.pyx"],
        include_dirs=pysam_get_include(),
    ),
    
    Extension(
        "indelpost.pileup",
        ["indelpost/pileup.pyx"],
        include_dirs=pysam_get_include(),
    ),

    Extension(
        "indelpost.varaln",
        ["indelpost/varaln.pyx"],
        include_dirs=pysam_get_include(),
    ),

    Extension(
        "indelpost.contig",
        ["indelpost/contig.pyx"],
        include_dirs=pysam_get_include(),
    ),

    Extension(
        "indelpost.softclip",
        ["indelpost/softclip.pyx"],
        include_dirs=pysam_get_include(),
    ),

    Extension(
        "indelpost.localn",
        ["indelpost/localn.pyx"],
        include_dirs=pysam_get_include(),
    ),

]



version = {}
with open("indelpost/version.py") as ver:
    exec(ver.read(), version)

setup(
    name="indelpost",
    version=version["__version__"],
    description="Postprocess indel calls",
    url="https://github.com/stjude/indelpost",
    author="Kohei Hagiwara",
    author_email="kohei.hagiwara@stjude.org",
    license="Apache License 2.0",
    packages=find_packages(exclude=["tests"]),
    cmdclass = {"build_ext": build_ext},
    ext_modules=cythonize(extensions, compiler_directives={'language_level' : "3"}),
)

# setup(ext_modules=cythonize("./variant_cy.pyx"))
# setup(ext_modules=cythonize("./utilities_cy.pyx"))
