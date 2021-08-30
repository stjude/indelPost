from setuptools import setup, find_packages
from setuptools.extension import Extension

def pip_install(pkg_name):
    import subprocess

    subprocess.check_call(
        ["python", "-m", "pip", "install", pkg_name], stdout=subprocess.DEVNULL
    )


try:
    from Cython.Build import cythonize
    from Cython.Distutils import build_ext
except ImportError:
    pip_install("cython")

    from Cython.Build import cythonize
    from Cython.Distutils import build_ext

try:
    from pysam import get_include as pysam_get_include
except ImportError:
    pip_install("pysam")
    from pysam import get_include as pysam_get_include

try:
    import numpy
except ImportError:
    pip_install("numpy")


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
        "indelpost.pileup", ["indelpost/pileup.pyx"], include_dirs=pysam_get_include(),
    ),
    Extension(
        "indelpost.varaln", ["indelpost/varaln.pyx"], include_dirs=pysam_get_include(),
    ),
    Extension(
        "indelpost.contig", ["indelpost/contig.pyx"], include_dirs=pysam_get_include(),
    ),
    Extension(
        "indelpost.local_reference", ["indelpost/local_reference.pyx"], include_dirs=pysam_get_include(),
    ),
    Extension(
        "indelpost.softclip",
        ["indelpost/softclip.pyx"],
        include_dirs=pysam_get_include(),
    ),
    Extension(
        "indelpost.localn", ["indelpost/localn.pyx"], include_dirs=pysam_get_include(),
    ),
    Extension(
        "indelpost.gappedaln",
        ["indelpost/gappedaln.pyx"],
        include_dirs=pysam_get_include(),
    ),
]


version = {}
with open("indelpost/version.py") as ver:
    exec(ver.read(), version)

setup(
    name="indelpost",
    version=version["__version__"],
    description="Python library for simple and complex indels",
    url="https://github.com/stjude/indelpost",
    author="Kohei Hagiwara",
    author_email="kohei.hagiwara@stjude.org",
    license="Apache License 2.0",
    packages=find_packages(exclude=["tests"]),
    cmdclass={"build_ext": build_ext},
    ext_modules=cythonize(extensions, compiler_directives={"language_level": "3"}),
    install_requires=["ssw-py", "numpy>=1.16.0", "pysam>=0.15.0", "cython>=0.29.12"],
)
