import os
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

# ssw-py: source files to include in installation for tar.gz
PACKAGE_PATH = os.path.abspath(os.path.dirname(__file__))
MODULE_PATH = os.path.join(PACKAGE_PATH, "ssw")
LIB_PATH = os.path.join(MODULE_PATH, "lib")
ssw_files = [
    os.path.relpath(os.path.join(root, f), MODULE_PATH)
    for root, _, files in os.walk(LIB_PATH)
    for f in files
    if ((".h" in f) or (".c" in f) or (".cpp" in f))
]

extra_compile_args = ["-Wno-unused-function"]
extensions = [
    Extension(
        "indelpost.variant",
        ["indelpost/variant.pyx"],
        include_dirs=pysam_get_include(),
        extra_compile_args=extra_compile_args,
    ),
    Extension(
        "indelpost.utilities",
        ["indelpost/utilities.pyx"],
        include_dirs=pysam_get_include(),
        extra_compile_args=extra_compile_args,
    ),
    Extension(
        "indelpost.pileup",
        ["indelpost/pileup.pyx"],
        include_dirs=pysam_get_include(),
        extra_compile_args=extra_compile_args,
    ),
    Extension(
        "indelpost.varaln",
        ["indelpost/varaln.pyx"],
        include_dirs=pysam_get_include(),
        extra_compile_args=extra_compile_args,
    ),
    Extension(
        "indelpost.contig",
        ["indelpost/contig.pyx"],
        include_dirs=pysam_get_include(),
        extra_compile_args=extra_compile_args,
    ),
    Extension(
        "indelpost.local_reference",
        ["indelpost/local_reference.pyx"],
        include_dirs=pysam_get_include(),
        extra_compile_args=extra_compile_args,
    ),
    Extension(
        "indelpost.softclip",
        ["indelpost/softclip.pyx"],
        include_dirs=pysam_get_include(),
        extra_compile_args=extra_compile_args,
    ),
    Extension(
        "indelpost.localn",
        ["indelpost/localn.pyx"],
        include_dirs=pysam_get_include(),
        extra_compile_args=extra_compile_args,
    ),
    Extension(
        "indelpost.gappedaln",
        ["indelpost/gappedaln.pyx"],
        include_dirs=pysam_get_include(),
        extra_compile_args=extra_compile_args,
    ),
    Extension(
        "indelpost.sswpy",
        sources=["indelpost/sswpy.pyx", "indelpost/ssw.c"],
        extra_compile_args=extra_compile_args,
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
    package_data={"ssw": ssw_files},
    cmdclass={"build_ext": build_ext},
    ext_modules=cythonize(extensions, compiler_directives={"language_level": "3"}),
    install_requires=["numpy>=1.16.0", "pysam>=0.15.0", "cython>=0.29.12"],
)
