try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup
import re
v_file = "cnv_intersect/version.py"
v_line = open(v_file, "rt").read()
v_re = r"^__version__ = ['\"]([^'\"]*)['\"]"
match = re.search(v_re, v_line, re.M)
if match:
    verstr = match.group(1)
else:
    raise RuntimeError("Unable to find version string in {}.".format(v_file))

setup(
    name="cnv_intersect",
    packages=["cnv_intersect"],
    version=verstr,
    description="Find overlapping CNVs",
    author="David A. Parry",
    author_email="david.parry@igmm.ed.ac.uk",
    url='https://github.com/david-a-parry/cnv_intersect',
    download_url='https://github.com/david-a-parry/cnv_intersect/archive/{}.tar.gz'.format(verstr),
    license='MIT',
    install_requires=['pysam'],
    scripts=["bin/cnv_intersect"],
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        ],
)
