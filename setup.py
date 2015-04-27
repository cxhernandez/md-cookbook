"""
 md-cookbook: Get cookin' with MD simulations!

 md-cookbook is a python script library that allows users to easily
 setup molecular dynamics (MD) trajectories using OpenMM.
"""

from setuptools import setup, find_packages
from glob import glob


classifiers = """\
    Development Status :: 3 - Alpha
    Intended Audience :: Science/Research
    License :: OSI Approved :: Apache Software License
    Programming Language :: Python
    Programming Language :: Python :: 2.6
    Programming Language :: Python :: 2.7
    Operating System :: Unix
    Operating System :: MacOS
    Operating System :: Microsoft :: Windows
    Topic :: Scientific/Engineering
    Topic :: Scientific/Engineering :: Information Analysis"""

setup(
    name="mdcookbook",
    version="0.1",
    packages=find_packages(),
    scripts=glob('./scripts/*'),
    zip_safe=True,
    platforms=["Windows", "Linux", "Mac OS-X", "Unix"],
    classifiers=[e.strip() for e in classifiers.splitlines()],
    author="Carlos Xavier Hernandez",
    author_email="cxh@stanford.edu",
    description="Get cooking with MD simulations!",
    license="MIT",
    keywords="molecular dynamics simulation",
    url="http://github.com/cxhernandez/md-cookbook"
)
