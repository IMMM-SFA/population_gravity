"""
Spatial population downscaling model

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files
"""


class VersionError(Exception):
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)


try:
    from setuptools import setup, find_packages
except ImportError:
    print("Must have setuptools installed to run setup.py. Please install and try again.")
    raise


def readme():
    with open('README.md') as f:
        return f.read()


def get_requirements():
    with open('requirements.txt') as f:
        return f.read().split()


setup(
    name='spatial_population_downscaling_model',
    version='1.0.0',
    packages=find_packages(),
    url='https://github.com/IMMM-SFA/spatial_population_downscaling_model',
    license='BSD 2-Clause',
    author='Hamidreza Zoraghein',
    author_email='Hamidreza.Zoraghein@du.edu',
    description='A model to allocate urban and rural populations for a defined region to a grid',
    python_requires='>=3.3.*, <4'
)
