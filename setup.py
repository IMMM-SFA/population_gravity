"""
Spatial population downscaling model

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files
"""

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
    name='population_gravity',
    version='1.0.2',
    packages=find_packages(),
    url='https://github.com/IMMM-SFA/population_gravity',
    include_package_data=True,
    license='BSD 2-Clause',
    author='Hamidreza Zoraghein, Chris R. Vernon',
    author_email='hzoraghein@popcouncil.org, chris.vernon@pnnl.gov',
    description='A model to downscale state-level urban and rural populations to a 1 km grid',
    python_requires='>=3.6.*'
)
