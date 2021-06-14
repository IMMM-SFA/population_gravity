"""
Spatial population downscaling model

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files
"""

from setuptools import setup, find_packages


setup(
    name='population_gravity',
    version='1.1.0',
    packages=find_packages(),
    url='https://github.com/IMMM-SFA/population_gravity',
    license='BSD 2-Clause',
    author='Hamidreza Zoraghein, Chris R. Vernon',
    author_email='hzoraghein@popcouncil.org, chris.vernon@pnnl.gov',
    description='A model to downscale state-level urban and rural populations to a 1 km grid',
    python_requires='>=3.6.*',
    install_requires=[
        'rasterio~=1.1.5',
        'simplejson~=3.17.0',
        'numpy~=1.19.5',
        'pandas~=1.0.5',
        'xarray~=0.16.2',
        'scipy~=1.5.1',
        'pathos~=0.2.6',
        'PyYAML~=5.3.1'
    ],
    include_package_data=True
)
