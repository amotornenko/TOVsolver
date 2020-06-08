from __future__ import absolute_import

from setuptools import setup, find_packages
from codecs import open


with open('README.md', encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='tovsolver',
    version='0.1.0',
    author='Anton Motornenko (FIAS)',
    author_email='motornenko@fias.uni-frankfurt.de',
    description='TOV solver',
    long_description=long_description,
    long_description_content_type='text/x-rst',
    url='https://github.com/amotornenko/TOVsolver',
    license='GPLv3',
    packages=find_packages(exclude=['doc', 'test']),
    include_package_data=True,
    package_data = {
    '' : ['*.dat'],
    },
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
    ],
    extras_require = {
    },
    #data_fies = [('tovsolver', [[os.path.join(data_dir, "data1")]
    classifiers=[
        'Development Status :: 1 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Natural Language :: English',
        'Operating System :: OS Independent',
    ])
