#!/usr/bin/env python
# -*- coding: utf-8 -*-

from os import path
import setuptools, datetime
CMD='tobacco'
NAME='tobacco'
SHORT_CMD='tobacco'

readme_file = path.join(path.dirname(path.abspath(__file__)), 'README.md')
try:
    from m2r import parse_from_file
    readme = parse_from_file(readme_file)
except ImportError:
    with open(readme_file) as f:
        readme = f.read()

today = datetime.date.today().strftime("%b-%d-%Y")
with open(path.join(NAME, '_date.py'), 'w') as fp :
    fp.write('date = \'%s\'' % today)

install_requires=['networkx==2.2','numpy==1.19.1', 'monty==3.0.4']

setuptools.setup(
    name=NAME,
    use_scm_version={'write_to': 'tobacco/_version.py'},
    setup_requires=['setuptools_scm'],
    author="Ryther Anderson,Yamil Colón, Diego Gómez-Gualdrón",
    maintainer="haidi Wang",
    description="MOF generator",
    long_description=readme,
    long_description_content_type="text/markdown",
    url="",
    python_requires="~=3.6",
    packages=['tobacco'],
    #packages_data=['tobacco/RCSRnets-2019-06-01.cgd'],
    classifiers=[
        "Programming Language :: Python :: 3.6"
    ],
    keywords='MOF',
    install_requires=install_requires,    
        entry_points={
          'console_scripts': [
              SHORT_CMD+'= tobacco.main:main']
   }
)
