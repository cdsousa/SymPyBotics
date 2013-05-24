from os.path import exists
from setuptools import setup

setup(
name='sympybotics',
version='0.3',
author='Cristovao D. Sousa',
author_email='crisjss@gmail.com',
description='Robot dynamic symbolic model generator',
license='BSD',
keywords='robot dynamics',
url='http://github.com/cdsousa/sympybotics',
packages=['sympybotics', 'sympybotics.dynident', 'sympybotics.tools'],
long_description=open('README.md').read() if exists('README.md') else '',
classifiers=[
  'Development Status :: 4 - Beta',
  'Intended Audience :: Developers',
  'License :: OSI Approved :: BSD License',
  'Operating System :: OS Independent',
  'Programming Language :: Python',
  'Intended Audience :: Education',
  'Intended Audience :: Manufacturing',
  'Intended Audience :: Science/Research',
  'Topic :: Scientific/Engineering :: Physics',
  'Topic :: Software Development :: Code Generators',
  ],
)