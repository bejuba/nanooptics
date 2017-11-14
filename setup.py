# setup.py
from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize
import numpy as np

setup(name='nanooptics',
      version='0.2',
      description='Python module for data evaluation and presentation of numerous quantum optical experiments.',
      author='Bernd Sontheimer',
      author_email='bernd.sontheimer@physik.hu-berlin.de',
      license='GPLv3',
      classifiers=[
            'Development Status :: 3 - Alpha',
            'Intended Audience :: Developers :: Physicists',
            'Topic :: Physics :: Optics :: Quantum',
            'License :: OSI Approved :: GPLv3',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.2',
            'Programming Language :: Python :: 3.3',
            'Programming Language :: Python :: 3.4',
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6',
      ],
      python_requires='>=3',
      install_requires=['cython', 'numpy', 'lmfit'],
      packages=find_packages(exclude=['examples*']),
      include_dirs=[np.get_include()],
      ext_modules=cythonize([Extension('*', ['nanooptics/*.pyx'])]),
      )