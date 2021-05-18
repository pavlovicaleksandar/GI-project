from distutils.core import setup
import setuptools
setup(
    name='giseed',
    packages=['genomics'],
    py_modules=['genomics.main'],
    version='1.1.0',
    license='MIT',
    description='Tool developed for GI course on ETF Belgrade',
    author='Marko Aleksandar Marko',
    author_email='',
    url='https://github.com/accchi/GI-project',
    keywords=['bwt', 'genomics', 'seed&extend'],
    install_requires=[
      'biopython==1.78',
      'numpy==1.20.3',
      'pytest==6.2.4',
      'click==8.0.0'
    ],
    entry_points={
      'console_scripts': [
        'giseed=genomics.main:cli'
      ]
    }
)
