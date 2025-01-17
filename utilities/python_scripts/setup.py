from setuptools import setup

setup(
    name='python_scripts',
    version='0.1.0',
    packages=['python_scripts', 'python_scripts.a2a', 'python_scripts.nanny', 'python_scripts.nanny.runio',
              'python_scripts.nanny.runio.hadrons', 'python_scripts.processing'],
    package_dir={'': 'src'},
    url='',
    license='',
    author='Michael Lynch',
    author_email='michaellynch628@gmail.com',
    description='Nanny, postprocessing, and A2A contraction scripts'
)
