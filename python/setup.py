from setuptools import setup, find_packages

setup(
    name='aeolis',
    version='0.0',
    author='Bas Hoonhout',
    author_email='b.m.hoonhout@tudelft.nl',
    packages=find_packages(),
    description='AeoLiS toolbox',
    long_description=open('README.txt').read(),
    install_requires=[
        'matplotlib',
        'netCDF4',
        'pandas',
        'numpy'
    ],
    entry_points={'console_scripts': [
        '{0} = aeolis:cmd'.format(
            'aeolis'),
    ]},
)
