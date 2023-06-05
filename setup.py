from setuptools import setup, find_packages

setup(
    name="HADDOCK_plotter",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        'numpy>=1.19.2',
        'pandas>=1.1.3',
        'matplotlib>=3.3.2',
        'seaborn>=0.11.0'
    ],
)
