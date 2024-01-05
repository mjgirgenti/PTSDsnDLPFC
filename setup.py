from setuptools import setup, find_packages

setup(
    name='ptsd',
    description='PTSD single cell',
    packages=find_packages(),
    install_requires=[
            'matplotlib',
            'numpy',
            'scikit-learn',
            'scipy',
            'pandas'
        ]
)
