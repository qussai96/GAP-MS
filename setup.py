from setuptools import setup, find_packages

with open("README.md", "r") as f:
    description=f.read()
setup(
    name='gapms',
    version='0.1.3',
    packages=find_packages(),
    install_requires=[
        'pandas',
        'numpy',
        'biopython',
        'matplotlib',
        'seaborn',
        'scipy',
        'scikit-learn',
        'xgboost',
        'shap',
        'psauron',
    ],
    entry_points={
        'console_scripts': [
            'gapms=gapms.main:main',
        ],
    },
    long_description=description,
    long_description_content_type="text/markdown",
)
