from setuptools import setup, find_packages

setup(
    name='BIP',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'pandas',
        'biopython',
    ],
    entry_points={
        'console_scripts': [
            'run_default_prediction=BIP.get_supported_and_highly_supported_proteins_genes:main',
            'run_filter_prediction=BIP.filter_prediction:main',
            'run_find_and_apply_score_cutoff=BIP.find_and_apply_score_cutoff:main',
            'run_unite_gtf=BIP.unite_gtf:main',
        ],
    },
)
