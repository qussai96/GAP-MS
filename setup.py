from pathlib import Path

from setuptools import setup, find_packages
from setuptools.command.build_py import build_py


class BuildPyWithModel(build_py):
    """Place the repository-level pretrained model inside the built package."""

    def run(self):
        super().run()
        source = Path(__file__).parent / "xgb_protein_classifier.pkl"
        destination = Path(self.build_lib) / "gapms" / source.name
        self.mkpath(str(destination.parent))
        self.copy_file(str(source), str(destination))

    def get_outputs(self, include_bytecode=1):
        outputs = super().get_outputs(include_bytecode)
        outputs.append(str(Path(self.build_lib) / "gapms" / "xgb_protein_classifier.pkl"))
        return outputs

with open("README.md", "r") as f:
    description=f.read()
setup(
    name='gapms',
    version='0.2.1',
    packages=find_packages(),
    package_data={
        'gapms': ['*.pl'],
        'gapms._peplocus': ['LICENSE'],  # Include all Perl scripts in the gapms package
    },
    include_package_data=True,
    cmdclass={'build_py': BuildPyWithModel},
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
        'pyahocorasick>=2.0',
        'pysam>=0.22',
    ],
    entry_points={
        'console_scripts': [
            'gapms=gapms.cli:main',
        ],
    },
    long_description=description,
    long_description_content_type="text/markdown",
)
