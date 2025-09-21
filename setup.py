from setuptools import setup, find_packages, Command
from setuptools.command.install import install
from setuptools.command.develop import develop
import sys
import subprocess
import os

this_dir = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_dir, "README.md"), encoding="utf-8") as f:
    long_description = f.read()


def initialize_protos_data():
    """Initialize Protos data directory after installation."""
    try:
        # Import here to ensure package is installed
        from protos.cli.init_data import init_data_directory
        from protos.io.paths.path_config import ProtosPaths
        
        # Get the default data directory
        paths = ProtosPaths()
        data_root = paths.data_root
        
        print(f"\nInitializing Protos data directory at: {data_root}")
        
        # Check if already initialized
        from pathlib import Path
        data_dir = Path(data_root)
        entity_registry = data_dir / "entity_registry.json"
        
        if entity_registry.exists():
            print("Data directory already initialized. Skipping initialization.")
        else:
            # Run initialization
            stats = init_data_directory(force=True)
            
            print(f"\n✅ Protos data directory initialized successfully!")
            print(f"   Location: {data_root}")
            print(f"   Directories created: {stats.get('directories_created', 0)}")
            print(f"   Registries created: {stats.get('registries_created', 0)}")
            print(f"   Reference files copied: {stats.get('reference_files_copied', 0)}")
            
    except Exception as e:
        print(f"\n⚠️  Warning: Could not initialize data directory: {e}")
        print("You can manually initialize it later with: python -m protos.cli.init_data")


class PostInstallCommand(install):
    """Post-installation command to initialize data directory."""
    def run(self):
        install.run(self)
        initialize_protos_data()


class PostDevelopCommand(develop):
    """Post-develop command to initialize data directory."""
    def run(self):
        develop.run(self)
        initialize_protos_data()


setup(
    name="protos",
    version="0.1.0",
    author="flurinh",
    author_email="hidberf@gmail.com",
    description="A protein data analysis library",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="[Repository URL]",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    include_package_data=True,
    package_data={
        "protos": [
            "reference_data/**/*",  # Include all reference data files
            "reference_data/README.md",
            "reference_data/*/registry.json",
            "reference_data/structure/mmcif/*.cif",
            "reference_data/structure/structure_dataset/*.json",
            "reference_data/structure/structure_dataset/standard/*.json",
            "reference_data/sequence/fasta/*.fasta",
            "reference_data/grn/ref/*.csv",
            "reference_data/grn/tables/*.csv",
            "reference_data/grn/configs/*.json",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.10",
    install_requires=[
        "typer",
        "numpy>=1.22.0",
        "pandas>=1.4.0",
        "matplotlib>=3.5.0",
        "scipy>=1.8.0",
        "plotly>=5.6.0",
        "biopython>=1.79",
        "requests>=2.27.0",
        "h5py>=3.6.0",
        "pytest>=7.0.0",
        "tqdm>=4.63.0",
        "scikit-learn>=1.0.2",
        "networkx>=2.7.0",
        "python-dotenv>=0.20.0",
        "gemmi>=0.5.5",
        "seaborn>=0.11.2",
        "openpyxl>=3.1.5",
        "click>=8.0.0",  # Added for CLI commands
    ],
    extras_require={
        "embedding": [
            "transformers>=4.20.0",
        ],
        "gpu": [
            "transformers>=4.20.0",
            "accelerate>=0.20.0",
            "sentencepiece>=0.1.96",
        ],
        "dev": [
            "black>=22.0.0",
            "isort>=5.10.0",
            "mypy>=0.950",
            "flake8>=4.0.0",
        ],
        "all": [
            "transformers>=4.20.0",
            "accelerate>=0.20.0",
            "sentencepiece>=0.1.96",
            "black>=22.0.0",
            "isort>=5.10.0",
            "mypy>=0.950",
            "flake8>=4.0.0",
        ]
    },
    cmdclass={
        'install': PostInstallCommand,
        'develop': PostDevelopCommand,
    },
    entry_points={
        'console_scripts': [
            'protos-init=protos.cli.init_data:init_data',
            'protos-cleanup=protos.cli.cleanup_data:cleanup_data',
            'protos = protos.cli.cli:app'
        ],
    },
)