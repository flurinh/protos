from setuptools import setup, find_packages

setup(
    name="protos",
    version="0.1.0",
    author="flurinh",
    author_email="hidberf@gmail.com",
    description="A protein data processing library",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="[Repository URL]",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    include_package_data=True,  # Include non-Python files
    package_data={
        "protos": [
            "reference_data/README.md",
            "reference_data/*/registry.json",
            "reference_data/structure/mmcif/*.cif",
            "reference_data/sequence/fasta/*.fasta",
            "reference_data/grn/tables/*.csv",
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
    ],
)