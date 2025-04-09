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
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy",
        "pandas",
        "matplotlib",
        "scipy",
        "plotly",
        "biopython",
        "requests", 
        "h5py",
        "importlib-resources; python_version < '3.9'",  # For importlib.resources compatibility
    ],
)