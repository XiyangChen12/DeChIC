# setup.py
from setuptools import setup, find_packages
from pathlib import Path

HERE = Path(__file__).resolve().parent

def read_readme():
    p = HERE / "README.md"
    return p.read_text(encoding="utf-8") if p.exists() else ""

setup(
    name="SCENE",
    version="1.0",
    author="Xiyang Chen",
    author_email="chenxiyang12@outlook.com",
    description="A snakemake workflow for single-cell DNA Deaminase based chromatin immuno-conversion sequencing data analysis",
    long_description=read_readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/XiyangChen12/SCENE",
    packages=find_packages(),  # will find the SCENE/ package
    include_package_data=True, # include files from MANIFEST.in
    python_requires=">=3.7",
    install_requires=[
        "snakemake>=7.19.1",
    ],
    entry_points={
        "console_scripts": [
            "SCENE=SCENE.run:main",  # run.py is now inside the SCENE/ package
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ],
)
