"""Legacy editable-install entry point for pip versions before PEP 660 support."""

from setuptools import find_packages, setup

setup(
    name="evotrace",
    version="0.1.0",
    description="Reproducible evolutionary sequence analysis",
    package_dir={"": "src"},
    packages=find_packages("src"),
    python_requires=">=3.9",
    entry_points={"console_scripts": ["evotrace=evotrace.cli:main"]},
)
