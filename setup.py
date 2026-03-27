"""Setup for helm-notation-toolkit. Primary config is in pyproject.toml."""
from setuptools import setup, find_packages

setup(
    name="helm-notation-toolkit",
    version="2.0.0",
    description="Python toolkit for HELM (Hierarchical Editing Language for Macromolecules) notation",
    python_requires=">=3.9",
    packages=find_packages(include=["helm_notation*"]),
    package_data={
        "helm_notation.resources": ["*.xml", "*.xsd", "*.properties"],
    },
    install_requires=[
        "rdkit",
        "lxml>=4.9",
    ],
    extras_require={
        "dev": ["pytest>=7.0", "pytest-cov"],
    },
    entry_points={
        "console_scripts": [
            "helm-toolkit=helm_notation.tools.cli:main",
        ],
    },
)
