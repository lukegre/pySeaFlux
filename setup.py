from os.path import exists

from setuptools import find_packages, setup

if exists("README.rst"):
    with open("README.rst") as f:
        long_description = f.read()
else:
    long_description = ""

CLASSIFIERS = [
    "Development Status :: 3 - Alpha",
    "License :: OSI Approved :: MIT License",
    "Operating System :: OS Independent",
    "Intended Audience :: Science/Research",
    "Programming Language :: Python :: 3 :: Only",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
]

setup(
    name="SeaFlux",
    author="Luke Gregor",
    author_email="lukegre@gmail.com	",
    description="Calculate sea fluxes",
    keywords="SeaFlux",
    license="GNUv3",
    classifiers=CLASSIFIERS,
    url="https://github.com/luke-gregor/SeaFlux",
    use_scm_version={"version_scheme": "guess-next-dev", "local_scheme": "dirty-tag"},
    long_description=long_description,
    long_description_content_type="text/x-rst",
    packages=find_packages(),
    install_requires=["numpy", "pandas", "xarray", "pooch"],
    test_suite="tests",
    tests_require=["pytest-cov"],
    setup_requires=[
        "wheel",
        "setuptools_scm",
        "setuptools>=30.3.0",
        "setuptools_scm_git_archive",
    ],
)
