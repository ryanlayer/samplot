import re

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("src/__init__.py", "r") as fd:
    version = re.search(
        r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]', fd.read(), re.MULTILINE
    ).group(1)

with open("requirements.txt", "r") as f:
    requires = f.read().splitlines()


setuptools.setup(
    name="samplot",
    version=version,
    description="plotting package for genomic structural variation",
    long_description=long_description,
    long_description_content_type='text/markdown',
    author="Jonathan Belyeu",
    author_email="jrbelyeu@gmail.com",
    url="https://github.com/ryanlayer/samplot.git",
    packages=["samplot"],
    package_data={"": ["LICENSE", "README.md"]},
    package_dir={"samplot": "src"},
    include_package_data=True,
    install_requires=requires,
    license="MIT",
    zip_safe=False,
    entry_points={"console_scripts": ["samplot = samplot.__main__:main"]},
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
