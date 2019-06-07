import setuptools
import os

if os.path.isfile("README.md"):
    with open("README.md", "r") as fh:
        long_description = fh.read()

if os.path.isfile("requirements.txt"):
    with open("requirements.txt", 'r') as f:
        install_requires = f.read().splitlines()
print(install_requires)

setuptools.setup(
     name='samplot',  
     version='1.0.3',
     scripts=['src/samplot.py', "src/samplot_vcf.py"] ,
     author="Jonathan Belyeu",
     author_email="jrbelyeu@gmail.com",
     description="plotting package for genomic structural variation",
     long_description=long_description,
     long_description_content_type="text/markdown",
     url="https://github.com/ryanlayer/samplot.git",
     install_requires=install_requires,
     packages=setuptools.find_packages(),
     classifiers=[
         "Programming Language :: Python :: 3",
         "Programming Language :: Python :: 2",
         "License :: OSI Approved :: MIT License",
     ],
 )
