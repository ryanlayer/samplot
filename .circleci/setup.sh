#!/bin/bash

set -exo pipefail

WORKSPACE=$(pwd)

# Set path
echo "export PATH=$WORKSPACE/anaconda/bin:$PATH" >> $BASH_ENV
source $BASH_ENV

## Passed from .circleci/config.yml (Only 2 or 3 permited)
pythonversion=$1
if (( $pythonversion != 2 && $pythonversion != 3 ))
then
    echo -e "\nERROR: Python 2 or 3 designation required. Python version $pythonversion was supplied. Please correct and run again\n"
    exit 1   
fi 

# setup conda and dependencies 
if [[ ! -d $WORKSPACE/anaconda ]]; then
    mkdir -p $WORKSPACE

    # step 1: download and install anaconda
    if [[ $OSTYPE == darwin* ]]; then
        tag="MacOSX"
        tag2="darwin"
    elif [[ $OSTYPE == linux* ]]; then
        tag="Linux"
        tag2="linux"
    else
        echo "Unsupported OS: $OSTYPE"
        exit 1
    fi  

    curl -O https://repo.continuum.io/miniconda/Miniconda$pythonversion-latest-$tag-x86_64.sh
    sudo bash Miniconda$pythonversion-latest-$tag-x86_64.sh -b -p $WORKSPACE/anaconda/
    sudo chown -R $USER $WORKSPACE/anaconda/

    mkdir -p $WORKSPACE/anaconda/conda-bld/$tag-64

    # step 2: setup channels
    conda config --system --add channels defaults
    conda config --system --add channels r
    conda config --system --add channels bioconda
    conda config --system --add channels conda-forge

    # step 3: install Samplot requirements
    conda install -y --file requirements.txt

fi






