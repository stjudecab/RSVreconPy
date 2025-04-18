#!/bin/bash

# Define variables
ENV_NAME="RSVreconEnv"
YML_FILE="RSV_env.yml"

# 1. Check if conda/mamba is installed
if ! command -v mamba &> /dev/null; then
    if ! command -v conda &> /dev/null; then
        echo "Error: Neither conda nor mamba found. Please install Anaconda/Miniconda first."
        exit 1
    else
        CONDA_CMD="conda"
        echo "Warning: Mamba not found, will use slower conda instead."
    fi
else
    CONDA_CMD="mamba"
fi

# 2. Create environment
echo "Creating conda environment '$ENV_NAME'..."
$CONDA_CMD env create -n $ENV_NAME -f $YML_FILE || {
    echo "Failed to create conda environment"
    exit 1
}

# 3. Activate environment
echo "Environment created successfully!"
echo -e "\nUse this command to activate the environment:"
echo "conda activate $ENV_NAME"
echo -e "\nTo remove this environment, use:"
echo "conda env remove -n $ENV_NAME"