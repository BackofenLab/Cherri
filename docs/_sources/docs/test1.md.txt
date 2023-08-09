# Welcome

CheRRI detects functional RNA-RNA interaction (RRI) sites, by evaluating if an interaction site most likely occurs in nature.
It helps to filter interaction sites generated either experimentally or by an RRI prediction algorithm by removing false positive interactions.

It is an open source project [hosted on GitHub](https://github.com/BackofenLab/Cherri).

## CheRRI's workflow

CheRRI can be run in two modes, the model generation or **train** mode, or the RRI evaluation or **eval** mode. 
For the evaluation of a given set of RRI sites, a model must be specified in CheRRI's **eval** mode. Here pre-trained models can be applied or the user trains a model using CheRRI's **train** mode.
To train a novel model, an RNA-RNA interactome dataset specifying all RRI sites should be provided. CheRRI makes use of replicate data by checking if an RRI site can be found in all replicates within an overlap threshold. This is how CheRRI builds the set of trusted RRIs. In **eval** mode, the interaction positions are reformatted in the same way as the trusted RRIs. In both modes, CheRRI uses the same core method to generate a feature set which will be used to select a model in **train** mode, and in **eval** mode for the evaluation of the biological relevance of the submitted RRI sites.


[<img src="./assets/images/CheRRI-workflow.png" width="70%" />](./plots/CheRRI-workflow.png)




# Installation

CheRRI was developed in Linux and tested on Ubuntu (18.04 LTS). Conda is required to install CheRRI.


## Install Conda

If you do not have Conda yet, you can e.g. install miniconda, a free + lightweight Conda installer. Get miniconda [here](https://docs.conda.io/en/latest/miniconda.html), choose the newest Python 3 Miniconda3 Linux 64-bit installer and follow the installation instructions. In the end, Conda should be accessed on the command line via (note that your version can be different):

```
$ conda --version
conda 4.10.3
```



## Create environment manually

To manually install CheRRI, first create a Conda environment:

```
conda create -n cherri python=3.8 -c conda-forge -c bioconda
conda activate cherri
```
Inside the environment, you need to install the following dependencies:


```
conda install -c conda-forge scikit-learn
conda install -c conda-forge networkx
conda install -c bioconda bedtools
conda install -c conda-forge biopython
conda install -c conda-forge interlap
conda install -c bioconda intarna
conda install -c conda-forge numpy
conda install -c conda-forge pandas
conda install -c conda-forge eden-kernel
conda install -c conda-forge biofilm
conda install -c conda-forge python-wget
```

Or create the environment with all dependencies at once:

```
conda create -n cherri -c conda-forge -c bioconda scikit-learn networkx numpy bedtools biopython interlap pandas intarna eden-kernel biofilm python-wget
conda activate cherri
```


You additionally need to set a fixed python hash seed within the conda environment:

```
conda env config vars set PYTHONHASHSEED=31337
```
Or set it just for your current session:

```
export PYTHONHASHSEED=31337
```
After setting the environment variable, reactivate your environment:
```
conda deactivate
conda activate cherri
```

## Manual installation

To install the tool itself, simply clone the repository and execute the installation script inside the cloned folder:

```
git clone https://github.com/BackofenLab/Cherri.git 
```

Make sure you are inside the CheRRI's conda environment and you are inside the tools folder. 

```
conda acivate cherri
cd Cherri
```

Than you can install CheRRI.

```
python -m pip install . --ignore-installed --no-deps -vv 
```


Now you can run CheRRI from any given folder:

```
cherri -h
```

## Install using pip
If you don't want to download the CheRRI's git folder you can also use the pipy package. 
```
pip install cherri
```

## Install CheRRI Conda package

```
conda install -c bioconda cherri

```
