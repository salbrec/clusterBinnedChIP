# Bin ChIP-seq peaks before applying hierarchical clustering!

According to the results presented in our research paper (link is coming) the clustering quality can significantly be improved when applying further preprocessing to peak profiles from ChIP-sequencing. This preprocessing, in particular, is called binning and describes the process of cutting the genome into adjacent non-overlapping windows *bins* that are annotated if there is an overlapping peak. 

On a set of 40 CTCF profiles from the ENCODE data portal, we could show that a bin-based clustering results in a three-fold higher adjusted mutual information, a clustering evaluation metric, compared to using the peaks.

The code for running the clustering based on the default peaks or applying the binning prior to clustering is provided by this repository. 

Below you’ll find example runs for reproducing the clustering for the CTCF example. Furthermore, at the bottom of this README, we discuss in detail the installation options you have to make this code run. In short, having a Python 3 installation together with bedtools should do it.

## Running examples ...

BED FILES must be SORTED!!!

## Installation Guidelines

### Anaconda environment

The code is written in Python for Linux systems. Indeed, everything was developed under WSL (Windows Subsystem for Linux) which we recommend for Windows users.

If you like to use Anaconda, you might simply install bedtools (we tested the code on version 2.30.0, however other versions probably work, too). 

However, apply the following steps to create your own environment and add the required packages. This was working perfectly fine on Ubuntu 20.04.5 LTS with Anaconda 22.11.1 

For installing the most recent Anaconda, which we recommend anyways, we refer to the straighforward 3-step installation described on the official [Anaconda website](https://docs.anaconda.com/anaconda/install/linux/)

```
conda create -n chip python=3.9 anaconda pandas scikit-learn

conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda

conda install -c bioconda bedtools
conda install -c anaconda seaborn
```

Or why don't you simply use the *yml* provided in this repository:

```
conda env create -f chip.yml
```

Note that users who already use scikit-learn and seaborn, might only need to install *bedtools* to an excisting environment. 



