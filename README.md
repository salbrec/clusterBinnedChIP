# Bin ChIP-seq peaks before applying hierarchical clustering!

According to the results presented in our research paper (link is coming) the clustering quality can significantly be improved when applying further preprocessing to peak profiles from ChIP-sequencing. This preprocessing, in particular, is called binning and describes the process of cutting the genome into adjacent non-overlapping windows *bins* that are annotated if there is an overlapping peak. 

On a set of 40 CTCF profiles from the ENCODE data portal, we could show that a bin-based clustering results in a three-fold higher adjusted mutual information, a clustering evaluation metric, compared to using the peaks.

<img src="utils/readme/AMI.gif" width="1000">

The code for running the clustering based on the default peaks or applying the binning prior to clustering is provided by this repository. 

Below you’ll find example runs for reproducing the clustering for the CTCF example. Furthermore, at the bottom of this README, we discuss in detail the installation options you have to make this code run. In short, having a Python 3 installation together with bedtools should do it.

## Running Examples

Run ```python binNclusterChIP.py --help``` to print detailed information about the parameters.

The script takes a table in tsv format in which the profiles are described together with a file path for each profile in bed file-format. This table must have the columns **ID**, **filePath** and **name**. The ID is used e.g. to access the binned profile files that are created by the script. The name is shown in the clustermap and the filePath is used to read in a certain profile. An example for CTCF is provided by "CTCF_files.tsv". The column **label** is optional. However, if you can label your profiles, then the agglomerative clustering is also validated at the end of the scritp.

With the following line one can create a clustermap based on the peaks. Jaccard-Distances are automatically computed. 

```
python binNclusterChIP -p CTCF_files.tsv -t peaks -o ./example_CTCF/
```

When using bins, one must also define the reference genome and the bin size:

```
binNclusterChIP -p CTCF_files.tsv -t bins -g GRCh38 -s 5kb -o ./example_CTCF/
```

**BED files must ALWAYS be SORTED before running the script!** 

You can use *bedtools* to sort your peak files (see [here](https://bedtools.readthedocs.io/en/latest/content/tools/sort.html))

## Runtime

It can be time consuming to get the binned profiles. Based on the CTCF example provided here, we got observed average approx. 30 seconds per profile for a 5kb binning. For 2kb resolution the time increases to approx. 1 minute. However, once the bins are derived from the peaks, the Jaccard-Distance is quickly computed based on the binned profiles. When using the peaks, there is no need to spend computational time on the binning, however, the Jaccard-Distance computation takes more time on the peak profiles.

## Installation Guidelines

Most of the Python packages required for this tool, are commonly used in the NGS community. Additionally, one might has to install *bedtools*. However, we provide some tips for setting up an Anaconda environment that has evrything needed to run the code.

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



