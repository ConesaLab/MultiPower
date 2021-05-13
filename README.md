# MultiPower and MultiML: statistical power studies for multi-omics experiments

## MultiPower

The MultiPower R method performs statistical power studies for multi-omics experiments, 
and is designed to assist users in experimental design as well as in the evaluation of already-generated multi-omics datasets. 
More details on the method can be found in our manuscript [[1]](#1) and in the 
[MultiPower User’s Guide](https://github.com/ConesaLab/MultiPower/blob/master/MultiPower_UsersGuide.pdf).

MultiPower is available as an R package, and can be installed as follows:

```
install.packages(“devtools”)
devtools::install_github(“ConesaLab/MultiPower”)
```

### Installing MultiPower dependencies

Some dependencies are required before running MultiPower. The following can be installed from R via 
`install.packages()` and loaded with `library()`:

- slam
- Lpmodeler
- Rsymphony

MultiPower also requires the *RnaSeqSampleSize* package, which can be installed following 
the instructions in the corresponding [Bioconductor repository](https://bioconductor.org/packages/release/bioc/html/RnaSeqSampleSize.html).

**NOTE for Linux users:** when working on a Linux system, installing *Rsymphony* will require some additional steps. 
First, the last version of SYMPHONY needs to be installed:

```
$ svn checkout https://projects.coin-or.org/svn/SYMPHONY/releases/5.6.16
SYMPHONY-5.6.16
$ cd SYMPHONY-5.6.16
$ ./configure
$ make
$ make install

```

Additional required libraries can then be installed as follows:

```
$ sudo apt-get install coinor-libcgl-dev coinor-libclp-dev
coinor-libcoinutils-dev coinor-libosi-dev
$ sudo apt-get install coinor-libsymphony-dev
$ sudo apt-get install autotools-dev
```

Next, download the R package from the [CRAN repository](https://cran.r-project.org/src/contrib/Rsymphony_0.1-26.tar.gz) 
and install from R:
 
```
install.packages(“Rsymphony_0.1-26.tar.gz”, repos = NULL)
```

After installation is complete, please refer to MultiPower’s 
[User Guide](https://github.com/ConesaLab/MultiPower/blob/master/MultiPower_UsersGuide.pdf) for instructions on 
how to use the different functions included in the package.

## MultiML

The MultiML method is included as a complementary tool to MultiPower, 
and is designed to help users determine the optimal sample size required to control for 
classification error rates when using one or more omics datasets. 
Details on the MultiML algorithm and its applications can be found in our manuscript [[1]](#1). 

If you are interested in using MultiML for your research, please see this folder(link) 
for scripts and instructions. For detailed information on how to run the tool, please read 
[MultiML's  User Guide](https://github.com/ConesaLab/MultiPower/blob/master/MultiPower_UsersGuide.pdf)



## Citation

If you are using MultiPower or MultiML in your research, please cite the following publication:

<a id="1">[1]</a>
Tarazona, S., Balzano-Nogueira, L., Gómez-Cabrero, D. et al. 
Harmonization of quality metrics and power calculation in multi-omic studies. 
Nat Commun 11, 3092 (2020). https://doi.org/10.1038/s41467-020-16937-8 

