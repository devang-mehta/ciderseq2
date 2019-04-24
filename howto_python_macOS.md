# How to install Python 3 and required modules on macOS X

A detailed guide on our favorite way to locally install Python and the required packages on your local macOS X machine. 

## Install Homebrew

Homebrew is a package manager for macOS, much like `apt-get` on Linux. Please install this. 

To install, simply open your Terminal and type:
```
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
```

## Install Python 3

Once you have Homebrew installed, open Terminal and type:

```
brew install python3
```

## Install MUSCLE
* Download Muscle from:www.drive5.com/muscle/downloads.htm
* Extract the file and copy to your cider-seq folder.
* If placed in another folder, make sure you edit the path address in the `ciderseq_config.json` file.
* For MAC OS, please edit the name of the MUSCLE file in the `ciderseq_config.json` file to muscle3.8.31_i86darwin64 or whichever version of MUSCLE you downloaded. 

## Install BLAST
* Download the latest blast .dmg file from: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
* Double click on the downloaded file to install (required administrator account)

## Install modules

Open Terminal and type:

```
pip install [modulename]
```

for example:

```
pip install biopython
```

will install Biopython.

Please repeat and use pip to install all the required modules for ciderseq (see [Readme](README.md))


**Remember, if you already have Python 2 installed, use:**

```
pip3 install [module]
```

to install the module on your Python 3 installation. 

Similarly, to run python commands using Python 3, type:

```
python3 [command]
```

for example:

```
python3 ciderseq.py --format fastq examples/ciderseq_config.json examples/example1.fastq
```

will run `ciderseq` using Python 3, even if you have Python 2 also installed. 
