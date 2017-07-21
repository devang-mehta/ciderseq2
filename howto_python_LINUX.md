# How to install Python 3 and required modules on Linux


## Install Python 3

First check and update any existing Python installations:

```
sudo apt-get install
sudo apt-get update
sudo apt-get -y upgrade
```
Check which version of Python3 you have:
```
python3 -v
```

If you do not have python3, please install using:
```
sudo apt-get install python3
```

#Install pip

pip is a package manager for python. To install:

```
sudo apt-get install -y python3-pip
```

## Install modules

Open Terminal and type:

```
pip3 install [modulename]
```

for example:

```
pip3 install biopython
```
will install Biopython.

Please repeat and use pip3 to install all the required modules for ciderseq (see [Readme](README.md))

## Install MUSCLE
* Download Muscle from:www.drive5.com/muscle/downloads.htm
* Extract the file and copy to your cider-seq folder.
* Rename the file to 'muscle'
* If placed in another folder, make sure you edit the path address in the `ciderseq_config.json` file.

## Install BLAST
* Download the latest blast file from: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
* Extract the tar.gz file to any location on your machine:
```
tar -xzvf [blast filename].tar.gz -C [/target/directory]
```
* Export the BLAST executables to your path by using the following commands:
```
export PATH=$PATH:[/target/directory]
```
