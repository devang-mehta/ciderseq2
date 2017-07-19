# How to install PYTHON 3 and Required Modules on macOS X

A detailed guide on our favorite way to locally install Python and the required packages on your macOS X machine. 

## Install Homebrew

Homebrew is a package manager for macOS, much like `apt-get` on Linux.

To install, simply open your Terminal and type:
```
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
```

## Install Python 3

Once you have Homebrew installed, open Terminal and type:

```
brew install python3
```

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
