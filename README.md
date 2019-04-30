# CIDER-Seq Data Analysis Software

This repository contains data analysis software for CIDER-Seq (Circular DNA Enrichment Sequencing) including an implementation
of the **DeConcat** algorithm for sequence de-concatenation.

## Please cite:

Mehta D, Hirsch-Hoffmann M, Were M, Patrignani A, Were H, Gruissem W, Vanderschuren H (2019) A new full-length circular DNA sequencing method for viral-sized genomes reveals that RNAi transgenic plants provoke a shift in geminivirus populations in the field. ***Nucleic Acids Research*** , Volume 47, Issue 2; https://doi.org/10.1093/nar/gky914


CODE: [![DOI](https://zenodo.org/badge/94900933.svg)](https://zenodo.org/badge/latestdoi/94900933)


## Table of contents

  * [Prerequisites](#prerequisites)
    + [Python Modules](#python-modules)
    + [Standard Modules:](#standard-modules)
  * [Installation](#cider-seq-installation)
  * [Structure](#structure)
  * [Usage](#usage)
    + [For small datasets, or manual job handling](#for-small-datasets--or-manual-job-handling)
      - [Options](#options)
    + [For large datasets](#for-large-datasets)
      - [COMMAND](#command)
      - [Options](#options-1)
      - [Typical Usage](#typical-usage)
  * [Input Files](#input-files)
  * [Output Files](#output-files)
  * [Config File](#config-file)
    + [Explanation of name/value pairs](#explanation-of-name-value-pairs)
  * [License](#license)
  * [Disclamer](#disclamer)
  * [References](#references)


## Prerequisites
(All of the following are essential to run CIDER-Seq)

See [howto python macOS](howto_python_macOS.md) or [howto_python_LINUX.md](howto_python_LINUX.md) for a brief guide to installing the following:

* Python 3
* [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE:BlastDocs&DOC_TYPE:Download "BLAST download page")
* [MUSCLE](http://www.drive5.com/muscle/ "MUSCLE homepage")

### Python Modules

* [NumPy](http://www.numpy.org/)
* [Biopython](https://www.biopython.org)
* [Click](http://click.pocoo.org/5/)
* [Matplotlib](https://matplotlib.org/)

### Standard Modules:
(Should come pre-installed with Python)
* sys
* os
* logging
* tempfile
* uuid
* json

## CIDER-Seq Installation
The following are minimal steps to install the package after prerequisites are fulfilled:

1. Clone/Download this Git repository onto your machine.
2. Create a copy of configuration file `examples/ciderseq_config.json`.
3. Edit the configuration file and set all necessary parameters (see [Config File](#config-file))
4. Run `ciderseq` example (see [Usage](#usage))

## Structure

The primary `ciderseq.py` script runs a pipeline on input sequence data consisting of the following process modules in order:

* **Separation** (`cider/separate.py`): Performs a BLASTn of input sequences against a user-defined dataset to bin reads into categories based on closed hit in the dataset. Non-hits are binned separately. This is useful in order to only process reads related to certain target sequences such as Virus A, Virus B and to filter out reads that do not match a target sequence(s). If switched off all the reads will be processed.

* **Align** (`cider/align.py`): An optional step which aligns reads against a reference sequence in order to trim sequence ends. Switch off if no reference sequence is available or if you wish to process all the reads in your input file (rather than just reads belonging to a specific target).

* **DeConcat** (`cider/deconcat.py`): The main de-contenation module. For a detailed description of the algorithm refer to *Mehta et al., 2017, bioRxiv*.

* **Annotation** (`cider/annotate.py`): This module annotates reads with gene/ORFs based on a reference dataset of proteins. It performs a tBLASTn to define the start and end sites of annotated genes/ORFs. It is capable of annotating circular sequences with genes crossing the sequence breakpoints.

* **Phasing** (`cider/phase.py`): Phasing ensures that all output sequences start at similar sites. This is essential for proper phylogenetic analysis. `phase.py` uses annotated sequences produced by `annotate.py` to phase genomes based on a user-defined `phaseto` gene. Hence, only annotated reads can be phased by this module.

*Note that `annotate.py` and `phase.py` are designed to deal with circular DNA sequences. The modules cannot be executed directly but only via `ciderseq.py` script.*


## Usage

For regular users we recommend you ask your SMRT sequencing service provider to install CIDER-Seq on their computing cluster.
If you have access to your own cluster or wish to process only a small dataset read on:


### For small datasets, or manual job handling

The primary run command is:

```
python3 ciderseq.py [options] CONFIGFILE INPUTFILE
```

#### Options

The options allow the user to skip certain processes in the `ciderseq.py` pipeline:

```
--no-separation
--no-alignment
--no-deconcatenation
--no-annotation
--no-phasing
```

Additionally, the option `--format` allows you to choose the format of your input file from either fasta, fastq, tab or gb. The default option is fasta.

 For example, for processing the files in `examples/` we run:

 ```
 python3 ciderseq.py --format fastq examples/ciderseq_config.json examples/example1.fastq
 ```

If we run:

```
python3 ciderseq.py --format fastq --no-separation --no-alignment examples/ciderseq_config.json examples/example1.fastq
```

the Separation and Alignment steps [described above](#structure) will be skipped.

Run `python3 ciderseq.py --help` for a brief description of usage and options.

### For large datasets

We provide `cstools.py` in order to process large sequence datasets on a computing cluster (similar to the one your SMRT Analysis software
is installed on). `cstools.py` is basically a file-handling program which will allow you to `split` your input sequences into batches
and then `join` the outputs into a single results directory. Command `plot` will create a set of different plots based on the processed 
results, separation- and phase-results are required. 

To run:

```
python3 cstools.py COMMAND [options] CONFIGFILE INPUTFILE
```

The `CONFIGFILE` and `INPUTFILE` are the same as the ones used by `ciderseq.py`. See [below](#input-files) for details.

#### COMMAND

`cstools.py` has two primary actions: `split` and `join` which split input files into several jobs and then join the split jobs once
processed. A third action called `plot` can be run after processing to produce DeConcat statistics. See `examples/plots` for examples of each chart.


#### Options

```
--format
--numjobs
--cluster
--clean
```

> `format`: is the format of the input file, default FASTA

>`numjobs`: is the number of processing jobs to create. the sequences equally distributed to the provides number of jobs. The default is 1.

>`cluster`: contains cluster submission parameters, e.g. `"bsub -n 4"` in a LSF environment.

>`clean`: if absent, the folders with the split data will not be deleted. Useful for debugging.


#### Typical Usage

The typical run command is:

```
python3 cstools.py split --format fastq --numjobs 8 --cluster "bsub -n 1" examples/ciderseq_config.json examples/example1.fastq
```

>The `split` command outputs all necessary execution commands and ask confirmation before executing them.


followed by (for joining):
```
python3 cstools.py join --clean examples/ciderseq_config.json examples/example1.fastq
```

>The `clean` option will remove the folder and files crated during **split**.

And (for plotting):
```
python3 cstools.py plot examples/ciderseq_config.json examples/example1.fastq
```

## Input Files

There are only two input files required:

* `CONFIGFILE`: file containing editable parameters to customise your CIDER-Seq run. (see [below](#config-file) for details)

* `INPUTFILE`: file containing your raw sequence data in a Seq.IO compatible format (fasta, fastq, tab, gb)

## Output Files

All output files are saved in the directory specified in the `CONFIGFILE`.
We recommend you create separate folders for each step defined in the `CONFIGFILE`.

See our example file `examples/ciderseq_config.json` for details.

## Config File

`ciderseq.py` reads configuration parameters from the provided `CONFIGFILE`.
see `examples/ciderseq_config.json` for an example.

Please read on for a detailed explanation of how to edit the `CONFIGFILE` file:

**We recommend using `examples/ciderseq_config.json` and editing the values for your analysis.**

There are several options available to modify in the `CONFIGFILE`.

Essential changes that you will almost certainly need to make are **highlighted in green** in the images below.

**Overall Settings**


![alt-text][image1]


**Configuring `separate.py`:**


![alt-text][image2]


**Configuring `align.py`:**


![alt-text][image3]


**Configuring `deconcat.py`:**


![alt-text][image4]


**Configuring `annotate.py`:**


![alt-text][image5]


**Configuring `phase.py`:**


![alt-text][image6]


**Configuring Phasing:**


![alt-text][image7]


### A more detailed explanation of name/value pairs

* loglevel : "DEBUG", will write log information into output dir
* outputdir : "destination of log information"
* separate (settings used in `cider/separate.py`)
  * outputdir : "output directory of results"
  * blastinit : "python code to initialize blast environment, if necessary"
  * blastexe : "blast executable"
  * blastndb : "location of the indexed blast database"
  * evalue : "used evalue for blastn", suggested="1"

* align (settings used in `cider/align.py`)
  * outputdir : "output directory of results"
  * muscleinit : "python code to initialize muscle environment, if necessary"
  * muscleexe : "muscle executable"
  * targets : (array of align target DNA fasta files)
     * "name of genome" : "path to DNA fasta file"
  * windowsize : "size of the sliding window", suggested="10"

* deconcat (settings used in `cider/deconcat.py`)
  * outputdir : "output directory of results"
  * muscleinit : "python code to initialize muscle environment, if necessary"
  * muscleexe : "muscle executable"
  * fragmentsize : "size of de-concatenation fragment", recommended="30"
  * statistics : (write statistic output) ("1"=on, "0"=off)

* annotate (settings used in `cider/annotate.py`)
  * outputdir : "output directory of results"
  * blastinit : "python code to initialize blast environment, if necessary"
  * blastexe : "blast executable"
  * tblastndb : "location of the indexed blast database"
  * evalue : "used evalue for tblastn", suggested="0.01"

* phase (settings used in `cider/phase.py`)
  * outputdir : "output directory of results"
  * outputformat : (format of results as array, use SeqIO valid formats) e.g. ["genbank","fasta"]=outputs both .gb and .fasta formats
  * phasegenomes : (array of genomes for phasing)
     * "name of genome" : (array of phasing parameters)
         * proteins : array of proteins in genome
             * "name of protein" : (protein strand information)
                 * "strand" : "1"=forward strand, "-1"=reverse strand"  
         * phaseto : "name of protein to set sequence start position"
         * offset : (offset to start before protein position) if "10", phaseto protein will start at position 10

## License
Copyright &copy; 2017, Matthias Hirsch-Hoffmann and Devang Mehta, ETH Zurich. CIDER-Seq is licensed under the GNU General Public License, version 3 (https://www.gnu.org/licenses/gpl-3.0.en.html) or see [LICENSE.txt](LICENSE.txt).

## Disclamer
This software is supplied 'as is' without any warranty or guarantee of support. ETH Zurich is not responsible for its use, misuse, or functionality. In no event shall the authors or copyright holders be liable for any claim, damages, or other liability arising from, out of, or in connection with this software.

## References
Hunter JD (2007) Matplotlib: A 2D graphics environment. Computing in Science & Engineering 9:3 90-95
doi:[10.5281/zenodo/573577](https://zenodo.org/record/573577#.WWjXm9N96L4)

[image1]: https://github.com/hirschhm/ciderseq/blob/master/config-images/ciderseq_ciderseq_config_json_at_master_%C2%B7_hirschhm_ciderseq.png
[image2]: https://github.com/hirschhm/ciderseq/blob/master/config-images/ciderseq_ciderseq_config_json_at_master_%C2%B7_hirschhm_ciderseq%202.png
[image3]: https://github.com/hirschhm/ciderseq/blob/master/config-images/ciderseq_ciderseq_config_json_at_master_%C2%B7_hirschhm_ciderseq%203.png
[image4]: https://github.com/hirschhm/ciderseq/blob/master/config-images/ciderseq_ciderseq_config_json_at_master_%C2%B7_hirschhm_ciderseq%204.png
[image5]: https://github.com/hirschhm/ciderseq/blob/master/config-images/ciderseq_ciderseq_config_json_at_master_%C2%B7_hirschhm_ciderseq%205.png
[image6]: https://github.com/hirschhm/ciderseq/blob/master/config-images/ciderseq_ciderseq_config_json_at_master_%C2%B7_hirschhm_ciderseq%206.png
[image7]: https://github.com/hirschhm/ciderseq/blob/master/config-images/ciderseq_ciderseq_config_json_at_master_%C2%B7_hirschhm_ciderseq%207.png
