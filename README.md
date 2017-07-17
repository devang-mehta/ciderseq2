<H1>CIDER-Seq Data Analysis Software</H1>

This repository contains data analysis software for CIDER-Seq (Circular DNA Enrichment Sequencing) including an implementation of the **DeConcat** algorithm for sequence de-concatenation.

<H4>Reference:</H4>

Mehta D#, Hirsch-Hoffmann M, Patrignani A, Gruissem W, Vanderschuren H (2017) CIDER-Seq: unbiased virus enrichment and single-read, full length genome sequencing. ***bioRxiv***. doi:**insert bioRxiv link here**

<H2>System requirements</H2>

* Python 3
* [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE:BlastDocs&DOC_TYPE:Download "BLAST download page")
* [MUSCLE](http://www.drive5.com/muscle/ "MUSCLE homepage")

<H4>Python Modules</H4>

* [Biopython](https://www.biopython.org)
* [Click](http://click.pocoo.org/5/)
* [Matplotlib](https://matplotlib.org/)

<H4>Standard modules:</H4>

* sys
* os
* logging
* tempfile
* uuid
* json

<H2>Installation</H2>

<H2>Structure</H2>

The primary `ciderseq.py` script runs a pipeline on input sequence data consisting of the following process modules in order:
* **Separation** (`cider/separate.py`): Performs a BLASTn of input sequences against a user-defined dataset to bin reads into categories based on closed hit in the dataset. Non-hits are binned separately. This is useful in order to only process reads related to certain target sequences such as Virus A, Virus B and to filter out reads that do not match a target sequence(s). If switched off all the reads will be processed.
* **Align** (`cider/align.py`): An optional step which aligns reads against a reference sequence in order to trim sequence ends. Switch off if no reference sequence is available or if you wish to process all the reads in your input file (rather than just reads belonging to a specific target).
* **DeConcat** (`cider/deconcat.py`): The main de-contenation module. For a detailed description of the algorithm refer to *Mehta et al., 2017, bioRxiv*.
* **Annotation** (`cider/annotate.py`): This module annotates reads with gene/ORFs based on a reference dataset of proteins. It performs a tBLASTn to define the start and end sites of annotated genes/ORFs. It is capable of annotating circular sequences with genes crossing the sequence breakpoints.
* **Phasing** (`cider/phase.py`): Phasing ensures that all output sequences start at similar sites. This is essential for proper phylogenetic analysis. `phase.py` uses annotated sequences produced by `annotate.py` to phase genomes based on a user-defined `phaseto` gene. Hence, only annotated reads can be phased by this module.

*Note that `annotate.py` and `phase.py` are designed to deal with circular DNA sequences*

We also provide an optional `cider/cstools.py` script which assists with parallel computing and also generates performance plots.

<H2>Usage</H2>
The primary run command is:

`python ciderseq.py [options] [configfile] [inputfile]`

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

 `python ciderseq.py --format='fastq' /examples/ciderseq_config.json /examples/example1.fastq`

If we run:

`python ciderseq.py --format='fastq' --no-separation --no-alignment /examples/ciderseq_config.json /examples/example1.fastq`

the Separation and Alignment steps [described above](#structure) will be skipped.

Run `python ciderseq.py --help` for a brief description of usage and options.

<H3>Input Files</H3>

There are only two input files required:

* inputfile: file containing your raw sequence data in a Seq.IO compatible format (fasta, fastq, tab, gb)

* configfile: file containing editable parameters to customise your CIDER-Seq run. (see [below](#config-file) for details)


<H3>Output Files</H3>

<H2>Config File</H2>

see `examples/ciderseq_config.json` for an example.
We recommend using `examples/ciderseq_config.json` and editing the values for your analysis.

<H3>Explanation of name/value pairs:</H3>
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
			* "name of protein" : (protein strand information) "strand" : "1"=forward strand, "-1"=reverse strand"

 * phaseto : "name of protein to set sequence start position"
 * offset : (offset to start before protein position) if "10", phaseto protein will start at position 10

<H2>References</H2>
Hunter JD (2007) Matplotlib: A 2D graphics environment. Computing in Science & Engineering 9:3 90-95

doi:[10.5281/zenodo/573577](https://zenodo.org/record/573577#.WWjXm9N96L4)
