# ciderseq
CIDERSEQ - Circular DNA Enrichment Sequencing - de-concatenation package

<H2>System requirements</H2>

<ul>
	<li>Python 3</li>
	<li>Blast - <a href="https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download" target=_blank>https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download</a> </l>
	<li>Muscle - <a href="http://www.drive5.com/muscle/" target=_blank>http://www.drive5.com/muscle/</a> </li>
</ul>

<H2>Python Modules</H2>
<ul>
	<li>Biopython - </li>
	<li>Click - </li>
	<li>Standard modules:</li>
	<ul>
		<li>sys</li>
		<li>os</li>
		<li>logging</li>
		<li>tempfile</li>
		<li>uuid</li>
		<li>json</li>
	</ul>
</ul>

<H2>Config File</H2>

<p>see examples/ciderseq_config.json<br><br>
explanation of name/value pairs:<br>
<ul>
<li>loglevel = DEBUG, will write log information into output dir</li>
<li>outputdir = destination of log information</li>

<li>separate = array of subvalues (s. following)</li>
<ul>
	<li>outputdir = output directory of results</li>
	<li>blastinit = python code to initialize blast environment, if necessary</li>
	<li>blastexe = blast executable </li>
	<li>blastndb = location of the indexed blast database</li>
	<li>evalue = used evalue for blastn</li>
</ul>

<li>align = array of subvalues (s. following)</li>
<ul>
	<li>outputdir = output directory of results</li>
	<li>muscleinit = python code to initialize muscle environment, if necessary</li>
	<li>muscleexe = muscle executable </li>
	<li>targets = array of align target DNA fasta files (s. following example)</li>
	<ul>
		<li><name of genome> = location and name of genome DNA fasta file</li> 
	</ul>
	<li>windowsize = used evalue for blastn</li>
</ul>

<li>deconcat = array of subvalues (s. following)</li>
<ul>
	<li>outputdir = output directory of results</li>
	<li>muscleinit = python code to initialize muscle environment, if necessary</li>
	<li>muscleexe = muscle executable </li>
	<li>fragmentsize = 	size of deconcatenation fragment, e.g. 30</li>
	<li>statistics = write statistic output (1 on, 0 off)</li>
</ul>
	
<li>annotate = array of subvalues (s. following)</li>
<ul>
	<li>outputdir = output directory of results</li>
	<li>blastinit = python code to initialize blast environment, if necessary</li>
	<li>blastexe = blast executable </li>
	<li>tblastndb = location of the indexed blast database</li>
	<li>evalue = used evalue for tblastn</li>
</ul>
	
<li>phase = array of subvalues (s. following)</li>
<ul>
	<li>outputdir = output directory of results</li>
	<li>outputformat = format of results as array, valid format s. SeqIO, e.g. ["genbank","fasta"] </li>
	<li>phasegenomes = array of genomes for phasing (s. following)</li>
	<ul>
		<li><name of genome> = array of phasing parameters (s. following)</li>
		<ul>
			<li>proteins = array of proteins in genome (s. following)</li>
			<ul>
				<li><name of protein = array of protein details (s. following)</li>
				<ul>
					<li>stand = 1 : forward strand, -1 reverse strand</li>
				</ul>
				<li>phaseto = name of protein to set 0 position</li>
				<li>offset = offset to start before protein position, e.g. 10, phaseto protein will start at position 10</li>
			</ul>
		</ul>
	</ul>
</ul>
</ul>
</p>

<H2>Examples</H2>


