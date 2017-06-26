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

</ul>

<H2>Config File</H2>

<p>see examples/ciderseq_config.json<br><br>
explanation of name/value pairs:<br>
<ul>
<li>loglevel = DEBUG, will write log information into output dir</li>
<li>outputdir = destination of log information</li>
<li>separate</li>
<ul>
	<li>outputdir = output directory of results</li>
	<li>blastinit = python code to initialize blast if necessary</li>
	<li>blastexe = blast executable </li>
	<li>blastndb = location of the indexed blast database</li>
	<li>evalue = used evalue for blastn</li>
</ul>
</ul>

</p>

<H2>Examples</H2>


