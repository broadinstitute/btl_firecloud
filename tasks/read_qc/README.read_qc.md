# README for read_qc

 To run:

1) Create a tar.gz file containing the genome reference.  Ensure there is no directory structure inside the tarfile, and that there is only a single fasta file.  To reduce bulk, ensure there are no core files or BadCoverage.cache directories.
eg 
cd /location/of/genome/reference
tar cvfz my_genome_build_v2.tar .
gzip --best my_genome_build_v2.tar

2) upload the fasta files and ref data to the bucket
eg
gsutil cp * gs://fc-a4c39031-98b7-4e1f-9937-eec981330527/test_data
*** Directory naming conventions?

3) Create a tsv loadfile containing the metadata.  Include the paths from the previous step
Example lives in test_data/test_data.read_qc.participant.loadfile.txt

4) Load the metadata into Firecloud
Go to the workspace, go to the data tab, and press the import metadata button.  
Press import from file
Choose the file you just made.
Press Upload

5) Run the workflow
Go to the Workspace, go to the Method Configurations tab, select the read_qc link
Press the Launch analysis button
Select the sample you want to run
Press the Launch button

6) Monitor the workflow
Go to the workspace, go to the Monitor tab
Drill down into a specific submission if you are interested

7a) Check results - by job
Go to the workspace, go to the monitor tab, select the view link next to the job you want to check
Continue from there with a lot of clicking...

7b) Check results - by sample
Go to the workspace, go to the data tab, ensure Participants is selected
Click on the links in the columns corresponding to the outputs you want to view.
With the Open link in the modal dialog that pops up, either rightclick-save link as... or rightclick-open new tab.

7c) Check results - via commandline
*** choose output file naming convention to make this straightforward
gsutil cp <bucket paths> <local path>
gsutil rsync -r <bucket paths> <local path>
 
 
 -----------
 
        $0  [-h] [-s <parent_ssf_ticket> ] [-l <labset_ssf_ticket> ] [-i <input_table>] [-t  <data_type> ] [-N] [-O] 


        Required:  

	-s <parent_ssf_ticket>		Parent SSF ticket. Examples include SSF-1002, SSFDEV-300, BTLDEV-195

	-l <labset_ssf_ticket>		Labset SSF ticket. Examples include SSF-2270, SSFDEV-301, BTLDEV-196

	-t <data_type>			Read data type:  pcr-free|jump|low-input

	-i <input_table>		Tab-separated input file with the following columns (showing header, and example):
					A:	sample_id	SM-G8KZ1	(required)
					B:	genus		Escherichia	(required)
					C:	species		coli		(required)
					D:	strain		K12 MG1655	(required)
					E:	specimen_id	SM-G8KZ1	(required only if exists and differs from strain)
					F:	gnumber		G90000
					G:	ref_path	/seq/ref/Ecoli/Ecoli.fasta	
					H:	fastq_r1	/path/to/fastq_r1
					I:	fastq_r2	/path/to/fastq_r2
					J:	bam		/path/to/bam	required if not fastq given
	
	Optional:

	-h		           	Print this help message.

	-N				No run.  Only print the commands to be run.

	-O				Overwrite.  Run from beginning even if intermediate/final files exist.
