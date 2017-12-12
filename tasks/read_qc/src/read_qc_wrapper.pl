#!/usr/bin/env perl 
use strict;

#use Getopt::Long;
use Getopt::Long qw(:config no_ignore_case);
use vars qw/ %opt /;
use Tie::IxHash;
#use FindBin;
#FindBin->again;
#use lib "$FindBin::Bin/../lib/tools";
use lib "/gsap/assembly_analysis/tools/lib/";
use AAGnuplot;						# AAGnuplot, LaTeX_aa, scalar_manip_aa to generate plots 
use LaTeX_aa;
use scalar_manip_aa;

sub usage() {
    print <<HELP;

        $0  [-h] [-b <bam> ] [-fq1 <read1_fq> ] [-fq2 <read2_fq>] [-l <list_file> ] [-a <adapter_fasta>] [-s <read_subset_count>] -t <num_threads> ]  [-R <rRNA_reference_file>] [-A <aligned_bam_reference] [-z] [-r] [-C] [-S] [-B] [-N] [-T] [-P] 


        Required:  

	-b|--bam		Input BAM file.  Filename must contain "bam".  SAM format is not supported.

	OR

	-fq1			Read 1 fastq file.  Filename must contain either "fastq" or "fq".
	-fq2 			(optional) Read 2 fastq file.  Filename must contain either "fastq" or "fq".

	OR	

	-l|--list_input		List of input BAM or fastq files.  One BAM or fastq (or set of fastq files) per line.
				Paired f1 and f2 fastq filenames should be space-separated.  Optionally the user may add a label to 				    each line to describe the data.  After the bam or fastq filenames there should be a space then the 
				label. Example:  

				/seq/picard/A0U19/C1-210_2012-05-08_2014-01-28/1/Solexa-93208/A0U19.1.bam	kapa_frags 
				/seq/picard/A0U44/C1-308_2012-05-10_2014-01-28/1/Solexa-93523/A0U44.1.bam	pcr-free
				/btl/projects/terry/cheapseq.r1.fastq  /btl/projects/terry/cheapseq.r2.fastq	cheapseq
	
	Optional:

	-h|help           	Print this help message.
	-s| --subset		Number of reads to subsample for read blast to nt [ default:  -s 10000 ] 
	-a| --adapter           Fasta file of adapter sequences to identify in reads.  By default 6 adapters sequences searched.
	-o| --out		Output PDF filename [default="read_qc.pdf"]
	-t| --threads           Number of threads for read blast [ default = 8 ]
	-A| --aligned_ref	Reference file used in aligned BAM input [ suggested if aligned BAM header without reference ] 
	-R| --rrna_ref		rRNA reference file [example:  /cil/shed/sandboxes/mbusby/SSF/scripts/QuanitfyRRNA/References/Human/rRNA.fasta ] 		
	-z| --lorenz            Run the Lorenz curve module. [ default = do not run Lorenz curve module ] 
	-r| --read_ln_histo     Run the read length histogram. [ default = do not run read length histogram ]
	-B| --no_read_blast     Do not run the read blast module [ default = run read blast module ]
	-C| --no_complexity	Do not run library complexity module.  [default = run complexity module]
	-J| --no_jump_eval	Do not run picard CollectJumpingLibraryMetrics.  [default = run odule]
	-S| --no_coverage       Do not run the coverage module. [default = run coverage module for aligned BAM input ] 
	-N| --no_overwrite      Do not overwrite existing files.  Use existing files for analysis.
	-T|--no_tmp_dir		Do not use tmp dir space in /broad/hptmp/<user> to write intermediate files.  Write output in cwd.
	-P|--no_plots		Do not plot charts, print only metrics.  Useful for many input lines (P1) [ default = plot charts ]


HELP
    exit;
}

if ($opt{h}) { usage(); exit; }

my ($input_list, $input_bam, $input_fq1, $input_fq2, $rrna_ref, $adapter_file, $help, $no_tmp_space, $aligned_ref);
my ($no_complexity, $no_coverage, $no_blast, $run_lorenz, $run_read_ln_histo, $no_plots, $no_overwrite, $no_jump);
my ($out_prefix) = ("read_qc"); 
my ($threads, $subset) = (8, 10000);
#my $short_blast_rd_ln 50;					# found that blastn-short on 10,000 reads takes too long, turn off

GetOptions(
	'b|bam:s'	   => \$input_bam,
	'fq1:s'		   => \$input_fq1,
	'fq2:s'		   => \$input_fq2,
	'l|list_input:s'   => \$input_list,
	's|subset'         => \$subset,
	't|threads'        => \$threads,
	'a|adapter:s'      => \$adapter_file,
	'R|rrna_ref:s'	   => \$rrna_ref,
	'A|aligned_ref:s'  => \$aligned_ref,
	'o|out:s'          => \$out_prefix,   
	'J|no_jump_eval'   => \$no_jump,
	'C|no_complexity'  => \$no_complexity,
	'S|no_coverage'    => \$no_coverage,
	'P|no_plots'	   => \$no_plots,
	'N|no_overwrite'   => \$no_overwrite,
	'z|lorenz'         => \$run_lorenz,
	'r|read_ln_histo'  => \$run_read_ln_histo,  
	'B|no_blast'       => \$no_blast,
	'T|tmp_space'      => \$no_tmp_space,
	'h|help'	   => \$help,
	);

usage() and exit if ($help) ;

##############################################        Tools	  ############################################

my $bwa                 = "/seq/software/picard/current/3rd_party/bwa/bwa";
my $fastqc              = "/cil/shed/apps/external/LibraryQC/FastQC/fastqc";
my $samtools            = "/broad/software/groups/gtba/software/samtools_0.1.18/bin/samtools";
#my $picard              = "/seq/software/picard/current/bin/";
my $picard              = "/seq/software/picard-public/current/";
my $aa_picard           = "/gsap/assembly_analysis/bin/picard/dist/"; 
my $blastn              = "/broad/software/groups/gtba/software/ncbi-blast-2.2.25+/bin/blastn";
my $tax_from_blast      = "/gsap/assembly_analysis/bin/taxonomy_from_blast.pl"; 
my $lorenz_tool         = "/home/unix/tshea/dev/lorenz_plot.py";		# modifcation of aaron tool
#my $complexity_by_start = "/cil/shed/apps/internal/bam_utilities/ComplexityByStartPos/ComplexityByStartPos";         
#my $complexity_by_start = "/cil/shed/apps/internal/RNA_utilities/ComplexityByStartPos/ComplexityByStartPos";

##############################################      Variables     ############################################

my @pdf_info_type = ("read group number", "read group name", "label", "total reads (r1)", "total reads (r2)", "passing filter (PF) reads", "% PF reads", "mean read length (r1)", "mean read length (r2)", "mean base quality (r1)", "mean base quality (r2)", "% Q2 bases (r1)", "% Q2 bases (r2)", "% Q10- bases (r1)", "% Q10- bases (r2)", "% Q20+ bases (r1)", "% Q20+ bases (r2)", "% Q30+ bases (r1)", "% Q30+ bases (r2)", "mean % GC (r1)", "mean % GC (r2)", "% duplication (de novo)", "estimated library size", "% rRNA", "% adapter (pair)", "% adapter (r1)", "% adapter (r2)", "percent aligned (pair)", "percent aligned (r1)", "percent aligned (r2)", "mean insert size (FR)", "median insert size (FR)", "mean insert size (RF)", "median insert size (RF)", "mean aligned coverage", "mean jump size", "percent jumps", "percent non-jumps", "percent chimera");
my @info_type     = qw(short_id label r1_total_reads r2_total_reads r1_pf_reads  pct_pf r1_mean_rd_ln r2_mean_rd_ln r1_mean_qual r2_mean_qual r1_pct_q2 r2_pct_q2 r1_pct_q10 r2_pct_q10 r1_pct_q20 r2_pct_q20 r1_pct_q30 r2_pct_q30 r1_mean_gc r2_mean_gc pct_duplication estimate_library_size pct_rrna pct_adapter r1_pct_adapter r2_pct_adapter pct_poly_a pct_poly_t pct_aligned r1_pct_aligned r2_pct_aligned fr_mean_insert fr_median_insert rf_mean_insert rf_median_insert mean_cvg jump_mean pct_jump pct_nonjump pct_chimera);
my @pdf_data;
my @toc_data;								# store read group number, full ID, short ID
my @toc_header = ("read input number", "read input name short", "label", "read input name full");

tie my %input, "Tie::IxHash"; 
my %orgs_rd_blast;
my %org_count; 							# hash of orgs in read blast, val is total # reads (for order)	
my @org_list;							
my @read_blast_data;
my @read_blast_headers;

my $pdf_file  = $out_prefix  . ".pdf";
my $stats_txt = $out_prefix . ".metrics.txt";

my %qual_by_cycle_r1;							# store mean base quality across read length (r1)
my %qual_by_cycle_r2;							# store mean base quality across read length (r2)
my ($max_read_ln_r1, $max_read_ln_r2) = (0,0);
my %qual_histo_r1;							# store number of bases at each Q score (r1)
my %qual_histo_r2;							# store number of bases at each Q score (r2)
my $max_qual_score = 0;							# highest Q score seen across any data sets in input

my %rd_ln_histo_r1;							# store the number of reads at binned rd ln (r1)
my %rd_ln_histo_r2;							# store the number of reads at binned rd ln (r2)
my %read_ln_histo_bins;							# store the bins used for rd ln histo

my %gc_histo_r1;							# store the number of reads at each %GC (r1)
my %gc_histo_r2;							# store the number of reads at each %GC (r2)

my $plot_insert_size = 0;						# by default no plot until an aligned data set found	
my %insert_size_histo_rf;						# store the # of pairs at given insert size (RF orientation)
my %insert_size_histo_fr;						# store the # of pairs at given insert size (FR orientation)
my $max_fr_size_to_plot = -1;
my $max_rf_size_to_plot = -1;

my $plot_complexity = 0;						# by default no plot until complexity data found
#my $plot_complexity_by_start = 0;					# by default no plot until complexity data found
my %complexity_histo;							# histogram of number of pairs seen # of times
#my %complexity_by_start_histo;						# to store the complexity by stat rarefactionData.txt 
my $max_complexity_count = 0;

my %coverage_plot;							# store the coverage at each window in reference
my %coverage_histo;							# store the coverage at each window in reference
my $max_histo_coverage = 0;
my %gc_bias_histo;							# store the normalized coverage at each GC value
my %gc_bias_ref_windows;						# store the number of ref windows at each GC value
my %references;								# structure to store info about ref-based plots
my %window_hash;

my $adapter_string;

my $user = `whoami`;
chomp($user);
my $qc_main_subdir;							# subdirectory for all files created by script 
if ($no_tmp_space) {
	$qc_main_subdir = "read_qc"; 				
} else {
	$qc_main_subdir = "/broad/hptmp/$user";
}
my $qc_fastqc_subdir       = "$qc_main_subdir/fastqc_data";		# create separate subdirectory for all Fastqc output
my $qc_read_data_subdir    = "$qc_main_subdir/read_data";
my $qc_metric_data_subdir  = "$qc_main_subdir/metrics_data";
my $qc_chart_data_subdir   = "$qc_main_subdir/chart_data"; 

my ($qual_by_cycle_r1_chartobj, $qual_by_cycle_r2_chartobj, $qual_histo_r1_chartobj, $qual_histo_r2_chartobj);
#my ($mean_gc_r1_chartobj, $mean_gc_r2_chartobj, $complexity_chartobj, @complexity_by_start_chartobj);
my ($mean_gc_r1_chartobj, $mean_gc_r2_chartobj, $complexity_chartobj);
my ($insert_size_histo_fr_chartobj, $insert_size_histo_rf_chartobj, $coverage_histo_chartobj, @coverage_chartobj);
my ($rd_ln_histo_r1_chartobj, $rd_ln_histo_r2_chartobj);
my (@gc_bias_chartobj);

##############################################        Main        ############################################

&checkUsage;							# subroutine to verify input provided correctly
&loadDotkits;
&verifyBamOrFqInput;						# subroutine to get data input from --bam or --fq1/--fq2 
&getListInput if ($input_list);
&evaluatePairing;
&makeOutputDirs;
&makeSeparateBams;

$adapter_string = &buildAdapterCommand($adapter_file) if ($adapter_file);

&runFastQcAndPicardModules;					# Runs FastQC and Picard modules (MeanQualityByCycle, etc..)

&getPercentPfFromFastq;						# with %PF not reported in most bams try to get from fastqs

&getPercentAlignedByRead;					# for paired+aligned input get % aligned by r1 and r2	

&getQualityYieldMetrics;

&getQualityByPosition;
if (! $no_plots) {
	($qual_by_cycle_r1_chartobj) = &plotQualityByPosition("read1", \%qual_by_cycle_r1, $max_read_ln_r1);
	($qual_by_cycle_r2_chartobj) = &plotQualityByPosition("read2", \%qual_by_cycle_r2, $max_read_ln_r2) if ($max_read_ln_r2 > 0);
}

&getQualityHisto;
if (! $no_plots) {
	($qual_histo_r1_chartobj) = &plotQualityHisto("read1", \%qual_histo_r1, $max_qual_score) ;
	($qual_histo_r2_chartobj) = &plotQualityHisto("read2", \%qual_histo_r2, $max_qual_score) ;
}

&getGcHisto;
if (! $no_plots) {
	($mean_gc_r1_chartobj) = &plotGcHisto("read1", \%gc_histo_r1);
	($mean_gc_r2_chartobj) = &plotGcHisto("read2", \%gc_histo_r2) if ($max_read_ln_r2 > 0);
}

&getAdapterMetrics;

&getInsertSizeMetrics;

&rRnaAlign($rrna_ref) if ($rrna_ref) ;

($max_fr_size_to_plot, $max_rf_size_to_plot) = &getInsertSizeHisto;

if ($plot_insert_size == 1 && ! $no_plots) {
	($insert_size_histo_rf_chartobj) = &plotInsertSizeHisto("RF", \%insert_size_histo_rf, $max_rf_size_to_plot);
	($insert_size_histo_fr_chartobj) = &plotInsertSizeHisto("FR", \%insert_size_histo_fr, $max_fr_size_to_plot);
}

unless ($no_complexity) {
	&getComplexityMetrics;
}
unless ($no_jump) {
	&getJumpMetrics;
}

unless ($no_coverage) {
	&getCoverageMetrics;
	if (! $no_plots) {
		($coverage_histo_chartobj) = &plotCoverageHisto(\%coverage_histo,$max_histo_coverage);
		(@coverage_chartobj) = &plotCoverageAlongRef(\%coverage_plot);
		$no_coverage = 1 if (scalar @coverage_chartobj < 1 ) ;
		&getGcBiasHisto;
		(@gc_bias_chartobj) = &plotGcBiasHisto(\%gc_bias_histo, \%gc_bias_ref_windows);
		&runLorenzCurve if ($run_lorenz);
	}
}	

unless ($no_blast) {
	&readBlast;
	&prepareBlastTable;
}

if ($run_read_ln_histo) {
	&getReadLnHisto;
	if (! $no_plots) {
		($rd_ln_histo_r1_chartobj) = &plotReadLnHisto("read1", \%rd_ln_histo_r1);
		($rd_ln_histo_r2_chartobj) = &plotReadLnHisto("read2", \%rd_ln_histo_r2) if ($max_read_ln_r2 > 0);
	}	
}
&prepareSummaryTable;
&createPdf;

##############################################    Subroutines     ############################################

sub checkUsage {
	if ($input_fq2 && ! $input_fq1) {
		print "\nThe --fq2 (read2) option requires that --fq1 (read1) also be used.\n\n";
		usage() and exit;
	}
	if (! $input_bam && ! $input_fq1 && ! $input_list) {
		print "\nProvide either -b <bam_file> or -fq1 <fastq_file> or -l <list_input_file>\n\n";
		usage() and exit;
	}
	if (($input_bam ||$input_fq1)  && $input_list) {				# if user gives too much input
		print "\nThe -l list option should not be used with -fq1/-fq2 or -b\n\n";
		usage() and exit;
	}
	if ($input_bam && $input_fq1) {
		print "\nThe --bam and --fq1 options may not be used together.  Please select BAM or fastq as input type\n\n";
		usage() and exit;
	}
}

sub loadDotkits {
	my $tmp_file = "run_dk_cmds.tmp";
	open(DK,">$tmp_file") || die "Cannot write $tmp_file file:  $!\n\n";

	print DK "eval `/broad/tools/dotkit/init`\n";
	print DK "source /broad/software/scripts/useuse\n";
	print DK "reuse group=gtba\n" if ($run_lorenz);			# only load next 3 dk if lorenz_curve option set
	print DK "reuse GAEMR\n"  if ($run_lorenz);
	print DK "reuse .matplotlib-1.1.1rc-python-2.7.1-sqlite3-rtrees\n"  if ($run_lorenz);
	print DK "reuse .gnuplot-4.4.0\n";
	print DK "reuse R-2.15\n";
	print DK "reuse Java-1.8\n";
	print DK "reuse Perl-5.10\n";
	print DK "reuse BLAST+" unless ($no_blast);
	close DK;
	print "Loading dotkits...\n\n";
	system "chmod 755 $tmp_file";
	system "./$tmp_file";
	#unlink $tmp_file;
}

sub verifyBamOrFqInput {
	my ($id);
	if ($input_bam) {
		if (! -e $input_bam) {							# if BAM file is not found
			print "\nCannot locate $input_bam BAM file\n\n";
			usage() and exit;
		}
		if ($input_bam !~ /bam/) {						# if BAM file mis-named
			print "\nCheck $input_bam BAM filename to verify it includes \"bam\" in filename\n\n";
			usage() and exit;
		}
		$id = $input_bam;
		$id =~ s/.*\///;			# strip off path to filename
		$id =~ s/\.bam$//;			# strip .bam extension off the ID
		$input{$id}{'bam'} = $input_bam;
		$input{$id}{'id'}  = $id;
	} elsif ($input_fq1) {
		if (! -e $input_fq1) { 							# if fq file not found
			print "\nCannot locate $input_fq1 fastq file\n\n";
			usage() and exit;
		}
		if ($input_fq1 !~ /fastq|fq/) {						# if fq file mis-named
			print "\nCheck $input_fq1 fastq filename to verify it includes \"fastq\" \"fq\" in filename\n\n";
			usage() and exit;
		}
		$id = $input_fq1;
		$id =~ s/.*\///; 			# strip off path to filename
		$id =~ s/\.fastq|\.fq$//;		# strip extension off the ID
		$id =~ s/\.r*\d$//i if ($input_fq2);	# for paired fq1 and fq2 files, strip off any "r1" or "R1" suffix
		$input{$id}{'fq1'} = $input_fq1;
		$input{$id}{'id'}  = $id;
		$input{$id}{'aligned'} = "false";
		$input{$id}{'pct_aligned'} = "na";
		if ($input_fq2) {
			if (! -e $input_fq2) { 							# if fq file not found
				print "\nCannot locate $input_fq2 fastq file\n\n";
				usage() and exit;
			}
			if ($input_fq2 !~ /fastq|fq/) {						# if fq file mis-named
				print "\nCheck $input_fq2 fastq filename to verify it includes \"fastq\" \"fq\" in filename\n\n";
				usage() and exit;
			}
			$id = $input_fq2;
			$id =~ s/.*\///; 			# strip off path to filename
			$id =~ s/\.fastq|\.fq//;		# strip extension off the ID
			$id =~ s/\.r*\d$//i;			# strip off any "r2" or "R2" suffix
			$input{$id}{'fq2'} = $input_fq2;
			$input{$id}{'paired'} = "true";
		} else {
			$input{$id}{'paired'} = "false";
		}
	} else {
		# list file will be evaluated in getListInput subroutine
	}
}
sub getListInput {
	my ($id, $bam, $unpaired_fq, $fq1, $fq2);
	my $label = " - ";
	open(IN,$input_list) || die "Cannot open the $input_list input list file:  $!\n\n";
	while(<IN>){
		my $line = $_;
		chomp($line);
		my @line = split(/\s+/,$line);
		if ($line =~ /(\S+)\s+(\S+)/ && $line !~ /bam/) {	# input in list is r1 and r2 fastq and there is a tag
			($fq1, $fq2) = ($1, $2);
			$label = $1 if ($line =~ /\S+\s+\S+\s+(\S+)/); 
			if ($fq1 !~ /fastq|fq/ || $fq2 !~ /fastq|fq/) {
				print "\nCheck $input_list to confirm that the two fastq files are named appropriately\n\n";
				usage() and exit;
			}
			die "Cannot locate $fq1 fq1 fastq file\n\n" if (! -e $fq1); 	# if fq1 file is not found
			die "Cannot locate $fq2 fq2 fastq file\n\n" if (! -e $fq2); 	# if fq2 file is not found
			$id = $fq1;
			$id =~ s/.*\///; 			# strip off path to filename
			$id =~ s/\.fastq|\.fq//;		# strip extension off the ID
			$id =~ s/\.r*\d$//i;                     # strip off any "r2" or "R2" suffix
			$input{$id}{'fq1'}     = $fq1;
			$input{$id}{'fq2'}     = $fq2;
			$input{$id}{'paired'}  = "true";
			$input{$id}{'aligned'} = "false";
			$input{$id}{'id'}      = $id;
			$input{$id}{'label'}   = $label;
		} elsif ($line =~ /bam/) {
			$label = $1 if ($line =~ /\S+\s+(\S+)/); 
			$bam = $line[0];
			$bam =~ s/\s+//g;
			die "Cannot locate $bam BAM file\n\n" if (! -e $bam); 	# if BAM file is not found
			$id = $bam;
			$id =~ s/.*\///;			# strip off path to filename
			$id =~ s/\.bam$//;			# strip .bam extension off the ID

			####  Do a check to see if multiple BAM files in list have same name. If yes, add label to ID  ###
			$id = $id . "_" . $label if (exists $input{$id}{'bam'} && $label ne " - "); 
			if (exists $input{$id}{'bam'}) {
				die "\nFatal error in $input_list:\n\nCombination of BAM file and label not unique.  Either rename BAM files to be unique or else choose uniqe label names for identically-named BAM files\n\n";
			}
			####################  Proceed if now have unique identifier for hash key  ########################

			$input{$id}{'bam'}   = $bam;
			$input{$id}{'id'}    = $id;
			$input{$id}{'label'} = $label;
		} elsif ($line =~ /fastq|fq/) {			# if only one fastq file given in command assume it is unpaired
			$label = $1 if ($line =~ /\S+\s+(\S+)/); 
			$unpaired_fq = $line;
			die "Cannot locate $unpaired_fq fastq file\n\n" if (! -e $unpaired_fq); 	# if fq file is not found
			$id = $unpaired_fq;
			$id =~ s/.*\///; 
			$id =~ s/\.fastq|\.fq//;
			$id =~ s/\.r\d$//i;                     # strip off any "r2" or "R2" suffix 
			$input{$id}{'fq1'}     = $unpaired_fq;
			$input{$id}{'fq2'}     = "na";
			$input{$id}{'paired'}  = "false";
			$input{$id}{'aligned'} = "false";
			$input{$id}{'id'}      = $id;
		} else {
			print "\nCheck filenames in $input_list.  BAMs must be named with \"bam\" and fastq files named with either \"fastq\" or \"fq\".\n\n";
		}
	}
	close IN;
}

sub evaluatePairing {
	print "\nEvaluating if BAM files are paired/unpaired and aligned/unaligned:\n";
	foreach my $id (keys %input) {
		my $bam = $input{$id}{'bam'};
		my ($total, $mapped, $pct_aligned) = (0, 0, "na");
		my ($paired, $r1_reads, $r2_reads) = (0, 0, 0);
		next if (! $bam);
		#print "$id\t$bam\n";
		open (FLAG,"$samtools flagstat $bam | ");
		while(<FLAG>){
			my $line = $_;
			chomp($line);
			if ($line =~ /(\d+)\s+\+\s+(\d+)\s+in\s+total\s+/) {
				$total = $1 + $2;
			} elsif ($line =~ /(\d+)\s+\+\s+(\d+)\s+paired\s+in\s+sequencing/) {
				$paired = 1 if ($1 + $2 > 0);;
			} elsif ($line =~ /(\d+)\s+\+\s+(\d+)\s+read1/) {
				$r1_reads = 1 if ($1 + $2 > 0);
		 	} elsif ($line =~ /(\d+)\s+\+\s+(\d+)\s+read2/) {
				$r2_reads = 1 if ($1 + $2 > 0);
			} elsif ($line =~ /(\d+)\s+\+\s+(\d+)\s+mapped/) {
				if ($1 > 0 || $2 > 0) {
					$input{$id}{'aligned'} = "true";
					$mapped = $1 + $2;
					$input{$id}{'pct_aligned'} = sprintf("%.1f",(($mapped/$total)*100));
				} else {
					$input{$id}{'aligned'} = "false";
					$input{$id}{'pct_aligned'} = "na";
				}	
			}
		}
		close FLAG;
		if ($paired == 1 && $r1_reads > 0 && $r2_reads > 0) {
			$input{$id}{'paired'} = "true";
			$plot_insert_size = 1;
		} else {
			$input{$id}{'paired'} = "false";
		}	
		print "$bam:\tPaired = $input{$id}{'paired'}\tAligned = $input{$id}{'aligned'}\n\n";
	}
}

sub makeOutputDirs {
	`mkdir $qc_main_subdir`        if (! -e $qc_main_subdir);
	`mkdir $qc_read_data_subdir`   if (! -e $qc_read_data_subdir);
	`mkdir $qc_metric_data_subdir` if (! -e $qc_metric_data_subdir);
	`mkdir $qc_chart_data_subdir`  if (! -e $qc_chart_data_subdir);
}

sub makeSeparateBams {
	print "Creating separate read1 and read2 BAM files\n\n";
	my ($out_bam_r1, $out_bam_r2, $out_bam_pair, $fq2sam_r1_cmd, $fq2sam_r2_cmd, $fq2sam_pair_cmd, $sam2fq_cmd, $sam2fq_r1_cmd, $sam2fq_r2_cmd);
	my ($make_r1_bam_cmd, $make_r2_bam_cmd);
	foreach my $id (keys %input) {
		####################  unpaired / fq1  ########
		if ($input{$id}{'paired'} eq "false" && defined $input{$id}{'fq1'} && not defined $input{$id}{'bam'}) {
			# create unpaired BAM from fq1 fastq file
			$out_bam_r1 = "$qc_read_data_subdir/$id.bam";
			$fq2sam_r1_cmd = "java -jar $picard/picard.jar FastqToSam  VALIDATION_STRINGENCY=SILENT F1=$input{$id}{'fq1'} O=$out_bam_r1 SAMPLE_NAME=sample";
			#print "ID $id -   Unpaired R1 fastq file is input:\n$fq2sam_r1_cmd\n";
			system "$fq2sam_r1_cmd" unless ($no_overwrite && -e "$out_bam_r1");
			$input{$id}{'r1bam'} = "$qc_read_data_subdir/$id.bam";
			$input{$id}{'r2bam'} = "na";
		} 
		####################  unpaired / bam ########
		elsif ($input{$id}{'paired'} eq "false" && not defined $input{$id}{'fq1'} && defined $input{$id}{'bam'}) {
			# nothing to do but keep track of bam as working bam
			#print "ID $id -   Unpaired BAM file is input.  No BAM creation needed.\n\n";
			$input{$id}{'r1bam'} = $input{$id}{'bam'};
			$input{$id}{'r2bam'} = "na";
			if ($input{$id}{'aligned'} eq "true") {		# if aligned make sure sorted by coord and get ref IDs
				$input{$id}{'sorted_aligned_bam'}  = &sortAndFindRef($id, "unpaired");
			}
		}
		####################  paired / fq1 and fq2  ########
		elsif ($input{$id}{'paired'} eq "true" && defined $input{$id}{'fq1'} && defined $input{$id}{'fq2'}) {
			# create r1 BAM and r2 BAM from r1 fq and r2 fq

			#print "Paired R1 fastq and R2 fastq is input.  Created R1 BAM and R2 BAM:\n";
			$out_bam_r1   = "$qc_read_data_subdir/$id.r1.bam";
			$out_bam_r2   = "$qc_read_data_subdir/$id.r2.bam";
			$out_bam_pair = "$qc_read_data_subdir/$id.paired.bam";
			$fq2sam_r1_cmd   = "java -jar $picard/picard.jar FastqToSam  VALIDATION_STRINGENCY=SILENT F1=$input{$id}{'fq1'} O=$qc_read_data_subdir/$id.r1.bam SAMPLE_NAME=sample";
			$fq2sam_r2_cmd   = "java -jar $picard/picard.jar FastqToSam  VALIDATION_STRINGENCY=SILENT F1=$input{$id}{'fq2'} O=$qc_read_data_subdir/$id.r2.bam SAMPLE_NAME=sample";
			$fq2sam_pair_cmd = "java -jar $picard/picard.jar FastqToSam  VALIDATION_STRINGENCY=SILENT F1=$input{$id}{'fq1'} F2=$input{$id}{'fq2'} O=$qc_read_data_subdir/$id.paired.bam SAMPLE_NAME=sample";
			system "$fq2sam_r1_cmd"   unless ($no_overwrite && -e "$out_bam_r1");
			system "$fq2sam_r2_cmd"   unless ($no_overwrite && -e "$out_bam_r2");
			system "$fq2sam_pair_cmd" unless ($no_overwrite && -e "$out_bam_pair");
			$input{$id}{'r1bam'} = "$qc_read_data_subdir/$id.r1.bam";
			$input{$id}{'r2bam'} = "$qc_read_data_subdir/$id.r2.bam";
			$input{$id}{'bam'}   = "$qc_read_data_subdir/$id.paired.bam";
		}
		####################  paired / bam ########
		elsif ($input{$id}{'paired'} eq "true" && ! exists $input{$id}{'fq1'} && ! exists $input{$id}{'fq2'} &&  exists $input{$id}{'bam'}) {
			# create r1 BAM and r2 BAM from paired BAM
			#print "ID $id - Paired BAM file is unput.  Create r1 BAM and r2 BAM with samtools view -f <flag>\n";
			$make_r1_bam_cmd = "$samtools view -h -f 64  $input{$id}{'bam'} -b -o $qc_read_data_subdir/$id.r1.bam"; 
			$make_r2_bam_cmd = "$samtools view -h -f 128 $input{$id}{'bam'} -b -o $qc_read_data_subdir/$id.r2.bam"; 
			system "$make_r1_bam_cmd\n" unless ($no_overwrite && -e "$qc_read_data_subdir/$id.r1.bam");
			system "$make_r2_bam_cmd\n" unless ($no_overwrite && -e "$qc_read_data_subdir/$id.r2.bam");
			$input{$id}{'r1bam'} = "$qc_read_data_subdir/$id.r1.bam";
			$input{$id}{'r2bam'} = "$qc_read_data_subdir/$id.r2.bam";
			if ($input{$id}{'aligned'} eq "true") { 				# sort by coordinate and find ref
				$input{$id}{'sorted_aligned_bam'}  = &sortAndFindRef($id, "paired");
			}
		}

	}
}


sub sortAndFindRef {
	my ($id, $paired) = ($_[0], $_[1]);
	my ($ref, $out_bam);
	open (VIEW,"$samtools view -h $input{$id}{'bam'} | egrep \"\@\" | ");
	while(<VIEW>) {
		my $line = $_;
		chomp($line);
		if ($line =~ /^\@HD.*SO\:(\S+)$/) {
			if ($1 ne "coordinate") {
				$out_bam = "$qc_read_data_subdir/$id.$paired.sorted.bam";
				my $sort_sam_cmd = "java -jar $picard/picard.jar SortSam VALIDATION_STRINGENCY=SILENT I=$input{$id}{'bam'} O=$qc_read_data_subdir/$id.$paired.sorted.bam VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate";
				print "Running SortSam to sort the aligned $paired BAM into coordinate sort order...\n\n";
				system "$sort_sam_cmd" unless ( $no_overwrite && -e "$qc_read_data_subdir/$id.$paired.sorted.bam");
				$input{$id}{'sorted_aligned_bam'} = "$qc_read_data_subdir/$id.$paired.sorted.bam";
			} else {
				$out_bam = $input{$id}{'bam'};
			}
		} elsif ($line =~ /^\@SQ.*UR\:(\S+)/)  {
			$ref = $1;
			$ref =~ s/^file\://;			# some BAM files have format of UR:file:/path/to/reference
			$input{$id}{'reffile'} = $ref;
			$ref =~ s/.*\///;
			$input{$id}{'ref'} = $ref;
			$references{$ref}{'reffile'} = $input{$id}{'reffile'}; 
			$references{$ref}{'id'} = 1 if (! exists $references{$ref});
			last;		# found reference now get out of here since Tophat might be below and turn to relative path
		} elsif ($line =~ /ID\:TopHat\s+VN\S+\s+CL\:.*\s+(\S+)\s+\S+fastq\s+\S+fastq/) {
			$ref = $1 . ".fa";
			$input{$id}{'reffile'} = $ref;
			$ref =~ s/.*\///;
			$input{$id}{'ref'} = $ref;
			$references{$ref}{'id'} = 1 if (! exists $references{$ref});
			last;
		} else {
			#print "Warning:  Reference not found $line\n";
		}
	}
	if (! defined $ref) {
		if (defined $aligned_ref) {
			print "Using user defined -A <ref>\n";
			$ref = $aligned_ref;
			$ref =~ s/.*\///;
			$input{$id}{'ref'} = $ref;
			$input{$id}{'reffile'} = $aligned_ref;
			$references{$ref}{'id'} = $aligned_ref;
		}
	}
	close VIEW;
	return($out_bam);
}

sub getPicardColumn {
	my ($line,$match) = @_;		 # get our line array ref and what to match on

    	my @line = @{$line};		 # dereference 

    	for (my $i = 0; $i <= $#line; $i++) {
       	 	if ($line[$i] eq $match) {
            	return $i;		 # go through array and return the index of the match
        }
        # go through array and return the index of the match
    	}
    	return -1;		 #otherwise we didn't find a match

}
sub runFastQcAndPicardModules {
	foreach my $id (keys %input) {
		my ($bam1, $bam2, $paired_bam, $sorted_bam) = ($input{$id}{'r1bam'}, $input{$id}{'r2bam'}, $input{$id}{'bam'}, $input{$id}{'sorted_aligned_bam'});
		my $out_bam1       = $bam1;
		my $out_bam2       = $bam2;
		my $out_paired_bam = $paired_bam;
		$out_bam1          =~ s/.*\///;
		$out_bam2          =~ s/.*\///;
		$out_paired_bam    =~ s/.*\/// if ($paired_bam);
		#########################################                 FastQC                   #################################
		print "Running FastQC...\n";
		`mkdir $qc_fastqc_subdir` if (! -e $qc_fastqc_subdir);
		my ($r1_fastqc_outdir, $r2_fastqc_outdir) = ($out_bam1, $out_bam2);
		$r1_fastqc_outdir =~ s/\.bam$/\_fastqc/;
		$r2_fastqc_outdir =~ s/\.bam$/\_fastqc/;
		my $r1_fastqc_cmd     = "$fastqc  --extract --out $qc_fastqc_subdir $bam1";
		if ($no_overwrite && -e "$qc_fastqc_subdir/$r1_fastqc_outdir/Images/per_base_quality.png") {
			print "The -N no_overwrite option set and read1 FastQC output exists.  Skipping FastQC...\n\n";
		} else {
			print "\t$r1_fastqc_cmd....\n\n";
			system "$r1_fastqc_cmd" unless ($no_overwrite && -e "$r1_fastqc_outdir/Images/per_base_quality.png");
		}
		if ($input{$id}{'paired'} eq "true")  {
			if ($no_overwrite && -e "$qc_fastqc_subdir/$r2_fastqc_outdir/Images/per_base_quality.png") {
				print "The -N no_overwrite option set and read2 FastQC output exists.  Skipping FastQC...\n\n";
			} else {
				my $r2_fastqc_cmd     = "$fastqc  --extract --out $qc_fastqc_subdir $bam2";
				print "\t$r2_fastqc_cmd\n\n";
				system "$r2_fastqc_cmd" unless ($no_overwrite && -e "$r2_fastqc_outdir/Images/per_base_quality.png");
			}
		}
		$input{$id}{'r1_fastqc_dir'} = "$qc_fastqc_subdir/$r1_fastqc_outdir";
		$input{$id}{'r2_fastqc_dir'} = "$qc_fastqc_subdir/$r2_fastqc_outdir";
		#########################################      Collect Quality Yield Metrics       #################################
		print "Running CollectQualityYieldMetrics...\n\n";
		my ($r1_q_yield_metrics)  = ("$qc_metric_data_subdir/$out_bam1.quality_yield.metrics");
		my ($r2_q_yield_metrics)  = ("$qc_metric_data_subdir/$out_bam2.quality_yield.metrics");
		my $r1_qual_yield_cmd = "java -jar $picard/picard.jar CollectQualityYieldMetrics VALIDATION_STRINGENCY=SILENT I=$bam1 O=$r1_q_yield_metrics";
		print "$r1_qual_yield_cmd\n\n";
		system "$r1_qual_yield_cmd"  unless ($no_overwrite && -e "$r1_q_yield_metrics");
		if ($input{$id}{'paired'} eq "true")  {
			my $r2_qual_yield_cmd = "java -jar $picard/picard.jar CollectQualityYieldMetrics VALIDATION_STRINGENCY=SILENT I=$bam2 O=$r2_q_yield_metrics";
			print "$r2_qual_yield_cmd\n\n";
			system "$r2_qual_yield_cmd"   unless ($no_overwrite && -e "$r2_q_yield_metrics");

		}
		$input{$id}{'r1_qual_yield_metrics'} = $r1_q_yield_metrics;
		$input{$id}{'r2_qual_yield_metrics'} = $r2_q_yield_metrics;
		#########################################  Mean Base Quality Along Read Position   #################################
		print "Running MeanQualityByCycle...\n\n";
		my ($r1_metrics, $r1_chart)  = ("$qc_metric_data_subdir/$out_bam1.quality_by_cycle.metrics", "$qc_metric_data_subdir/$out_bam1.quality_by_cycle.pdf");
		my ($r2_metrics, $r2_chart)  = ("$qc_metric_data_subdir/$out_bam2.quality_by_cycle.metrics", "$qc_metric_data_subdir/$out_bam2.quality_by_cycle.pdf");
		my $r1_qual_cmd = "java -jar $picard/picard.jar MeanQualityByCycle VALIDATION_STRINGENCY=SILENT I=$bam1 O=$r1_metrics CHART_OUTPUT=$r1_chart";
		system "$r1_qual_cmd" unless ($no_overwrite && -e "$r1_metrics");
		if ($input{$id}{'paired'} eq "true")  {
			my $r2_qual_cmd = "java -jar $picard/picard.jar MeanQualityByCycle VALIDATION_STRINGENCY=SILENT I=$bam2 O=$r2_metrics CHART_OUTPUT=$r2_chart";
			system "$r2_qual_cmd"  unless ($no_overwrite && -e "$r2_metrics");

		}
		$input{$id}{'r1_qual_metrics'} = $r1_metrics;
		$input{$id}{'r2_qual_metrics'} = $r2_metrics;
		#################################              Quality Score Distribution          #################################
		print "Running QualityScoreDistribution...\n\n";
		my ($r1_qual_hist_metrics, $r1_qual_hist_chart)  = ("$qc_metric_data_subdir/$out_bam1.quality_histo.metrics", "$qc_metric_data_subdir/$out_bam1.quality_histo.pdf");
		my ($r2_qual_hist_metrics, $r2_qual_hist_chart)  = ("$qc_metric_data_subdir/$out_bam2.quality_histo.metrics", "$qc_metric_data_subdir/$out_bam2.quality_histo.pdf");
		my $r1_qual_hist_cmd = "java -jar $picard/picard.jar QualityScoreDistribution VALIDATION_STRINGENCY=SILENT I=$bam1 O=$r1_qual_hist_metrics CHART_OUTPUT=$r1_qual_hist_chart";
		system "$r1_qual_hist_cmd"  unless ($no_overwrite && -e "$r1_qual_hist_metrics");
		if ($input{$id}{'paired'} eq "true")  {
			my $r2_qual_hist_cmd = "java -jar $picard/picard.jar QualityScoreDistribution VALIDATION_STRINGENCY=SILENT I=$bam2 O=$r2_qual_hist_metrics CHART_OUTPUT=$r2_qual_hist_chart";
			system "$r2_qual_hist_cmd"  unless ($no_overwrite && -e "$r2_qual_hist_metrics");

		}
		$input{$id}{'r1_qual_hist_metrics'} = $r1_qual_hist_metrics;
		$input{$id}{'r2_qual_hist_metrics'} = $r2_qual_hist_metrics;
		##########################################                 Read GC                 #################################
		print "Running CollectReadGCMetrics...\n\n";
		my ($r1_gc_metrics, $r1_gc_chart)  = ("$qc_metric_data_subdir/$out_bam1.read_gc.metrics", "$qc_metric_data_subdir/$out_bam1.read_gc.pdf");
		my ($r2_gc_metrics, $r2_gc_chart)  = ("$qc_metric_data_subdir/$out_bam2.read_gc.metrics", "$qc_metric_data_subdir/$out_bam2.read_gc.pdf");
		my $r1_gc_cmd = "java -jar ${aa_picard}CollectReadGCMetrics.jar VALIDATION_STRINGENCY=SILENT I=$bam1 O=$r1_gc_metrics CHART_OUTPUT=$r1_gc_chart";
		#print "$id\t$r1_gc_cmd\n";
		system "$r1_gc_cmd" unless ($no_overwrite && -e "$r1_gc_metrics");
		if ($input{$id}{'paired'} eq "true")  {
			my $r2_gc_cmd = "java -jar ${aa_picard}CollectReadGCMetrics.jar VALIDATION_STRINGENCY=SILENT I=$bam2 O=$r2_gc_metrics CHART_OUTPUT=$r2_gc_chart";
			system "$r2_gc_cmd" unless ($no_overwrite && -e "$r2_gc_metrics");

		}
		$input{$id}{'r1_gc_metrics'} = $r1_gc_metrics;
		$input{$id}{'r2_gc_metrics'} = $r2_gc_metrics;
		#########################################     AlignmentSummaryMetrics         ######################################
		##########################  Single command for paired/unpaired.  does not work on separate r1 and r2 ###############
		print "Running AlignmentSummaryMetrics...\n\n";
		my ($asm_cmd, $asm_metrics) = ("", "");
		if ($input{$id}{'paired'} eq "true") {
			$asm_metrics = "$qc_metric_data_subdir/$out_paired_bam.alignment_summary.metrics";
			$asm_cmd = "java -jar $picard/picard.jar CollectAlignmentSummaryMetrics VALIDATION_STRINGENCY=SILENT I=$paired_bam O=$asm_metrics";
		} else {
			$asm_metrics = "$qc_metric_data_subdir/$out_bam1.alignment_summary.metrics";
			$asm_cmd = "java -jar $picard/picard.jar CollectAlignmentSummaryMetrics VALIDATION_STRINGENCY=SILENT I=$bam1 O=$asm_metrics";
		}	
		$asm_cmd .= " $adapter_string" if ($adapter_file);			# tack on the ADAPTER_SEQUENCE parameter
		system "$asm_cmd"  unless ($no_overwrite && -e "$asm_metrics");
		$input{$id}{'asm_metrics'} = $asm_metrics;
		#########################################      Library Complexity (two methods)    #################################
		#if (! $no_complexity && $input{$id}{'aligned'} eq "true") {
		#	print "Running ComplexityByStartPos...\n";
		#	my ($complexity_by_start_metrics, $complexity_by_start_cmd, $complexity_by_start_counts, $complexity_by_start_txt);
		#	if ($input{$id}{'paired'} eq "true") {
		#		$complexity_by_start_metrics = ("$qc_metric_data_subdir/$out_paired_bam.complexity_by_start.metrics");
		#		$complexity_by_start_cmd = "$complexity_by_start -bam $paired_bam -out $complexity_by_start_metrics ";
		#	} else {
		#		$complexity_by_start_metrics =  ("$qc_metric_data_subdir/$out_bam1.complexity_by_start.metrics");
		#		$complexity_by_start_cmd = "$complexity_by_start -single_end -bam $bam1 -out $complexity_by_start_metrics ";
		#	}
		#	$complexity_by_start_counts = $complexity_by_start_metrics . "_counts.txt";
		#	$complexity_by_start_txt    = $complexity_by_start_metrics . "_rarefactionData.txt";
		#	if ($no_overwrite && -e "$complexity_by_start_txt") {
		#		print "ComplexityByStartPos output exists and -N used - skipping ComplexityByStartPos...\n\n"
		#	} else {
		#		print "Command:\n$complexity_by_start_cmd\n\n";
		#		system "$complexity_by_start_cmd" unless  ($no_overwrite && -e "$complexity_by_start_metrics");
		#	}
		#	$input{$id}{'complexity_by_start_metrics'} = $complexity_by_start_txt ;
		#	unlink $complexity_by_start_counts;
		#}
		if (! $no_complexity && $input{$id}{'paired'} eq "true") {
			my ($bam_comp_metrics) = ("$qc_metric_data_subdir/$out_paired_bam.complexity.metrics");  
			print "Running EstimateLibraryComplexity...\n\n";
			my $complexity_cmd = "java -jar $picard/picard.jar EstimateLibraryComplexity VALIDATION_STRINGENCY=SILENT I=$paired_bam O=$bam_comp_metrics";
			system "$complexity_cmd" unless ($no_overwrite && -e "$bam_comp_metrics");
			$input{$id}{'r1_comp_metrics'}     = "na";
			$input{$id}{'r2_comp_metrics'}     = "na";
			$input{$id}{'paired_comp_metrics'} = "$bam_comp_metrics";

		}
		##################################         Jump  Metrics (paired+aligned only)      ###############################
		if ($input{$id}{'paired'} eq "true" && $input{$id}{'aligned'} eq "true") {
			my ($jump_metrics) = ("$qc_metric_data_subdir/$out_paired_bam.jump.metrics");  
			print "Running CollectJumplingLibraryMetrics...\n\n";
			my $jump_metrics_cmd = "java -jar $picard/picard.jar CollectJumpingLibraryMetrics I=$sorted_bam O=$jump_metrics VALIDATION_STRINGENCY=SILENT";
			system "$jump_metrics_cmd" unless ($no_overwrite && -e "$jump_metrics");
			$input{$id}{'r1_jump_metrics'}     = "na";
			$input{$id}{'r2_jump_metrics'}     = "na";
			$input{$id}{'jump_metrics'} = "$jump_metrics";

		}
		##################################     Insert Size Metrics (paired+aligned only)    ###############################
		if ($input{$id}{'paired'} eq "true" && $input{$id}{'aligned'} eq "true") {
			my ($bam_insert_size_metrics, $bam_insert_size_chart) = ("$qc_metric_data_subdir/$out_paired_bam.insert_size.metrics", "$qc_metric_data_subdir/$out_paired_bam.insert_size.pdf");  
			print "Running CollectInsertSizeMetrics...\n\n";
			my $insert_size_metrics_cmd = "java -jar $picard/picard.jar CollectInsertSizeMetrics W=15000 I=$sorted_bam O=$bam_insert_size_metrics H=$bam_insert_size_chart VALIDATION_STRINGENCY=SILENT";
			system "$insert_size_metrics_cmd" unless ($no_overwrite && -e "$bam_insert_size_metrics");
			$input{$id}{'r1_comp_metrics'}     = "na";
			$input{$id}{'r2_comp_metrics'}     = "na";
			$input{$id}{'paired_insert_size_metrics'} = "$bam_insert_size_metrics";

		}
		##################################         CoverageAlongReference (aligned only)     ###############################
		if (! $no_coverage && $input{$id}{'aligned'} eq "true") {
			my ($cvg_metrics, $cvg_chart, $cvg_cmd, $cwd_pdf, $cwd_histo_pdf);
			print "Running CoverageAlongReference...\n\n";
			if ($input{$id}{'paired'} eq "true") {
				($cvg_metrics, $cvg_chart) = ("$qc_metric_data_subdir/$out_paired_bam.coverage_along_reference.metrics", "$qc_metric_data_subdir/$out_paired_bam.coverage_along_reference.pdf");  
				$cvg_cmd = "java -jar $aa_picard/CoverageAlongReference.jar I=$paired_bam SC=$cvg_chart SEQ_OUTPUT=$cvg_metrics VALIDATION_STRINGENCY=SILENT";
			} else {
				($cvg_metrics, $cvg_chart) = ("$qc_metric_data_subdir/$out_bam1.coverage_along_reference.metrics", "$qc_metric_data_subdir/$out_bam1.coverage_along_reference.pdf");  
				$cvg_cmd = "java -jar $aa_picard/CoverageAlongReference.jar I=$bam1 SEQ_OUTPUT=$cvg_metrics SC=$cvg_chart VALIDATION_STRINGENCY=SILENT";
			}
			system "$cvg_cmd" unless ($no_overwrite && -e "$cvg_metrics");
			$input{$id}{'seq_cvg_metrics'} = "$cvg_metrics";
			$cwd_pdf = $cvg_chart;						# CoverageAlongReference has quirk ..
			$cwd_pdf =~ s/.*\///;						# in writing output to cwd even if cmd ..
			$cwd_histo_pdf = $cwd_pdf;					# directs it to some other dir 
			$cwd_histo_pdf =~ s/\.pdf$/\.histogram\.pdf/;			# So these files moved to keep clean
			system "mv $cwd_pdf $qc_metric_data_subdir/" if (-e $cwd_pdf);		
			system "mv $cwd_histo_pdf $qc_metric_data_subdir/" if (-e $cwd_histo_pdf);	

		}
		##################################         GC  Bias (aligned only)     ###############################
		if (! $no_coverage && $input{$id}{'aligned'} eq "true" && -e $input{$id}{'reffile'}) {
			my $ref = $input{$id}{'reffile'};
			my ($gc_bias_cmd, $gc_bias_metrics, $gc_bias_summary, $gc_bias_chart);
			print "Running CollectGcBiasMetrics...\n\n";
			if ($input{$id}{'paired'} eq "true") {
				($gc_bias_metrics, $gc_bias_chart, $gc_bias_summary) = ("$qc_metric_data_subdir/$out_paired_bam.gc_bias.metrics", "$qc_metric_data_subdir/$out_paired_bam.gc_bias.pdf", "$qc_metric_data_subdir/$out_paired_bam.gc_bias.summary_metrics");  
				$gc_bias_cmd = "java -jar $picard/picard.jar CollectGcBiasMetrics SUMMARY_OUTPUT=$gc_bias_summary REFERENCE_SEQUENCE=$ref I=$paired_bam CHART_OUTPUT=$gc_bias_chart O=$gc_bias_metrics VALIDATION_STRINGENCY=SILENT";
			} else {
				($gc_bias_metrics, $gc_bias_chart, $gc_bias_summary) = ("$qc_metric_data_subdir/$out_bam1.gc_bias.metrics", "$qc_metric_data_subdir/$out_bam1.gc_bias.pdf", "$qc_metric_data_subdir/$out_bam1.gc_bias.summary_metrics");  
				$gc_bias_cmd = "java -jar $picard/picard.jar CollectGcBiasMetrics SUMMARY_OUTPUT=$gc_bias_summary REFERENCE_SEQUENCE=$ref I=$bam1 CHART_OUTPUT=$gc_bias_chart O=$gc_bias_metrics VALIDATION_STRINGENCY=SILENT";
			}
			system "$gc_bias_cmd" unless ($no_overwrite && -e "$gc_bias_metrics");
			$input{$id}{'gc_bias_metrics'} = "$gc_bias_metrics";
		} elsif ($input{$id}{'aligned'} eq "true" && ! -e $input{$id}{'reffile'}) {
			print "Aligned BAM but cannot locate reference files .... unable to run CollectGcBiasMetrics...\n\n";
		} else {
			#print "$no_coverage $input{$id}{'aligned'} $input{$id}{'reffile'}\n";
		}
	}
}


sub getPercentPfFromFastq {
	my ($fq1, $pct_poly_a, $pct_poly_t, $pct_pf, $pf_count, $filt,$pf_metrics_out);
	my $total_count = 0;
	my $seq_line = 0;				# know that next line is sequence line
	my $poly_a_count = 0;
	my $poly_t_count = 0;
	my $a_string = "AAAAAAAAAAAAAAAAAAAA";
	my $t_string = "TTTTTTTTTTTTTTTTTTTT";
	my $name;
	foreach my $id (keys %input) {
		$fq1 = $input{$id}{'fq1'} ;		
		$pf_metrics_out = "$qc_metric_data_subdir/$id.pf_metrics.txt"; 
		if (! -e $pf_metrics_out) {
			print "Parsing fastq file to get PF counts...\n\n";
			open(PF,">$pf_metrics_out") || die "Cannot write to the $pf_metrics_out out file:  $!\n\n";
			print PF "total\tpf_reads\tpct_pf\tpoly_a_reads\tpct_poly_a\tpoly_t_reads\tpct_poly_t\n";
			next if (! -e $fq1);
			if ($fq1 =~ /\.gz$/) {
				open (IN,"zcat $fq1 | ") || warn "cannot zcat the $fq1 file\n\n";
			} else {
				open(IN,"$fq1") || warn "$fq1 not found\n\n";
			}
			while (<IN>) {
				my $line = $_;
				chomp($line);
				if ($line =~ /^\@.*\:(Y|N)\:/) {
					#print "FOUND\t$line\n";
					$filt = $1;
					$total_count++;	
					$pf_count++ if ($filt eq "N");
					$seq_line = 1;
					$name = $line;
				} elsif ($seq_line == 1)  {
					if  ($line =~ /$a_string/) {
						#print "POLY-A\t$name\t$line\n";
						$poly_a_count++;
						$seq_line = 0;
					} elsif ($line =~ /$t_string/) {
						print "POLY-T\t$name\t$line\n";
						$poly_t_count++;
						$seq_line = 0;
					} else {
						#print "NO_POLY\t$name\t$line\n";
						# poly-a and poly-t not found in sequence
					}	
				} else {		
					#print "TEST\n$line\n";
				}	
			}
			close IN;
			$pct_pf                     = 100*($pf_count/$total_count) if ($total_count > 0);
			$pct_poly_a                 = sprintf("%.1f",(100*($poly_a_count/$total_count))) if ($total_count > 0);
			$pct_poly_t                 = sprintf("%.1f",(100*($poly_t_count/$total_count))) if ($total_count > 0);
			$input{$id}{'pct_pf'}       = $pct_pf;
			$input{$id}{'pct_poly_a'}   = $pct_poly_a;
			$input{$id}{'pct_poly_t'}   = $pct_poly_t;
			print PF "$total_count\t$pf_count\t$pct_pf\t$poly_a_count\t$pct_poly_a\t$poly_t_count\t$pct_poly_t\n";
			close PF;
		} else {
			print "Found existing $pf_metrics_out file...\n\n";
			open(PF, "$pf_metrics_out") || warn "Cannot open the $pf_metrics_out file\n\n";
			while (<PF>) {
				my $line = $_;
				chomp($line);
				if ($line =~ /^\d+\s+\d+\s+(\S+)\s+\d+\s+(\S+)\s+\d+\s+(\S+)/) {
					$input{$id}{'pct_pf'}     = $1;
					$input{$id}{'pct_poly_a'} = $2;
					$input{$id}{'pct_poly_t'} = $3;
				}
			}
			close PF;
		}
	}
}

sub getPercentAlignedByRead {
	foreach my $id (keys %input) {
		#print "getPercentAlignedByRead $id\n";
		my ($bam1, $bam2) = ($input{$id}{'r1bam'}, $input{$id}{'r2bam'});
		#print "getPercentAlignedByRead $bam1 $bam2\n$input{$id}{'paired'} $input{$id}{'aligned'}\n";
		if ($input{$id}{'paired'} ne "true") {
			$input{$id}{'r2_pct_aligned'} = "na" if ($input{$id}{'paired'} ne "true");
			#next;					# continue in case unpaired but aligned
		}
		if ($input{$id}{'aligned'} ne "true") {
			#$input{$id}{'r1_pct_aligned'} = "na" if ($input{$id}{'paired'} ne "true");
			#$input{$id}{'r2_pct_aligned'} = "na" if ($input{$id}{'paired'} ne "true");
			$input{$id}{'r1_pct_aligned'} = "na";
			$input{$id}{'r2_pct_aligned'} = "na";
			next;
		} 
		foreach my $bam ($bam1, $bam2) {
			next if (! -e $bam || $bam eq "na");	# in case we have aligned BAM but unpaired
			my ($read, $total) = ( -1, -1);
			my $mapped = -1;

			open (FLAG,"$samtools flagstat $bam | ");
			while(<FLAG>){
				my $line = $_;
				chomp($line);
				if ($line =~ /(\d+)\s+\+\s+(\d+)\s+in\s+total\s+/) {
					$total = $1 + $2;
				} elsif ($line =~ /(\d+)\s+\+\s+(\d+)\s+mapped/) {
					$mapped = $1 + $2;
				} elsif ($line =~ /(\d+)\s+\+\s+(\d+)\s+read(\d+)/) {	# find if read1 or read2 where mapped > 0
					$read = 1 if ($1 + $2 > 0 && $3 == 1);		# example flagstat line: 500 + 10 read1
					$read = 2 if ($1 + $2 > 0 && $3 == 2);
					$input{$id}{'r1_pct_aligned'} = sprintf("%.1f",(($mapped/$total)*100)) if  ($read == 1 && $mapped > 0);
					$input{$id}{'r2_pct_aligned'} = sprintf("%.1f",(($mapped/$total)*100)) if  ($read == 2 && $mapped > 0);
				}
			}
			close FLAG;
		}
	}
}

sub getQualityYieldMetrics {
	foreach my $id (keys %input) {
		foreach my $metrics ($input{$id}{'r1_qual_yield_metrics'}, $input{$id}{'r2_qual_yield_metrics'}) {
			my ($read, $total_reads, $pf_reads, $pct_pf_reads, $read_ln, $total_bases);
			my ($q20_bases, $q30_bases);
			$input{$id}{'r2_pct_q20'}     = "na" if ($input{$id}{'paired'} ne "true");
			$input{$id}{'r2_pct_q30'}     = "na" if ($input{$id}{'paired'} ne "true");
			$input{$id}{'r2_rd_ln'}       = "na" if ($input{$id}{'paired'} ne "true");
			$input{$id}{'r2_total_reads'} = "na" if ($input{$id}{'paired'} ne "true");
			next if (($metrics =~ /\.r2\.bam/ || $metrics =~ /na.quality_yield.metrics/) && $input{$id}{'paired'} ne "true" );
			if ($metrics =~ /\.r1\.bam/) {
				$read = 1;
			} elsif ($metrics =~ /\.r2\.bam/) {
				$read = 2;
			} else {					# unpaired, so make this read1
				#print "CHECK $metrics\n";
				$read = 1;
			}
			my $start_count = 0;
			open(IN,$metrics) || die "Cannot open $metrics metrics file:$!\n\n";
			while(<IN>) {
				my $line = $_;
				chomp($line);
				my @line = split(/\s+/,$line);
				if ($line =~ /^TOTAL_READS/) {
					$start_count = 1;
					next;
				}
				if ($start_count != 1) {
					next;
				} else {
					($total_reads, $pf_reads, $read_ln, $total_bases, $q20_bases, $q30_bases) = ($line[0], $line[1], $line[2], $line[3], $line[5], $line[7]);
					$start_count = 0;	
				}	
			}
			close IN;
			$input{$id}{'r1_pct_q20'}     = sprintf("%.1f",(($q20_bases/$total_bases)*100)) if  ($read == 1);
			$input{$id}{'r2_pct_q20'}     = sprintf("%.1f",(($q20_bases/$total_bases)*100)) if  ($read == 2);
			$input{$id}{'r1_pct_q30'}     = sprintf("%.1f",(($q30_bases/$total_bases)*100)) if  ($read == 1);
			$input{$id}{'r2_pct_q30'}     = sprintf("%.1f",(($q30_bases/$total_bases)*100)) if  ($read == 2);
			$input{$id}{'r1_total_reads'} = $total_reads  if  ($read == 1);
			$input{$id}{'r2_total_reads'} = $total_reads  if  ($read == 2);
			$input{$id}{'r1_pf_reads'}    = $pf_reads     if  ($read == 1);
			$input{$id}{'r2_pf_reads'}    = $pf_reads     if  ($read == 2);		# should be same no. as r1 PF
			$input{$id}{'r1_mean_rd_ln'}  = $read_ln      if  ($read == 1);
			$input{$id}{'r2_mean_rd_ln'}  = $read_ln      if  ($read == 2);
			if (! defined $input{$id}{'pct_pf'} ) {				# % PF may already have been calculated
				$input{$id}{'pct_pf'} = sprintf("%.1f",(($input{$id}{'r1_pf_reads'}/$input{$id}{'r1_total_reads'})*100)) if  ($input{$id}{'r1_total_reads'} > 0) # PF refers to pairs so no need to break down r1 and r2 for pct pf;
			}

		}
	}
}

sub getQualityByPosition {
	foreach my $id (keys %input) {
		foreach my $metrics ($input{$id}{'r1_qual_metrics'}, $input{$id}{'r2_qual_metrics'}) {
			my ($read, $sum, $avg);
			my $count;
			$input{$id}{'r2_mean_qual'} = "na" if ($input{$id}{'paired'} ne "true");
			next if (($metrics =~ /\.r2\.bam/ || $metrics =~ /na.quality_by_cycle.metrics/) && $input{$id}{'paired'} ne "true" );
			if ($metrics =~ /\.r1\.bam/) {
				$read = 1;
			} elsif ($metrics =~ /\.r2\.bam/) {
				$read = 2;
			} else {					# unpaired, so make this read1
				#print "CHECK $metrics\n";
				$read = 1;
			}
			my $start_count = 0;
			open(IN,$metrics) || die "Cannot open $metrics metrics file:$!\n\n";
			while(<IN>) {
				my $line = $_;
				chomp($line);
				if ($line =~ /^CYCLE/) {
					$start_count = 1;
					next;
				}
				if ($start_count != 1) {
					next;
				} else {
					if ($line =~ /^(\d+)\s+(\S+)/) {	
						my ($cycle, $qual) = ($1, $2);
						$qual_by_cycle_r1{$id}{$cycle} = $qual if ($read == 1);
						$qual_by_cycle_r2{$id}{$cycle} = $qual if ($read == 2);
						$count++;
						$sum += $qual;
						$max_read_ln_r1 = $count if ($read == 1 && $count > $max_read_ln_r1);
						$max_read_ln_r2 = $count if ($read == 2 && $count > $max_read_ln_r2);
					}
				}	
			}
			close IN;
			$input{$id}{'r1_mean_qual'} =  sprintf("%.1f",($sum/$count)) if  ($read == 1);
			$input{$id}{'r2_mean_qual'} = sprintf("%.1f",($sum/$count)) if  ($read == 2);
			#print "AVG $read $input{$id}{'r1_mean_qual'}  $input{$id}{'r2_mean_qual'}\n";

		}
	}
}

sub getQualityHisto {
	foreach my $id (keys %input) {
		my $total_count = 0;
		foreach my $metrics ($input{$id}{'r1_qual_hist_metrics'}, $input{$id}{'r2_qual_hist_metrics'}) {
			my ($read, $sum, $qual, $q_count);
			my $q_2_sum  = 0 ;
			my $q_10_sum = 0 ;
			next if (($metrics =~ /\.r2\.bam/ || $metrics =~ /na.quality_histo.metrics/) && $input{$id}{'paired'} ne "true" );
			if ($metrics =~ /\.r1\.bam/) {
				$read = 1;
			} elsif ($metrics =~ /\.r2\.bam/) {
				$read = 2;
			} else {					# unpaired, so make this read1
				#print "CHECK $metrics\n";
				$read = 1;
			}
			my $start_count = 0;
			open(IN,$metrics) || die "Cannot open $metrics metrics file:$!\n\n";
			while(<IN>) {
				my $line = $_;
				chomp($line);
				if ($line =~ /^QUALITY/) {
					$start_count = 1;
					next;
				}
				if ($start_count != 1) {
					next;
				} else {
					if ($line =~ /^(\d+)\s+(\S+)/) {	
						($qual, $q_count) = ($1, $2);
						$qual_histo_r1{$id}{$qual} = $q_count if ($read == 1);
						$qual_histo_r2{$id}{$qual} = $q_count if ($read == 2);
						$sum      += $q_count;
						$q_2_sum  += $q_count if ($qual == 2);
						$q_10_sum += $q_count if ($qual < 11);
						$max_qual_score = $qual if ($qual > $max_qual_score);
						$total_count += $q_count;
					}
				}	
			}
			close IN;
			$input{$id}{'r1_pct_q2'}  = sprintf("%.1f",(($q_2_sum/$sum)*100)) if  ($read == 1);
			$input{$id}{'r2_pct_q2'}  = sprintf("%.1f",(($q_2_sum/$sum)*100)) if  ($read == 2);
			$input{$id}{'r1_pct_q10'} = sprintf("%.1f",(($q_10_sum/$sum)*100)) if  ($read == 1);
			$input{$id}{'r2_pct_q10'} = sprintf("%.1f",(($q_10_sum/$sum)*100)) if  ($read == 2);
			#print "AVG $read $input{$id}{'r1_pct_q20'}  $input{$id}{'r2_pct_q20'}\n";
			###########################  Now make all histo values be fraction of total bases  #########################
			for (my $i = 0; $i <= $qual; $i++) { 
				if ($read == 1) {
					if (exists $qual_histo_r1{$id}{$i} ) {
						$qual_histo_r1{$id}{$i} = $qual_histo_r1{$id}{$i} / $total_count;
					} else {
						$qual_histo_r1{$id}{$i} = 0;
					}	
				} else {
					if (exists $qual_histo_r2{$id}{$i} ) {
						$qual_histo_r2{$id}{$i} = $qual_histo_r2{$id}{$i} / $total_count;
					} else {
						$qual_histo_r2{$id}{$i} = 0;
					}	

				}
			}

		}
	}
}

sub getGcHisto {
	foreach my $id (keys %input) {
		foreach my $metrics ($input{$id}{'r1_gc_metrics'}, $input{$id}{'r2_gc_metrics'}) {
			$input{$id}{'r2_mean_gc'} = "na" if $input{$id}{'paired'} eq "false";
			my ($read, $sum, $total_reads, $avg);
			next if (($metrics =~ /\.r2\.bam/ || $metrics =~ /na.read_gc.metrics/) && $input{$id}{'paired'} ne "true" );
			if ($metrics =~ /\.r1\.bam/) {
				$read = 1;
			} elsif ($metrics =~ /\.r2\.bam/) {
				$read = 2;
			} else {					# unpaired, so make this read1
				#print "CHECK $metrics\n";
				$read = 1;
			}
			my $start_count = 0;
			open(IN,$metrics) || die "Cannot open $metrics metrics file:$!\n\n";
			while(<IN>) {
				my $line = $_;
				chomp($line);
				if ($line =~ /^GC/) {
					$start_count = 1;
					next;
				}
				if ($start_count != 1) {
					next;
				} else {
					if ($line =~ /^(\d+)\s+(\S+)/) {	
						my ($gc, $num) = ($1, $2);
						$gc_histo_r1{$id}{$gc} = $num if ($read == 1);
						$gc_histo_r2{$id}{$gc} = $num if ($read == 2);
						$total_reads += $num;
						$sum += ($gc*$num);
					}
				}	
			}
			close IN;
			for  (my $i = -1; $i <= 100; $i++) {                          # replace with fraction of total, start at -1
				if ($read == 1) {
					if (exists $gc_histo_r1{$id}{$i}  && $total_reads != 0) {
						$gc_histo_r1{$id}{$i} = $gc_histo_r1{$id}{$i} / $total_reads;
					} else {
						$gc_histo_r1{$id}{$i} = 0;
					}
				} else {
					if (exists $gc_histo_r2{$id}{$i}  && $total_reads != 0) {
						$gc_histo_r2{$id}{$i} = $gc_histo_r2{$id}{$i} / $total_reads;
					} else {
						$gc_histo_r2{$id}{$i} = 0;
					}
				}
			}
			$input{$id}{'r1_mean_gc'} = sprintf("%.1f",($sum/$total_reads)) if  ($read == 1);
			$input{$id}{'r2_mean_gc'} = sprintf("%.1f",($sum/$total_reads)) if  ($read == 2);

		}
	}
}

sub getComplexityMetrics {
	foreach my $id (keys %input) {
		my $metrics = $input{$id}{'paired_comp_metrics'};  
		my ($est_lib_size, $last_line) = (" - ", "", "");
		my $pct_duplication = " - ";
		my ($dup_column, $lib_size_column);
		if  ($input{$id}{'paired'} eq "false") {
			$input{$id}{'pct_duplication'} = "na";
			$input{$id}{'estimate_library_size'} = "na";
			next;
		}
		my $start_histo_count = 0;
		open(IN,$metrics) || warn "Cannot open $metrics metrics file:$!\n\n";	# change from die to warn 
		while(<IN>) {
			my $line = $_;
			chomp($line);
			my @line = split(/\t/,$line);
			if ($line =~ /^LIBRARY/) {
				print "TEST $last_line\n";
				$dup_column      = &getPicardColumn(\@line,"PERCENT_DUPLICATION");
				$lib_size_column = &getPicardColumn(\@line,"ESTIMATED_LIBRARY_SIZE");
				$last_line = $line;
			} elsif ($last_line =~ /^LIBRARY/) {
				$pct_duplication = $line[$dup_column] || 0;
				$est_lib_size    = $line[$lib_size_column] || 0;# found cases of no value for this field
				$last_line = $line;
			} else {
				$last_line = $line;
			}
			if ($line =~ /^duplication_group_count/) {
				$start_histo_count = 1;
				next;
			}
			if ($start_histo_count != 1) {
				next;
			} else {
				if ($line =~ /^(\d+)\s+(\S+)/) {	
					my ($dup_count, $num) = ($1, $2);
					$complexity_histo{$id}{$dup_count} = $num;
					$max_complexity_count = $dup_count if ($dup_count > $max_complexity_count);
				}
			}	
		}
		close IN;
		$max_complexity_count++	;			# increment this by one in case no duplicate reads 	
		$input{$id}{'pct_duplication'} = sprintf("%.3f",(100*$pct_duplication)) if ($pct_duplication) ne " - ";
		$input{$id}{'estimate_library_size'} = "$est_lib_size";
	}
}

sub getJumpMetrics {
	foreach my $id (keys %input) {
		my $last_line = " - ";
		my $metrics = $input{$id}{'jump_metrics'};  
		my ($pct_jump, $pct_nonjump, $pct_chimera, $jump_mean) = (" - ", " - ", " - ", " - ");
		my ($pct_jump_column, $pct_nonjump_column, $pct_chimera_column, $jump_mean_column);
		if  ($input{$id}{'paired'} eq "false") {
			$input{$id}{'pct_jump'}    = "na";
			$input{$id}{'pct_nonjump'} = "na";
			$input{$id}{'pct_chimera'} = "na";
			$input{$id}{'jump_mean'}   = "na";
			next;
		}
		open(IN,$metrics) || warn "Cannot open $metrics metrics file:$!\n\n";	# change from die to warn 
		while(<IN>) {
			my $line = $_;
			chomp($line);
			my @line = split(/\t/,$line);
			if ($line =~ /^JUMP_PAIRS/) {
				$pct_jump_column      = &getPicardColumn(\@line,"PCT_JUMPS");
				$pct_nonjump_column   = &getPicardColumn(\@line,"PCT_NONJUMPS");
				$pct_chimera_column   = &getPicardColumn(\@line,"PCT_CHIMERAS");
				$jump_mean_column     = &getPicardColumn(\@line,"JUMP_MEAN_INSERT_SIZE");
				$last_line = $line;
			} elsif ($last_line =~ /^JUMP_PAIRS/) {
				$pct_jump    = $line[$pct_jump_column]    || 0;
				$pct_nonjump = $line[$pct_nonjump_column] || 0;
				$pct_chimera = $line[$pct_chimera_column] || 0;
				$jump_mean   = $line[$jump_mean_column]   || 0;
				$last_line = $line;
			} else {
				$last_line = $line;
			}
		}
		close IN;
		$input{$id}{'pct_jump'}    = sprintf("%.1f",(100*$pct_jump))    if ($pct_jump) ne " - ";
		$input{$id}{'pct_nonjump'} = sprintf("%.1f",(100*$pct_nonjump)) if ($pct_nonjump) ne " - ";
		$input{$id}{'pct_chimera'} = sprintf("%.1f",(100*$pct_chimera)) if ($pct_chimera) ne " - ";
		$input{$id}{'jump_mean'}   = sprintf("%.1f",$jump_mean)         if ($jump_mean) ne " - ";
	}
}

#sub plotRarefactionData {
#	my (@chartobj);
#	my ($data_test, $data_test_max) = (-1, -1);
#	foreach my $id (keys %input) {
#		next if ($input{$id}{'aligned'} eq "false");
#		my $metrics = $input{$id}{'complexity_by_start_metrics'};
#		$plot_complexity_by_start = 1;
#		my(@data, @x_point, @legends, @values, $max_x, $max_y);
#		open(IN,$metrics) || die "Cannot open $metrics metrics file:$!\n\n";
#		while(<IN>) {
#			my $line = $_;
#			chomp($line);
#			my @line = split(/\s+/,$line);
#			next if ($line =~ /^[Bam|Reads|Pairs]/);
#			$complexity_by_start_histo{$id}{$line[1]} = $line[2];
#			#$data[$line[1]] = $line[2];
#			push(@data,$line[2]);
#			push(@x_point,$line[1]);
#		}
#		close IN;
#		push(@legends,$id);
#		push(@values,\@data);
#                $data_test = $#data ;
#                $data_test_max = $data_test if ($data_test > $data_test_max); 
#		my $output = "$qc_chart_data_subdir/$id" . "_rarefaction_data.eps";
#		$output =~ s/\./\_/g;				# pdflatex fails with any(or too many) dots in filename
#		$output =~ s/\_eps$/\.eps/;
#		my $short_id = substr($id,0,10);
#		my $title  = "Complexity by Start Position Rarefaction Plot - $short_id";
#        	my $chartobj = new AAGnuplot('title' => "$title",
#                                 'xlabel' => "Pairs/Reads Sequenced",
#                                 'ylabel' => "Pairs/Reads Aligned",
#                                 #'xrange' => "[1:$max_x]",
#                                 #'yrange' => "[0:$max_y]",
#                                 'output' => "$output",
#                                 'graph' => "lines");
#                                 #'graph' => "x:y");
#        	$chartobj->addOpts('size' => "1.4,1.4");
#        	$chartobj->plotAoA('data' => \@values,
#                                'xlabels' => \@x_point,
#                                'legend labels' => \@legends);
#		push(@chartobj, $chartobj);
#	}
#	if ($data_test_max == -1) {
#		print "No data to print\n";
#        	return("NA");                                 # this is to avoid attempt to plot nothing.
#	} else {
#        	return(@chartobj);
#
#	}
#}

sub buildAdapterCommand {
	my $in_file = $_[0];
	my ($id, $string);
	my %adapter;
	open(IN,$in_file) ||  die "Cannot open $in_file adapter fasta file:  $!\n\n";
	while(<IN>) {
		my $line = $_;
		chomp($line);
		if ($line =~ /^>(\S+)/) {
			$id = $1;
		} else {
			$adapter{$id} .= $line;
		}
	}
	close IN;
	foreach my $adapter_id (keys %adapter) {
		$string .= " ADAPTER_SEQUENCE=$adapter{$adapter_id}";
	}
	return($string);

}

sub getAdapterMetrics {
	foreach my $id (keys %input) {
		my $metrics = $input{$id}{'asm_metrics'};
		next if (! -e $metrics);
		my ($adapter_column, $r1_adapter_pct, $r2_adapter_pct, $paired_adapter_pct) = (-1, "na", "na", "na");
		open(ASM,$metrics) || die "Cannot open $metrics metrics file:$!\n\n";
		while(<ASM>) {
			my $line = $_;
			chomp($line);
			my @line = split(/\s+/,$line);
			if ($line =~ /PCT_ADAPTER/) {
				$adapter_column = &getPicardColumn(\@line,"PCT_ADAPTER");
				next if ($adapter_column == -1);
			} elsif ($line =~ /FIRST_OF_PAIR/) {
				$line[$adapter_column] = "0.000" if ($line[$adapter_column] == "0");
				$r1_adapter_pct = $line[$adapter_column];
			} elsif ($line =~ /SECOND_OF_PAIR/) {
				$line[$adapter_column] = "0.000" if ($line[$adapter_column] == 0);
				$r2_adapter_pct = $line[$adapter_column];
			} elsif ($line =~ /^PAIR\s+/) {
				$line[$adapter_column] = "0.000" if ($line[$adapter_column] == 0);
				$paired_adapter_pct = $line[$adapter_column];
			} elsif ($line =~ /UNPAIRED/) {
				$line[$adapter_column] = "0.000" if ($line[$adapter_column] == 0);
				$r1_adapter_pct = $line[$adapter_column];
			}	
		}
		close ASM;

		# picard defines adapter as:  %  PF reads unaligned & match to known adapter seqs right from the start of the read #
		$input{$id}{'r1_pct_adapter'} = $r1_adapter_pct*100 	if ($r1_adapter_pct ne "na");	
		$input{$id}{'r2_pct_adapter'} = $r2_adapter_pct*100	if ($r2_adapter_pct ne "na"); 
		$input{$id}{'pct_adapter'}    = $paired_adapter_pct*100 if ($paired_adapter_pct ne "na");
	}

}

sub getInsertSizeMetrics {
	foreach my $id (keys %input) {
		my $metrics  = $input{$id}{'paired_insert_size_metrics'};
		my ($fr_mean, $fr_median, $rf_mean, $rf_median) = ("", "", "", "");
		my ($mean_column, $median_column, $orientation_column, $start_check, $at_histo) = (-1, -1, -1, -1, -1);
		if  ($input{$id}{'aligned'} eq "false" || $input{$id}{'paired'} eq "false") {
			$input{$id}{'fr_mean_insert'}   = "na";
			$input{$id}{'rf_mean_insert'}   = "na";
			$input{$id}{'fr_median_insert'} = "na";
			$input{$id}{'rf_median_insert'} = "na";
			next;
		}
		open(ISM,$metrics) || die "Cannot open $metrics metrics file:$!\n\n";
		while(<ISM>) {
			my $line = $_;
			chomp($line);
			my @line = split(/\s+/,$line);
			if ($line =~ /MEAN_INSERT_SIZE/) {
				($mean_column)        = &getPicardColumn(\@line,"MEAN_INSERT_SIZE");
				($median_column)      = &getPicardColumn(\@line,"MEDIAN_INSERT_SIZE");
				($orientation_column) = &getPicardColumn(\@line,"PAIR_ORIENTATION");
			} elsif ($line =~ /\s+FR\s+/ || $line =~ /\s+RF\s+/) {
				$input{$id}{'fr_mean_insert'}     = sprintf("%.0f",($line[$mean_column]))  if $line[$orientation_column] eq "FR";
				$input{$id}{'fr_median_insert'}   = $line[$median_column]  if $line[$orientation_column] eq "FR";
				$input{$id}{'rf_mean_insert'}     = sprintf("%.0f",($line[$mean_column]))  if $line[$orientation_column] eq "RF";
				$input{$id}{'rf_median_insert'}   = $line[$median_column] if $line[$orientation_column] eq "RF";
			} else {
				next;
			}
		}
		close ISM;
	}

}

sub getInsertSizeHisto {
	my ($max_fr, $max_rf, $rf_size, $fr_size) = (0, 0, 0, 0);
	foreach my $id (keys %input) {
		my $metrics_file  = $input{$id}{'paired_insert_size_metrics'};
        	my $start_count = 0;
        	my ($total_rf_inserts, $total_fr_inserts) = (0,0);
        	my ($get_fr, $get_rf) = (-1, -1);
		my $size = 0;
		next if ($input{$id}{'paired'} eq "false" || $input{$id}{'aligned'} eq "false");
        	open(ISM,"$metrics_file") || die "Cannot open the $metrics_file insert_size_metrics file:  $!\n\n";
        	while(<ISM>) {
                	my $line = $_;
                	chomp($line);
                	my $line_num = $.;                      
                	my @line = split(/\t/, $line);
                	if ($line =~ /^insert_size/) {
                       		$start_count = 1;
                        	if ($#line == 1) {                      # if only 2 columns then get second column;
                               		$get_fr = 1;
                        	} else {
                                	foreach my $i (1 .. $#line) {                           # which column is correct direction 
                                       		$get_fr = $i if ($line[$i] =~ /fr_count/);
                                        	$get_rf = $i if ($line[$i] =~ /rf_count/);
                                	}
                        	}
                	}
                	if ($start_count !=  1) {
                       		next;
                	} elsif ($start_count == 1 && $line =~ /^\d/) {
                        	$size = $line[0];
                                #$insert_size_histo_rf{$id}{$size} = $line[$get_rf] if ($line[$get_rf] && $get_rf > -1);
                                #$insert_size_histo_fr{$id}{$size} = $line[$get_fr] if ($line[$get_fr] && $get_fr > -1) ;
                                $insert_size_histo_rf{$id}{$size} = $line[$get_rf] if ($get_rf > -1);
                                $insert_size_histo_fr{$id}{$size} = $line[$get_fr] if ($get_fr > -1) ;
                                $total_rf_inserts                 += $insert_size_histo_rf{$id}{$size} if ($insert_size_histo_rf{$id}{$size});
                                $total_fr_inserts                 += $insert_size_histo_fr{$id}{$size} if ($insert_size_histo_fr{$id}{$size});
                	} else {                                               # this should be last line in file (blank)
                        	next;
                	}
        	}
        	close ISM;
        	##############  Switch the value from number of inserts to fraction of total aligned (fr+rf) inserts ##############
        	for (my $j = 0 ; $j <= $size; $j++) {
                        if (exists $insert_size_histo_rf{$id}{$j}) {
                               	$insert_size_histo_rf{$id}{$j} = $insert_size_histo_rf{$id}{$j} / ($total_rf_inserts + $total_fr_inserts);
				$rf_size = $j;
                        } else {
                               	$insert_size_histo_rf{$id}{$j} = 0;
                        }
                        if (exists $insert_size_histo_fr{$id}{$j}) {
                               	$insert_size_histo_fr{$id}{$j} = $insert_size_histo_fr{$id}{$j} / ($total_fr_inserts + $total_rf_inserts);
				$fr_size = $j;
                        } else {
                               	$insert_size_histo_fr{$id}{$j} = 0;
                        }
        	}
        	#################################  Determine appropriate cutoff to plot  #######################################
        	my $mean_fr = $input{$id}{'fr_mean_insert'} || -1;
        	my $mean_rf = $input{$id}{'rf_mean_insert'} || -1;
        	if ($rf_size > ($mean_rf*3)) {                      # so not plotting size of upwards of 1 MB for chimerics
               		$rf_size = int($mean_rf*3);
			$max_rf = $rf_size if ($rf_size > $max_rf);

		} else {
			$max_rf = $rf_size ;
		}
        	if ($fr_size > ($mean_fr*3)) {                      # so not plotting size of upwards of 1 MB for chimerics
               		$fr_size = int($mean_fr*3);
			$max_fr = $fr_size if ($fr_size > $max_fr);
        	} else {
			$max_fr = $fr_size;
		}
	}
        return($max_fr, $max_rf);                                                  # return size to later compare for the max size
}

sub getCoverageMetrics {
	foreach my $id (keys %input) {
		my ($sum, $win_count) = (0, 0, 0);
		my $metrics_file  = $input{$id}{'seq_cvg_metrics'};
        	my $start_count = 0;
		my $ref = $input{$id}{'ref'};
		my ($win, $cvg, $rounded_cvg);
		next if ($input{$id}{'aligned'} eq "false");
		print "Collecting coverage metrics for $id....\n\n";
        	open(CVG,"$metrics_file") || die "Cannot open the $metrics_file coverage metrics file:  $!\n\n";
        	while(<CVG>) {
                	my $line = $_;
                	chomp($line);
                	my $line_num = $.;                      
                	my @line = split(/\t/, $line);
                	if ($line =~ /^Ref_Window/) {
                       		$start_count = 1;
			}
                	if ($start_count !=  1) {
                       		next;
                	} elsif ($start_count == 1 && $line =~ /^\d+/) {
                        	($win,$cvg) = ($line[0], $line[1]);
                                $coverage_plot{$id}{$win_count} = $cvg;
				$rounded_cvg = sprintf("%.0f",$cvg);
				$coverage_histo{$id}{$rounded_cvg}++;
				$max_histo_coverage = $rounded_cvg if ($rounded_cvg > $max_histo_coverage);
				$window_hash{$id}{$win_count} = $win;
				$win_count++;
				$sum += $cvg;
				$references{$ref}{'max_ref_cvg_win'} = $win_count;
                	} else {                                                # this should be last line in file (blank)
                        	next;
			}
		}
		close CVG;
		$input{$id}{'mean_cvg'} = sprintf("%.0f",($sum/$win_count));
	}

}

sub runLorenzCurve {
	foreach my $ref (keys %references) {
		my ($ref_file, $bam_list, $out, $cmd, $metrics, $orig_png, $png);
		$ref_file = $references{$ref}{'reffile'};
		$out = "$qc_chart_data_subdir/$ref";
		$out =~ s/\./\_/g;						# pdflatex can't handle extensions well
		$orig_png = $out . ".lorenz_plot.png";
		$png = $out . "_lorenz_plot.png";
		foreach my $lib (keys %input) {
			next if	($input{$lib}{'ref'} ne $ref); 
			if ($input{$lib}{'aligned'} eq "true") {
				$bam_list .= "-b $input{$lib}{'bam'} ";
			} else {
				$bam_list .= "-b $input{$lib}{'r1bam'} ";
			}
		}
		$cmd = "$lorenz_tool $bam_list -r $ref_file -o $out";		
		print "Running library_comparision.py to generate Lorenz curve...\n";
		if ($no_overwrite == 1 && -e $png) {
			print "Lorenz curve plot exists for reference and -N used.  Skipping library_comparision.py...\n\n";
		} else {
			print "$cmd\n\n";
			system "$cmd";
			`mv $orig_png $png` if ($orig_png);			# because pdflatex can't handled extensions well
		}
		$references{$ref}{'lorenz_plot'} = $png;
	}
}

sub readBlast {
	foreach my $id (keys %input) {
		my ($subset_fasta_r1, $subset_fasta_r2, $subset_fasta_paired, $subset_r1_cmd, $subset_r2_cmd, $subset_pair_cmd);
		my ($read_count, $tax_out_paired, $tax_out_r1, $tax_out_r2);
		my ($bam1, $paired_bam) = ($input{$id}{'r1bam'}, $input{$id}{'bam'});
		my $rd_ln = $input{$id}{'r1_mean_rd_ln'};
		my $out_subset_bam1       = $bam1;
		my $out_paired_subset_bam = $paired_bam;

		if ($input{$id}{'paired'} eq "true" && $paired_bam) {
			$out_paired_subset_bam =~ s/.*\///;
			$out_paired_subset_bam =~ s/^/read_qc\/read_data\//;
			$out_paired_subset_bam =~ s/\.bam/\.$subset\.bam/;
			$subset_pair_cmd = "java -jar $aa_picard/SubsetSam.jar I= $paired_bam O= $out_paired_subset_bam R=$subset";
			system "$subset_pair_cmd" unless ($no_overwrite && -e "$out_paired_subset_bam");

			$subset_fasta_paired = &bamToFasta($out_paired_subset_bam) if  (-e $out_paired_subset_bam);	
			$tax_out_paired = &runBlastAndTax($subset_fasta_paired,$rd_ln);
			&parseTax($tax_out_paired, $subset_fasta_paired, $id);		# give fasta to get total rd count 
		} else {
			$out_subset_bam1 =~ s/.*\///;
			$out_subset_bam1 =~ s/^/read_qc\/read_data\//;
			$out_subset_bam1 =~ s/\.bam/\.$subset\.bam/;
			### SubsetSam.jar does not work on paired reads - use DownsampleSam  ###
			#$subset_r1_cmd = "java -jar $aa_picard/SubsetSam.jar I= $bam1 O= $out_subset_bam1 R=$subset";
			$subset_r1_cmd = "java -jar $picard/picard.jar DownsampleSam I=$bam1 O= $out_subset_bam1 R=$subset";
			unless ($no_overwrite && -e "$out_subset_bam1") {
				print "Subsetting $bam1\n$subset_r1_cmd\n";	
				system "$subset_r1_cmd" if (-e $bam1);
			}
			$subset_fasta_r1 = &bamToFasta($out_subset_bam1) if  (-e $out_subset_bam1);
			$tax_out_r1 = &runBlastAndTax($subset_fasta_r1,$rd_ln);
			&parseTax($tax_out_r1, $subset_fasta_r1, $id);		

		}
	}

}

sub bamToFasta {
	my $bam = $_[0];
	my $fasta;
	$fasta = $bam;
	$fasta =~ s/\.bam/\.fasta/;
	open(OUT,">$fasta") || warn "Cannot write to the $fasta:  $!\n\n";
	open(IN,"$samtools view $bam | ") || warn "Cannot open the $bam file:  $!\n\n";
	while(<IN>) {
		my $line = $_;
		chomp($line);
		my @line = split(/\s+/,$line);
		print OUT ">$line[0]\n$line[9]\n";
	}
	close IN;
	close OUT;
	return $fasta;
}

sub runBlastAndTax {
	my ($fasta, $rd_ln) = ($_[0], $_[1]);
	my $blast_out = $fasta;
	$blast_out    =~ s/\.fasta/\.blast/;
	my $tax_out   = $blast_out . ".tax.out";
	my $task;
	#if ($rd_ln < $short_blast_rd_ln) {		# blastn-short not practical, back to dc-megablast for all reads
	#	$task = "blastn-short";
	#} else {
		$task = "dc-megablast";
	#}
	
	my $blast_cmd = "$blastn -query $fasta -db /broad/data/blastdb/nt/nt -out $blast_out -evalue 1e-10  -task $task -outfmt 6 -num_alignments 1 -num_descriptions 1 -max_target_seqs 1 -num_threads $threads";
	my $tax_cmd   = "$tax_from_blast -i $blast_out -o $tax_out  -k -y -g -s -m 50 -l";

	if ($no_overwrite && -e $blast_out) {
		print "The -N no_overwrite option set and blast output exists ... skipping read blast ...\n\n";
	} else {
		print "Running read blast using $task for $rd_ln length reads...\n\n";
		system "$blast_cmd";
	}
	if ($no_overwrite && -e $tax_out) {
		print "The -N no_overwrite option set and taxonomy from blast output exists ... skipping taxonomy from blast ...\n\n";
	} else {
		print "\nRunning taxonomy from blast...\n\n";
		system "$tax_cmd" unless (-z $blast_out);		# do not run tax from blast if blast out empty	
	}
	return $tax_out;
}

sub parseTax {
	my ($tax_out, $fasta, $id) = ($_[0], $_[1], $_[2]);
	my %read_blast;
	my ($hit_count, $read_count, $nohit_count, $fraction, $tax_sum_out) = (0, 0, 0, 0, "");
	$tax_sum_out = $fasta;
	$tax_sum_out =~ s/\.fasta/\.blast\.tax\.summary\.out/;
	$read_count = `grep -c ">" $fasta`;
	chomp($read_count);
	if (-e $tax_out) {
        	open(TAX,"$tax_out") || warn "Cannot open $tax_out taxonomy from blast file:  $!\n\n";
        	while(<TAX>){
               		my $line = $_;
        		chomp($line);
                	my @line = split(/\t+/,$line);
               		$hit_count++;
                	$read_blast{$line[3]}++;
		}
	}
        close TAX;
        $nohit_count = $read_count - $hit_count;
        $read_blast{"No_Hit"} = $nohit_count;
        print "Parsing taxonomy_from_blast.pl output ...\n\n";
        open(SUM,">$tax_sum_out") || die "Can't write to $tax_sum_out: $!\n\n";
        print SUM "Reads\tPct\tHit\n";
        foreach my $i (keys %read_blast) {
                if ($read_count > 0 ) {
                        $fraction = sprintf("%.3f",($read_blast{$i} / $read_count)*100);
                } else {
                        $fraction = " - ";
                }
                print SUM "$read_blast{$i}\t$fraction\t$i\n";
                $orgs_rd_blast{$i}{$id} = $fraction;
                $org_count{$i} += $read_blast{$i};
        }
        close SUM;

}

sub prepareBlastTable {
	my $org_count = 0;
	my $lib_count = 0;
	push(@read_blast_headers,"Organism");
	foreach my $i (keys %input) {
		$lib_count++;
		my $header = "Library " . $lib_count;
		push(@read_blast_headers,$header);
	}
	foreach my $j (sort { $org_count{$b} <=> $org_count{$a} } keys %org_count) {
       		my @blast_push;
        	push(@org_list,$j);
        	push(@blast_push,$j);
        	$org_count++;
        	last if ($org_count > 35);
        	foreach my $k (keys %input) {
                	my $val;
                	if (exists $orgs_rd_blast{$j}{$k}) {
                       	 	$val = $orgs_rd_blast{$j}{$k} . "%";
                	} else {
                       	 	$val = "0.0%";
                	}
                	push(@blast_push,$val);
	        }
        	push(@read_blast_data,[@blast_push]);
		print "\n";
	}
}

sub getGcBiasHisto {
	foreach my $id (keys %input) {
		my $metrics_file  = $input{$id}{'gc_bias_metrics'};
        	my $start_count = 0;
		my $ref = $input{$id}{'ref'};
		my ($gc, $cvg, $ref_win, $sum_ref_win, $gc_column, $win_column, $cvg_column);
		print "Collecting GC bias metrics for $id....\n\n";
		if ($input{$id}{'aligned'} eq "false" || ! -e $metrics_file || ! -e $input{$id}{'reffile'}) {
			print "Unable to run GC bias - either not aligned ref not found or metrics not found...\n\n";
			next;
		}
		($references{$ref}{'max_gc_bias_cvg'}, $references{$ref}{'max_gc_bias_ref_win'}) = (0, 0);
        	open(CVG,"$metrics_file") || die "Cannot open the $metrics_file coverage metrics file:  $!\n\n";
        	while(<CVG>) {
                	my $line = $_;
                	chomp($line);
                	my $line_num = $.;                      
                	my @line = split(/\t/, $line);
                	if ($line =~ /READ_STARTS/) {
                       		$start_count = 1;
				$gc_column  = &getPicardColumn(\@line,"GC");
				$win_column = &getPicardColumn(\@line,"WINDOWS");
				$cvg_column = &getPicardColumn(\@line,"NORMALIZED_COVERAGE");
			}
                	if ($start_count !=  1) {
                       		next;
                	} elsif ($start_count == 1 && ($line =~ /^All\s+Reads/ || $line =~ /^\d+/)) {
                        	($gc,$ref_win, $cvg) = ($line[$gc_column],  $line[$win_column], $line[$cvg_column]);
                                $gc_bias_histo{$id}{$gc} = $cvg;
				$gc_bias_ref_windows{$ref}{$gc} = $ref_win || 0 ;
				$sum_ref_win += $ref_win;
				$references{$ref}{'max_gc_bias_cvg'} = $cvg if ($cvg > $references{$ref}{'max_gc_bias_cvg'});
                	} else {                                                # this should be last line in file (blank)
                        	next;
			}
		}
		close CVG;
		#####################  update the reference windows to fraction of total windows  #######################
		for (my $i = 0 ;$i < 101 ; $i++	) {
			my $orig = $gc_bias_ref_windows{$ref}{$i};
			my $newval = $gc_bias_ref_windows{$ref}{$i} || 0;	
			$gc_bias_ref_windows{$ref}{$i} = ($newval / $sum_ref_win);
			$references{$ref}{'max_gc_bias_ref_win'} = $gc_bias_ref_windows{$ref}{$i} if ($gc_bias_ref_windows{$ref}{$i} > $references{$ref}{'max_gc_bias_ref_win'});
		}	
	}
}


sub plotInsertSizeHisto {
	my ($dir, $hash_ref, $max_x) = ($_[0], $_[1], $_[2]);
	my (@x_point, @legends, @values, @data, $k, $legend);
	my ($data_test, $data_test_max) = (-1, -1);                             # if there is no data to plot, gnuplot hangs
	my $max_y = 0;
	my %histo = %$hash_ref;
	my $output = "$qc_chart_data_subdir/insert_size_metric_histogram_" . $dir;
	print "Plotting $dir insert size metrics histogram...\n\n";
        foreach my $lib (keys %input) {
		next if ($input{$lib}{'aligned'} eq "false" || $input{$lib}{'paired'} eq "false");
                $legend = $lib;
                $legend =~ s/\s+/\_/g;
                my @data;       
                for (my $j=1; $j <= $max_x; $j++) {
                        $data[$j] = $histo{$lib}{$j} if (defined $histo{$lib}{$j});
			push(@x_point,$j);
			next if (! $data[$j]);
			$max_y = $data[$j] if ($data[$j] > $max_y);
                }
		shift(@data);							# to fix first value not at base position 1
                push(@legends,$legend);
                push(@values,\@data);
                $data_test = $#data ;
                $data_test_max = $data_test if ($data_test > $data_test_max); 
        }
        return("NA") if ($data_test_max == -1);                                 # this is to avoid attempt to plot nothing.
        my $chartobj = new AAGnuplot('title' => "Insert Size Histogram - Pairs Aligned in $dir Orientation ",
                                 'xlabel' => "Insert Size (bp)",
                                 'ylabel' => "Fraction of Aligned $dir Pairs",
                                 'xrange' => "[1:$max_x]",
                                 'yrange' => "[0:$max_y]",
                                 'output' => "$output.eps",
                                 'graph' => "lines");
                                 #'graph' => "x:y");
        $chartobj->addOpts('size' => "1.4,1.4");
        $chartobj->plotAoA('data' => \@values,
                                'xlabels' => \@x_point,
                                'legend labels' => \@legends);
        return($chartobj);
}

sub plotCoverageHisto {
	my ($hash_ref, $max_x) = ($_[0], $_[1]);
	my (@x_point, @legends, @values, @data, $k, $legend);
	my ($data_test, $data_test_max) = (-1, -1);                             # if there is no data to plot, gnuplot hangs
	my $max_y = 0;
	my %histo = %$hash_ref;
	my $output = "$qc_chart_data_subdir/coverage_histogram";
	print "Plotting coverage histogram...\n\n";
        foreach my $lib (keys %input) {
		next if ($input{$lib}{'aligned'} eq "false");
                $legend = $lib;
                $legend =~ s/\s+/\_/g;
                my @data;       
                for (my $j=0; $j <= $max_x+1; $j++) {
                        #$data[$j] = $histo{$lib}{$j} if (defined $histo{$lib}{$j});
                        $data[$j] = $histo{$lib}{$j} || 0;
			push(@x_point,$j);
			$max_y = $data[$j] if ($data[$j] > $max_y);
                }
		shift(@data);							# to fix first value not at base position 1
                push(@legends,$legend);
                push(@values,\@data);
                $data_test = $#data ;
                $data_test_max = $data_test if ($data_test > $data_test_max); 
        }
        return("NA") if ($data_test_max == -1);                                 # this is to avoid attempt to plot nothing.
        my $chartobj = new AAGnuplot('title' => "Aligned Coverage Histogram ",
                                 'xlabel' => "Depth of Coverage",
                                 'ylabel' => "Number of Bases",
                                 'xrange' => "[0:$max_x]",
                                 'yrange' => "[0:$max_y]",
                                 'output' => "$output.eps",
                                 'graph' => "lines");
                                 #'graph' => "x:y");
        $chartobj->addOpts('size' => "1.4,1.4");
        $chartobj->plotAoA('data' => \@values,
                                'xlabels' => \@x_point,
                                'legend labels' => \@legends);
        return($chartobj);

}

sub plotCoverageAlongRef {
	my ($hash_ref) = ($_[0]);
        my ($maxwin, $ref_name_brief) = (0, 0);
	my @chartobj;
        my ($max_median, $yaxis_max, $median) = (0, 0, 0);
        foreach my $ref (keys %references) {
		$maxwin = $references{$ref}{'max_ref_cvg_win'};
                next if ($maxwin == 0);                                         # if there was no coverage found, go to next file
                my (@x_point, @legends, @values, @data);
                foreach my $lib (keys %input) {
                        next if ($ref ne $input{$lib}{'ref'});
                        my @data;       
                        for (my $j=0; $j <= ($maxwin-1); $j++) {
                                push(@x_point,$window_hash{$lib}{$j});
                                $data[$j] = $coverage_plot{$lib}{$j} if (defined $coverage_plot{$lib}{$j});
                        }
                        push(@legends,$lib);
                        push(@values,\@data);

                }

                $ref_name_brief  = $ref;
                $ref_name_brief =~ s/.*(v\d+)\/(\S+)/$2_$1/;
                $ref_name_brief =~ s/\.fasta//;
                $ref_name_brief =~ s/\.//g;
		my $output = "$qc_chart_data_subdir/$ref_name_brief" . "_sequence_coverage.eps",

                my $title = "Read Coverage Along Reference " . $ref;
                my $chartobj = new AAGnuplot('title' => "$title",
                                         'xlabel' => "Reference Window",
                                         'ylabel' => "Coverage",
                                         'xrange' => "[0:$maxwin]",
                                         'yrange' => "[0:$yaxis_max]",
                                         'output' => "$output",
                                         'graph' => "lines");

                $median = $chartobj->getMedian('data' => \@values);
                my $multipler = $median <=20 ? 5 : 3;
                $yaxis_max = $median * $multipler;
                $chartobj->addOpts('yrange' => "[0:$yaxis_max]");
                $chartobj->addOpts('size' => "1.5,1");
                $chartobj->plotAoA('data' => \@values,
                                        'xlabels' => \@x_point,
                                        'legend labels' => \@legends);
                push(@chartobj, $chartobj);
        }
	return(@chartobj);
}

sub plotGcBiasHisto {
	my ($hash_ref, $ref_win_has_ref) = ($_[0], $_[1]);
	my ($max_gc_bias_cvg, $max_gc_bias_ref_win); 
	my @chartobj;
	my ($ref_win_factor, $ref_name_brief);					
        foreach my $ref (keys %references) {
		next if ( ! -e $references{$ref}{'reffile'} );
		$max_gc_bias_cvg     = $references{$ref}{'max_gc_bias_cvg'};
		$max_gc_bias_ref_win = $references{$ref}{'max_gc_bias_ref_win'};
		$ref_win_factor = 1/$max_gc_bias_ref_win;			# factor to get ref window curve to max at 1.0
                next if ($max_gc_bias_cvg < 1);                                 # if there was no coverage found, go to next file
		my @ref_data;
                my (@x_point, @legends, @values, @data);
		for (my $i=0; $i <= 100; $i++) {
			$ref_data[$i] = ($ref_win_factor*$gc_bias_ref_windows{$ref}{$i}) || 0 ;
		}
                push(@values,\@ref_data);
		push(@legends,"Reference_GC");
                foreach my $lib (keys %input) {
                        next if ($ref ne $input{$lib}{'ref'});
                        my @data;       
                        for (my $j=0; $j <= 100; $j++) {
                                push(@x_point,$j);
                                $data[$j] = $gc_bias_histo{$lib}{$j} || 0;
                        }
                        push(@legends,$lib);
                        push(@values,\@data);

                }
                $ref_name_brief  = $ref;
                $ref_name_brief =~ s/.*(v\d+)\/(\S+)/$2_$1/;
                $ref_name_brief =~ s/\.fasta//;
                $ref_name_brief =~ s/\.//g;

		my $output = "$qc_chart_data_subdir/$ref_name_brief" . "_gc_bias.eps",
                my $title = "GC Bias in Reference Coverage " . $ref;
                my $chartobj = new AAGnuplot('title' => "$title",
                                         'xlabel' => "%GC",
                                         'ylabel' => "Normalized Coverage",
                                         'xrange' => "[0:100]",
                                         'yrange' => "[0:$max_gc_bias_cvg]",
                                         'output' => "$output",
                                         'graph' => "lines");

        	$chartobj->addOpts('size' => "1.4,1.4");
                $chartobj->plotAoA('data' => \@values,
                                        'xlabels' => \@x_point,
                                        'legend labels' => \@legends);
                push(@chartobj, $chartobj);
        }
	return(@chartobj);
}

sub plotQualityByPosition {
	my ($read, $qual_hash_ref, $max_cycle) = ($_[0], $_[1], $_[2]);
	my (@x_point, @legends, @values, @data, $k, $legend);
	my ($data_test, $data_test_max) = (-1, -1);                             # if there is no data to plot, gnuplot hangs
	my %qual_by_cycle = %$qual_hash_ref;
	my $output = "$qc_chart_data_subdir/mean_base_quality_by_position_" . $read;
	print "Plotting $read mean quality by cycle...\n\n";
        foreach my $lib (keys %input) {
		next if ($read eq "read2" && $input{$lib}{'paired'} eq "false");
                $legend = $lib;
                $legend =~ s/\s+/\_/g;
                my @data;       
                for (my $j=1; $j <= $max_cycle; $j++) {
                        $data[$j] = $qual_by_cycle{$lib}{$j} if (defined $qual_by_cycle{$lib}{$j});
			push(@x_point,$j);
                }
		shift(@data);							# to fix first value not at base position 1
                push(@legends,$legend);
                push(@values,\@data);
                $data_test = $#data ;
                $data_test_max = $data_test if ($data_test > $data_test_max); 
        }
        return("NA") if ($data_test_max == -1);                                 # this is to avoid attempt to plot nothing.
        my $chartobj = new AAGnuplot('title' => "Mean Base Quality by Position - $read ",
                                 'xlabel' => "Position in Read (bp)",
                                 'ylabel' => "Quality",
                                 'xrange' => "[1:$max_cycle]",
                                 'yrange' => "[0:40]",
                                 'output' => "$output.eps",
                                 #'graph' => "lines");
                                 'graph' => "x:y");
        $chartobj->addOpts('size' => "1.4,1.4");
        $chartobj->plotAoA('data' => \@values,
                                'xlabels' => \@x_point,
                                'legend labels' => \@legends);
        return($chartobj);
}

sub plotQualityHisto {
	my ($read, $qual_hash_ref, $max_q_score) = ($_[0], $_[1], $_[2]);
	my (@x_point, @legends, @values, @data, $k, $legend);
	my $max_y = 0;
	my ($data_test, $data_test_max) = (-1, -1);                             # if there is no data to plot, gnuplot hangs
	my %qual_histo = %$qual_hash_ref;
	my $output = "$qc_chart_data_subdir/base_quality_histogram_" . $read;
	print "Plotting $read base quality histogram...\n\n";
        foreach my $lib (keys %input) {
		next if ($read eq "read2" && $input{$lib}{'paired'} eq "false");
                $legend = $lib;
                $legend =~ s/\s+/\_/g;
                my @data;       
                for (my $j=1; $j <= $max_q_score; $j++) {
                        #$data[$j] = $qual_histo{$lib}{$j} if (defined $qual_histo{$lib}{$j});
                        $data[$j]  = $qual_histo{$lib}{$j} || 0;
			push(@x_point,$j);
			$max_y = $data[$j] if ($data[$j] > $max_y);
                }
		shift(@data);							# to fix first value not at base position 1
                push(@legends,$legend);
                push(@values,\@data);
                $data_test = $#data ;
                $data_test_max = $data_test if ($data_test > $data_test_max); 
        }
        return("NA") if ($data_test_max == -1);                                 # this is to avoid attempt to plot nothing.
        my $chartobj = new AAGnuplot('title' => "Base Quality Histogram - $read ",
                                 'xlabel' => "Quality Score",
                                 'ylabel' => "Fraction of Total Bases",
                                 'xrange' => "[0:$max_q_score]",
                                 'yrange' => "[0:$max_y]",
                                 'output' => "$output.eps",
                                 #'graph' => "lines");
                                 'graph' => "bar");
        $chartobj->addOpts('size' => "1.4,1.4");
        $chartobj->plotAoA('data' => \@values,
                                'xlabels' => \@x_point,
                                'legend labels' => \@legends);
        return($chartobj);
}

sub plotGcHisto {
	my ($read, $gc_hash_ref) = ($_[0], $_[1]);
	my (@x_point, @legends, @values, @data, $k, $legend);
	my $max_y = 0;								# max y = largest # of reads at a given %GC
	my ($data_test, $data_test_max) = (-1, -1);                             # if there is no data to plot, gnuplot hangs
	my %gc = %$gc_hash_ref;
	my $output = "$qc_chart_data_subdir/gc_histogram_" . $read;
	print "Plotting $read gc histogram...\n\n";
        foreach my $lib (keys %input) {
		next if ($read eq "read2" && $input{$lib}{'paired'} eq "false");
                $legend = $lib;
                $legend =~ s/\s+/\_/g;
                my @data;       
                for (my $j=0; $j <= 100; $j++) {
                        $data[$j] = $gc{$lib}{$j} || 0;
			push(@x_point,$j);
			$max_y = $data[$j] if ($data[$j] > $max_y);
                }
		shift(@data);							# to fix first value not at base position 1
                push(@legends,$legend);
                push(@values,\@data);
                $data_test = $#data ;
                $data_test_max = $data_test if ($data_test > $data_test_max); 
        }
        return("NA") if ($data_test_max == -1);                                 # this is to avoid attempt to plot nothing.
        my $chartobj = new AAGnuplot('title' => "GC Content - $read ",
                                 'xlabel' => "% GC",
                                 'ylabel' => "Fraction of Total Reads",
                                 'xrange' => "[1:100]",
                                 'yrange' => "[0:$max_y]",
                                 'output' => "$output.eps",
                                 'graph' => "lines");
        $chartobj->addOpts('size' => "1.4,1.4");
        $chartobj->plotAoA('data' => \@values,
                                'xlabels' => \@x_point,
                                'legend labels' => \@legends);
        return($chartobj);
}

sub plotComplexityHisto {
	my ($complexity_hash_ref) = ($_[0]);
	my (@x_point, @legends, @values, @data, $k, $legend);
	my $max_y = 0;								# max y = largest # of reads at a given %GC
	my ($data_test, $data_test_max) = (-1, -1);                             # if there is no data to plot, gnuplot hangs
	my %complexity = %$complexity_hash_ref;
	my $output = "$qc_chart_data_subdir/library_complexity";
	print "Plotting library complexity histogram...\n\n";
        foreach my $lib (keys %input) {
		next if ($input{$lib}{'paired'} eq "false");
		$plot_complexity = 1;
                $legend = $lib;
                $legend =~ s/\s+/\_/g;
                my @data;       
                for (my $j=0; $j <= ($max_complexity_count+1); $j++) {
                        $data[$j] = $complexity{$lib}{$j} || 0;
			push(@x_point,$j);
			$max_y = $data[$j] if ($data[$j] > $max_y);
                }
		#shift(@data);							# to fix first value not at base position 1
                push(@legends,$legend);
                push(@values,\@data);
                $data_test = $#data ;
                $data_test_max = $data_test if ($data_test > $data_test_max); 
        }
        return("NA") if ($data_test_max == -1);                                 # this is to avoid attempt to plot nothing.
        my $chartobj = new AAGnuplot('title' => "Library Complexity Duplication Group Count",
                                 'xlabel' => "Duplication Count",
                                 'ylabel' => "Number of Pairs",
                                 'xrange' => "[1:$max_complexity_count]",
                                 'yrange' => "[0:$max_y]",
                                 'output' => "$output.eps",
                                 'graph' => "lines");
        $chartobj->addOpts('size' => "1.4,1.4");
        $chartobj->plotAoA('data' => \@values,
                                'xlabels' => \@x_point,
                                'legend labels' => \@legends);
        return($chartobj);
}

sub rRnaAlign {
	my $ref = $_[0];
	my ($out, $out_met, $total_pairs, $pct, $fq1, $fq2, $index, $r_index, $local_ref, $bam, $fai, $hit_count);
	my $s2fq_cmd;
	$local_ref = $ref;
	$local_ref =~ s/.*\///;					# will make copy in cwd then create index files locally (small file)
	$index   = $ref . ".bwt";
	$r_index = $ref . ".rbwt";
	$fai   = $ref . ".fai";
	foreach my $id (keys %input) {
		$total_pairs = $input{$id}{'r1_total_reads'};
		$fq1 = $input{$id}{'fq1'} || " - ";
		$fq2 = $input{$id}{'fq2'} || " - ";
		#print "TEST FQ $fq1 $fq2\n";
		$out = "read_qc/read_data/$id.rRNA";
		$out_met = "read_qc/metrics_data/$id.rRNA_metrics_flagstat.txt";
		if (! -e $fq1 || ! -e $fq2) {
		 	$bam = $input{$id}{'bam'};
			$fq1 = "read_qc/read_data/$id.1.fq";
			$fq2 = "read_qc/read_data/$id.2.fq";
			if ($input{$id}{'paired'} eq "true") {
				$s2fq_cmd = "java -jar $picard/picard.jar SamToFastq I=$bam F=$fq1 F2=$fq2 VALIDATION_STRINGENCY=SILENT";
			} else {
				$s2fq_cmd = "java -jar $picard/picard.jar SamToFastq I=$bam F=$fq1 VALIDATION_STRINGENCY=SILENT";
			}
			print "Running SamToFastq:\n$s2fq_cmd\n\n";
			system "$s2fq_cmd";
		}
		if ($no_overwrite && -e $out_met && ! -z $out_met ) {
			print "The -N no_overwrite option set and rRNA flagstat file exists ... skipping rRNA alignment ... \n\n";
		} else {
			print "\naligning $id to rRNA ref...\n";
			unless (-e $local_ref ) {
				print "\tcopying reference...\n";
				`cp -s $ref $local_ref` if (! -e $local_ref);
			}
			unless (-e $index && -e $r_index && ! -z $index && ! -z $r_index) {
				print "Building bwa index for $local_ref...\n";
				`$bwa index $local_ref`;
			}
			unless (-e $fai && ! -z $fai) {
				print "building samtools index for $local_ref...\n";
				`samtools faidx $local_ref`
			}
			print "\trunning bwa ...\n";

			print "$bwa aln $local_ref $fq1 > $out.1.sai\n";
			`$bwa aln $local_ref $fq1 > $out.1.sai`; 

			if ($input{$id}{'paired'} eq "true") {
				print "$bwa aln $local_ref $fq2 > $out.2.sai\n";
				`$bwa aln $local_ref $fq2 > $out.2.sai`; 
			}

			if ($input{$id}{'paired'} eq "true") {
				print "$bwa sampe $local_ref $out.1.sai $out.2.sai $fq1 $fq2 > $out.bwa.sam\n";
				`$bwa sampe $local_ref $out.1.sai $out.2.sai $fq1 $fq2 > $out.bwa.sam`;
			} else {
				print "$bwa samse $local_ref $out.1.sai $fq1 -f $out.bwa.sam\n";
				`$bwa samse $local_ref $out.1.sai $fq1 -f $out.bwa.sam`;

			}

			print "\tconverting to bam format...\n";
			print "$samtools view -bt $fai -o $out.bam $out.bwa.sam";
			`$samtools view -bt $fai -o $out.bam $out.bwa.sam`;

			print "Running flagstat...\n";
			print "$samtools flagstat $out.bam > $out_met\n";
			`$samtools flagstat $out.bam > $out_met`;	# previously used samtools view -F12 ; switch to flagstat
		}	
		$hit_count = `grep "with itself and mate mapped" $out_met`;
		$hit_count =~ s/^(\d+)\s+.*/$1/;
		$hit_count = $hit_count/2;
		$pct = sprintf("%.1f",(100*($hit_count/$total_pairs)));
		$input{$id}{'pct_rrna'} = $pct;
		foreach my $tmp_file ("$out.1.sai", "$out.2.sai", "$out.bwa.sam", "$out.bam") {
			`rm $tmp_file` if (-e $tmp_file);
		}
	}

}


sub getReadLnHisto {
	print "Getting read length histogram info from fastqc output ... \n\n";
	foreach my $id (keys %input) {
		foreach my $fastqc_dir ($input{$id}{'r1_fastqc_dir'}, $input{$id}{'r2_fastqc_dir'}) {
			my ($read, $metrics);
			$metrics = $fastqc_dir . "/fastqc_data.txt";
			if ($metrics =~ /\.r2\_fastqc/ || ($metrics =~ /\/na\//)) {
				$read = 2;
			} else {					
				$read = 1;
			}
			next if ($read == 2 && $input{$id}{'paired'} ne "true" );
			my $start_count = 0;
			open(IN,$metrics) || die "Cannot open $metrics metrics file:$!\n\n";
			while(<IN>) {
				my $line = $_;
				chomp($line);
				if ($line =~ /^\>\>Sequence\s+Length\s+Distribution/) {
					$start_count = 1;
					next;
				}
				if ($start_count != 1) {
					next;
				} elsif ($start_count == 1 && $line =~ /END_MODULE/) {
					last;
				} else {
					if ($line =~ /^(\d+)\-\d+\s+(\S+)/) {	
						my ($bin, $num) = ($1, $2);
						$num =~ s/\..*//;
						$rd_ln_histo_r1{$id}{$bin} = $num if ($read == 1);
						$rd_ln_histo_r2{$id}{$bin} = $num if ($read == 2);
						$read_ln_histo_bins{$bin} = 1 unless (exists $read_ln_histo_bins{$bin});
					}
				}	
			}
			close IN;
		}
	}
}

sub plotReadLnHisto {
	my ($read, $rd_ln_hash_ref) = ($_[0], $_[1]);
	my (@x_point, @legends, @values, @data, $k, $legend);
	my $max_y = 0;								# max y = largest # of reads at a given %GC
	my $max_x = 0;
	my ($data_test, $data_test_max) = (-1, -1);                             # if there is no data to plot, gnuplot hangs
	my %rd_ln = %$rd_ln_hash_ref;
	my $output = "$qc_chart_data_subdir/read_length_histogram_" . $read;
	print "Plotting $read read length histogram...\n\n";
        foreach my $lib (keys %input) {
		next if ($read eq "read2" && $input{$lib}{'paired'} eq "false");
                $legend = $lib;
                $legend =~ s/\s+/\_/g;
                my @data;       
		foreach my $bin (sort {$a <=>  $b} keys %read_ln_histo_bins) {
			#print "TEST BIN $bin\t$rd_ln{$lib}{$bin} \n";
                        my $val = $rd_ln{$lib}{$bin} || 0;
                        push(@data,$val);
			push(@x_point,$bin);
			#$max_x = $bin;
			$max_y = $rd_ln{$lib}{$bin} if ($rd_ln{$lib}{$bin} > $max_y);
                }
		#shift(@data);							# to fix first value not at base position 1
                push(@legends,$legend);
                push(@values,\@data);
                $data_test = $#data ;
                $data_test_max = $data_test if ($data_test > $data_test_max); 
        }
        return("NA") if ($data_test_max == -1);                                 # this is to avoid attempt to plot nothing.
        my $chartobj = new AAGnuplot('title' => "Read Length Histogram - $read ",
                                 'xlabel' => "Read Length (bases)",
                                 'ylabel' => "Fraction of Total Reads",
                                 #'xrange' => "[0:$max_x]",
                                 'yrange' => "[0:$max_y]",
                                 'output' => "$output.eps",
                                 'graph' => "lines");
        $chartobj->addOpts('size' => "1.4,1.4");
        $chartobj->plotAoA('data' => \@values,
                                'xlabels' => \@x_point,
                                'legend labels' => \@legends);
        return($chartobj);
}

sub prepareSummaryTable {
	my $id_count = 0;
	my ($short_id, $label);
	open(OUT,">$stats_txt") || die "Cannot write $stats_txt stats output file:  $!\n\n"; 
	print OUT "id\t";
	foreach my $stat (@info_type) {
		print OUT "$stat\t";
	}
	print OUT "\n";
	foreach my $id (keys %input) {
		print OUT "$id\t";
		my @push;
		my @push_toc;
		$id_count++;
		$short_id = $input{$id}{'short_id'} = substr($id,0,10);
		$label    = $input{$id}{'label'};
		push(@push,$id_count);
		foreach my $toc_val ($id_count, $short_id, $label, $id) {
			push(@push_toc,$toc_val);
		}
		foreach my $info_type (@info_type) {
			#print "TYPE $info_type\n";
			my $stat =  $input{$id}{$info_type} || " - ";
			#print "STAT $stat orig $input{$id}{$info_type}\n";
			$stat = commify($stat) if ($stat =~ /\d+/);
			push(@push,$stat); 
			print OUT "$stat\t";
		}	
		print OUT "\n";
		push(@pdf_data, [@push]);
		push(@toc_data,[@push_toc]);
	}
	print OUT "\n";
	close OUT;
		
}

sub createPdf {
	my $out_title = "Read QC";
	if ( defined $pdf_file ) {
		my $tex_file = $pdf_file;
		unless ($tex_file =~ s/\.[Pp][Dd][Ff]$/.tex/) {
			$tex_file .= '.tex';
		}

                open my $TEX_FH, ">", $tex_file or die "Error opening $tex_file for write: $!";
                print $TEX_FH tex_header(-2.29);
                print $TEX_FH "\\section*{$out_title}\n";

                print $TEX_FH "\\section{Read Input Legend}\n";                  		# Print out Stats Table
                print $TEX_FH tex_table("Legend", \@toc_header, \@toc_data, undef, undef, 1.0, undef, 0, 0);   
		print $TEX_FH "\n" . '\clearpage' . "\n";

                print $TEX_FH "\\section{Summary Metrics}\n";                          		# Print out Stats Table
                for (my $i = 0; $i < scalar @pdf_data; $i += 4) {                               # 4 read groups per page
                        my $top_slice = ($#pdf_data <= $i+3 ?  $#pdf_data :  $i +3) ;
                        my @sub_data = @pdf_data[$i ..$top_slice];
                        print $TEX_FH tex_table("Summary Metrics", \@pdf_info_type, \@sub_data, undef, undef, 1.0, undef, 1, 0);   
                         #if ($i % 36 == 0) { print $TEX_FH "\n" . '\clearpage' . "\n";}
                        if ($i % 36 == 0 && $i > 4) {
                                #print $TEX_FH "\n" . '\clearpage' . "\n";
                                print $TEX_FH "\n" . '\clearpage';
                        } elsif ($i % 36 == 0 && $i < 5) {      
                                print $TEX_FH "\n" ;
                        }
                                
                }
		print $TEX_FH "\n" . '\clearpage' . "\n";

		#################################      Print Read Blast Summary Table     ##########################################

		unless ($no_blast) {
			print $TEX_FH "\\section{Read Blast Results}\n";
			&printLongTable(\@read_blast_headers,\@read_blast_data,$TEX_FH);
		}

		#################################        Print Read Length Histogram      ##########################################
		if ($run_read_ln_histo && ! $no_plots) {
			print $TEX_FH "\\section{Read Length Histogram}\n";
			print $TEX_FH $rd_ln_histo_r1_chartobj->writeToLatex() if ($rd_ln_histo_r1_chartobj);
			print $TEX_FH $rd_ln_histo_r2_chartobj->writeToLatex() if ($rd_ln_histo_r2_chartobj);
			print $TEX_FH "\n" . '\clearpage' . "\n";
		}
		#################################      Print Base Quality and GC Plots    ##########################################
		if (! $no_plots) {
			print $TEX_FH "\\section{Mean Base Quality Along Read Position}\n";
			print $TEX_FH $qual_by_cycle_r1_chartobj->writeToLatex() if ($qual_by_cycle_r1_chartobj);
			print $TEX_FH $qual_by_cycle_r2_chartobj->writeToLatex() if ($qual_by_cycle_r2_chartobj);
			print $TEX_FH "\n" . '\clearpage' . "\n";

			print $TEX_FH "\\section{Base Quality Histogram}\n";
			print $TEX_FH $qual_histo_r1_chartobj->writeToLatex()    if ($qual_histo_r1_chartobj);
			print $TEX_FH $qual_histo_r2_chartobj->writeToLatex()    if ($qual_histo_r2_chartobj ne "NA");

			print $TEX_FH "\\section{Read GC Histogram}\n";
			print $TEX_FH $mean_gc_r1_chartobj->writeToLatex()       if ($mean_gc_r1_chartobj);
			print $TEX_FH $mean_gc_r2_chartobj->writeToLatex()       if ($mean_gc_r2_chartobj);
		}
		#######################################  Print Complexity Plots ####################################################
		if (! $no_complexity && ! $no_plots) {
			if ($plot_complexity == 1) {
				print $TEX_FH "\\section{Library Complexity:  Duplication}\n";
				print $TEX_FH $complexity_chartobj->writeToLatex()   if ($complexity_chartobj && $complexity_chartobj ne "NA");
				print $TEX_FH "\n" . '\clearpage' . "\n";
			}
			#if ($plot_complexity_by_start == 1) {
			#	print $TEX_FH "\\section{Library Complexity:  Rarefaction}\n";
			#	foreach my $c_obj (@complexity_by_start_chartobj) {
			#		print $TEX_FH $c_obj->writeToLatex()   if (@complexity_by_start_chartobj ne "NA");
			#		print $TEX_FH "\n" . '\clearpage' . "\n";
			#	}			
			#}	
		}
		#################################     PrintCvg and Insert Size Related Plots  ######################################
		if (! $no_coverage && ! $no_plots) {
			print $TEX_FH "\\section{Coverage Along Reference}\n";
			foreach my $c_obj (@coverage_chartobj) {
				print $TEX_FH $c_obj->writeToLatex();
			}
			print $TEX_FH "\n" . '\clearpage' . "\n";
			print $TEX_FH $coverage_histo_chartobj->writeToLatex()   if ($coverage_histo_chartobj ne "NA");
			print $TEX_FH "\n" . '\clearpage' . "\n";

			print $TEX_FH "\\section{GC Bias in Reference Coverage}\n";
			 foreach my $c_obj (@gc_bias_chartobj) {
			 	print $TEX_FH $c_obj->writeToLatex();
			}
			print $TEX_FH "\n" . '\clearpage' . "\n";

			if ($run_lorenz) {
				print $TEX_FH "\\section{Lorenz Curve Coverage Uniformity}\n";
				foreach my $ref (keys %references) {
					my $png = $references{$ref}{'lorenz_plot'};
					my $title = "Reference $ref";
					if (-e $png) {
						print $TEX_FH "\\begin{figure}[bp!]\n";
						print $TEX_FH "\t\\centering\n";
						print $TEX_FH "\t\\caption{$title}\n";
						print $TEX_FH "\t\\includegraphics[scale=0.75]{$png}\n";
						print $TEX_FH "\\end{figure}\n";
					}
				}	
				print $TEX_FH "\n" . '\clearpage' . "\n";
			}
		}
		if ($plot_insert_size == 1 && ! $no_plots) {
			print $TEX_FH "\\section{Insert Size Histogram}\n";
			print $TEX_FH $insert_size_histo_fr_chartobj->writeToLatex()   if ($insert_size_histo_fr_chartobj ne "NA");
			print $TEX_FH $insert_size_histo_rf_chartobj->writeToLatex()   if ($insert_size_histo_rf_chartobj ne "NA");
			print $TEX_FH "\n" . '\clearpage' . "\n";
		}	
		#####################################        Print FASTQC Plots     ################################################
		if (! $no_plots) {
			print $TEX_FH "\\section{Fastqc Base Quality Along Read Position}\n";
			foreach my $id (keys %input) {
				my $r1_title   = $id . " read 1 Per Base Quality";
				my $r2_title   = $id . " read 2 Per Base Quality";
				my $r1_fq_file = $input{$id}{'r1_fastqc_dir'} . "/Images/per_base_quality.png";
				my $r2_fq_file = $input{$id}{'r2_fastqc_dir'} . "/Images/per_base_quality.png";
				if (-e $r1_fq_file) {
					print $TEX_FH "\\begin{figure}[bp!]\n";
					print $TEX_FH "\t\\centering\n";
					print $TEX_FH "\t\\caption{$r1_title}\n";
					print $TEX_FH "\t\\includegraphics[scale=0.5]{$r1_fq_file}\n" if (-e $r1_fq_file);
					print $TEX_FH "\\end{figure}\n";
				}
				if (-e $r2_fq_file) {
					print $TEX_FH "\\begin{figure}[bp!]\n";
					print $TEX_FH "\t\\centering\n";
					print $TEX_FH "\t\\caption{$r2_title}\n";
					print $TEX_FH "\t\\includegraphics[scale=0.5]{$r2_fq_file}\n" if (-e $r2_fq_file);
					print $TEX_FH "\\end{figure}\n";
				}
			}	
		}
		print $TEX_FH "\n" . '\clearpage' . "\n";

		print $TEX_FH  tex_footer(-2.29);
		close $TEX_FH;

		if (tex2pdf($tex_file)) { 
		        print "PDF created\n";
		} else {
		  	print "PDF creation failed\n";
		}
	}
	my $pdf_final = $pdf_file;
	$pdf_final = $pdf_final . ".pdf" unless ($pdf_final =~ /.pdf$/);
	#system "cp $pdf_final ../";
}




sub printLongTable {
    my ($headers,$data, $TEX_FH) = @_;
    print $TEX_FH "\n\\begin{center}\n";
    my $string = "l" x scalar(@{$headers});
    print $TEX_FH "\\begin{longtable}{".$string."}\n";
    print $TEX_FH "\\caption[Read Blast Results]{Read Blast Results ($subset pairs or unpaired reads)} \\label{grid_mlmmh}\n";
    
    for (my $j = 0; $j <= 1; $j++) {
        
        $string = "";
        if ($j == 1) {
            $string .= "\\multicolumn{".($#{$headers} + 1)."}{c}{{\\tablename} \\thetable{} -- Continued}\\\\[0.5ex]\n";
        }

        $string .= "\\hline\\hline \\\\[-2ex]\n";
        for (my $i = 0; $i <= $#{$headers};$i++) {
            $string .= "   \\multicolumn{1}{c}{\\textbf{".${$headers}[$i]."}} &";
            if ($i == $#{$headers}) {
                $string .= "\\\\[0.5ex] \\hline\n\\\\[-1.8ex]\n\n";
            } else {
                $string .= "\n";
            }
        }
        $string .= "\\endfirsthead\n" if ($j == 0);
        $string .= "\\endhead\n" if ($j == 1);
        print $TEX_FH "$string\n";
    }

    print $TEX_FH "  \\\\[-1.8ex] \\hline\\hline\n";
    print $TEX_FH "   \\multicolumn{".($#{$headers} + 1)."}{r}{{Continued on Next Page\\ldots}} \\\\\n\\endfoot\n";
    print $TEX_FH "  \\\\[-1.8ex] \\hline\\hline\n\\endlastfoot\n";
    
    my $prev_scaff = "scaffold_0";
    for (my $i = 0; $i <= $#{$data}; $i++) {
        #$string = join(" & ",@{$data}[$i]);
        $string = "";
        if ($i % 2 == 0) {
                $string .= "\\rowcolor[gray]{0.8}[1\\tabcolsep] ";
        }
        my $scaf = ${${$data}[$i]}[1];
        if ($scaf ne $prev_scaff ) {
                print $TEX_FH "\\hline\n";
                print $TEX_FH "\\hline\n";
                print $TEX_FH "\\hline\n";
        }
        $string .= join(" & ",@{${$data}[$i]});
        $string =~ s/\%/\\%/g;
        $string =~ s/\_/\\_/g;

        print $TEX_FH $string." \\\\\n";
        $prev_scaff = $scaf;
    }
    print $TEX_FH "\\end{longtable}\\end{center}\n";

}
