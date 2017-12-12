#!/usr/bin/env perl 
use strict;

use Getopt::Std;
use vars qw/ %opt /;
use Tie::IxHash;
#use FindBin;
#FindBin->again;
#use lib "$FindBin::Bin/../lib/tools";
use lib "/gsap/assembly_analysis/tools/lib/";
use AAGnuplot;						# AAGnuplot, LaTeX_aa, scalar_manip_aa to generate plots 
use LaTeX_aa;
use scalar_manip_aa;
my $opt_string = 'hs:l:t:i:NO';
getopts( "$opt_string", \%opt) or die "Input wrong: $!\n";

sub usage() {
    print <<HELP;

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


HELP
    exit;
}

if ($opt{h}) { usage(); exit; }

my ($input_table, $labset_dir, $parent_ssf, $labset_ssf, $data_type);

if ($opt{t}) {
        $data_type  = $opt{t};
	if ($data_type !~ /^pcr/i && $data_type !~ /^jump/i && $data_type !~ /^low/i) {
		print "For -t select either \"fragment\" or \"jump\"\n";
        	usage() and exit;
	}
} else {
        print "User must supply -t <type> data type \"pcr-free\" or \"low-input\" or \"jump\"\n\n";
        usage() and exit;
}

if ($opt{i}) {
        $input_table  = $opt{i};
} else {
        print "User must supply -i <input_table>\n\n";
        usage() and exit;
}

if ($opt{s}) {
        $parent_ssf  = $opt{s};
	$parent_ssf =~ tr/a-z/A-Z/;
} else {
        print "User must supply -s <ssf> parent SSF ticket number\n\n";
        usage() and exit;
}

if ($opt{l}) {
        $labset_ssf  = $opt{l};
	$labset_ssf =~ tr/a-z/A-Z/;
} else {
        print "User must supply -l <labset_ssf> labset SSF ticket number\n\n";
        usage() and exit;
}

##############################################      Variables     ############################################

my $bwa      = "/seq/software/picard/current/3rd_party/bwa/bwa";
my $samtools = "/broad/software/groups/gtba/software/samtools_0.1.18/bin/samtools";
my $picard   = "/seq/software/picard-public/current";
my $gaemr    = "/gsap/assembly_analysis/gaemr-packages/GAEMR/bin";
my $qc_tool =  "/opt/src/read_qc_wrapper.pl";
my %info;
my $data_dir;
my $alignment_dir;
my $qc_dir;

##############################################        Main        ############################################

($labset_dir, $data_dir,$alignment_dir,$qc_dir) = &createDirs($parent_ssf,$labset_ssf,$data_type);

&parseTable($input_table);

&makeBams; 

&alignReads; 

&qcReads;

&parseQc;

&writeOut;


##############################################    Subroutines     ############################################

sub createDirs {
	my ($p, $l, $data) = @_;
	my ($p_dir, $l_dir, $d_dir, $a_dir, $q_dir, $type);
	if ($data =~ /jump/i) {
		$type = "Jumping";
	} elsif ($data =~ /^pcr/i)  {
		$type = "PCR-Free";
	} elsif ($data =~ /^low/i) {
		$type = "LowInputNextera";
	} else {
		print "The -d <data_type> must begin with either low or pcr or jump\n";
		usage() and exit;
	}
	$p_dir = "/btl/projects/SSF/$type/$p";
	$l_dir = "$p_dir/$l";
	$d_dir = "$l_dir/data";
	$a_dir = "$l_dir/alignment";
	$q_dir = "$l_dir/qc";
	#print "$p_dir\n$l_dir\n$d_dir\n$a_dir\n$q_dir\n";
	`mkdir $p_dir` unless ( -d $p_dir);
	`mkdir $l_dir` unless ( -d $l_dir);
	`mkdir $d_dir` unless ( -d $d_dir);
	`mkdir $a_dir` unless ( -d $a_dir);
	`mkdir $q_dir` unless ( -d $q_dir);
	my $local_table = $l_dir . "/$l.qc_input.tbl";
	`cp $input_table $local_table` if (! -e $local_table);
	return($l_dir,$d_dir,$a_dir,$q_dir);
}

sub parseTable {
	my $in = $_[0];
	open(IN,"$in") || die "Cannot open $in table:  $!\n\n";
	while(<IN>) {
		my $line = $_;
		$line =~ s/\*//g;
		chomp($line);
		my @line = split(/\t/,$line);
		next if ($. == 1 && $line =~ /genus/i);
		my $sample_id = $line[0];
		$info{$sample_id}{'genus'}     = $line[1];
		$info{$sample_id}{'species'}   = $line[2];
		$info{$sample_id}{'strain'}    = $line[3];
		$info{$sample_id}{'specimen'}  = $line[4];
		$info{$sample_id}{'gnum'}      = $line[5];
		$info{$sample_id}{'ref_path'}  = $line[6];
		$info{$sample_id}{'fq1'}       = $line[7] || "NA";
		$info{$sample_id}{'fq2'}       = $line[8] || "NA";
		$info{$sample_id}{'seq_bam'}   = $line[9] || "NA";
		if ($info{$sample_id}{'seq_bam'} eq "NA" && ($info{$sample_id}{'fq1'} eq "NA" || $info{$sample_id}{'fq2'} eq "NA")) {
			print "For $sample_id cannot locate either BAM or both R1/R2 fastq file. Check $line[7] $line[8] $line[9]\n\n";
			exit;
		}
	}
	close IN;
	foreach my $i (keys %info) {
		#print "$i\t$info{$i}{'genus'}\t$info{$i}{'fq1'}\t$info{$i}{'fq2'} \n";
	}
}

sub makeBams {
	chdir $data_dir;
	foreach my $sample_id (keys %info) {
		my $r1       = $info{$sample_id}{'fq1'};
		my $r2       = $info{$sample_id}{'fq2'};
		my $seq_bam  = $info{$sample_id}{'seq_bam'};
		my $bam = $data_dir . "/" . "$sample_id." . $data_type . ".unmapped.bam";

		if (-e $seq_bam  && $r1 !~/\//) {  	# if seq_bam exists and no path in fq1 entry
			`ln -s $seq_bam $bam`;		# link over the seq bam
		} else {				# otherwise create unmapped bam from fq files
			my $cmd = "java -jar $picard/picard.jar  FastqToSam F1=$r1 F2=$r2 O=$bam SAMPLE_NAME=sample \n";
			if (! -e $bam) {
				print "No $bam found.  Running:\n$cmd\n\n";
				system $cmd unless $opt{N};	
			} elsif (-z $bam) {
				print "Empty $bam found.  Running:\n$cmd\n\n";
				system $cmd unless $opt{N};	
			} else {
				print "$bam exists ... skipping FastqToSam \n\n";
			}
		}
		$info{$sample_id}{'unmapped_bam'}     = $bam;
	}
}

sub alignReads {
	foreach my $sample_id (keys %info) {
		chdir $alignment_dir;
		my $sample_dir = $alignment_dir . "/" . $sample_id;
		`mkdir $sample_id` if (! -e "$sample_id");
		chdir $sample_id;
		my $in   = $info{$sample_id}{'unmapped_bam'};
		my $ref = $info{$sample_id}{'ref_path'};
		my $out = $sample_dir . "/$sample_id." . $data_type . ".aligned";
		my $cmd = "$gaemr/align_reads.py -r $ref -i $in -o $out\n";
		my $pwd = `pwd`;
		if (! -e "$out.bam") {
			print "Running:\n$cmd\n\n";
			system $cmd  unless $opt{N};	
		} elsif (-z "$out.bam") {
			print "Empty $out.bam found.  Running:\n$cmd\n\n";
			system $cmd unless $opt{N};	
		} else {
			print "$out.bam exists ... skipping align_reads.py \n\n";
		}
		$info{$sample_id}{'aligned_bam'} = "$out.bam";
	}
}

sub qcReads {
	chdir $qc_dir;
	foreach my $sample_id (keys %info) {
		my $bam = $info{$sample_id}{'aligned_bam'};
		my $cmd = "$qc_tool -T -N -B -o $sample_id -b $bam"; 
		my $out = "$qc_dir/$sample_id.metrics.txt";
		if (! -e $out) {
			print "Running QC tool:\n$cmd\n\n";
			system $cmd  unless $opt{N};
		} elsif (-z $out) {
			print "Empty $out found.  Running:\n$cmd\n\n";
			system $cmd unless $opt{N};	
		} else {
			print "$out exists ... skipping read_qc for sample $sample_id\n\n";
		}
		$info{$sample_id}{'qc_out'} = "$out";
	}
}

sub parseQc {
	chdir $labset_dir;
	my $out = "$labset_dir/$labset_ssf.read_qc_summary.tbl";
	open(OUT,">$out") || die "Cannot write to $out file:  $!\n\n";
	foreach my $sample_id (keys %info) {
		my $in = $info{$sample_id}{'qc_out'};
		open(IN,"$in") || die "Cannot open the $in qc file:  $!\n\n";
		while(<IN>) {
			my $line = $_;
			chomp($line);
			my @line = split(/\t/,$line);
			next if ($line !~ /\S+/ || $line =~ /r1_mean_qual/);
			$info{$sample_id}{'num_pairs'}        = $line[3];
			$info{$sample_id}{'read_ln'}          = $line[7];
			$info{$sample_id}{'r1_mean_qual'}     = $line[9];
			$info{$sample_id}{'r2_mean_qual'}     = $line[10];
			$info{$sample_id}{'r1_mean_gc'}       = $line[19];
			$info{$sample_id}{'r2_mean_gc'}       = $line[20];
			$info{$sample_id}{'pct_duplication'}  = $line[21];
			$info{$sample_id}{'est_lib_size'}     = $line[22];
			$info{$sample_id}{'pct_adapter'}      = sprintf("%.3f",$line[24]);
			$info{$sample_id}{'pct_aligned'}      = $line[29];
			$info{$sample_id}{'fr_mean_insert'}   = $line[32];
			$info{$sample_id}{'fr_median_insert'} = $line[33];
			$info{$sample_id}{'rf_mean_insert'}   = $line[34];
			$info{$sample_id}{'rf_median_insert'} = $line[35];
			$info{$sample_id}{'mean_cvg'}         = $line[36];
			$info{$sample_id}{'pct_jump'}         = $line[38];
		}
		close IN;
	}
	print OUT "metric\tmean\tmin\tmax\n";
	foreach my $i ("num_pairs","read_ln","r1_mean_qual","r2_mean_qual","r1_mean_gc","r2_mean_gc","pct_duplication","est_lib_size","pct_adapter","pct_aligned","fr_mean_insert","fr_median_insert","rf_mean_insert","rf_median_insert","mean_cvg","pct_jump") {
		next if ($data_type =~ /pcr/ && ($i eq "rf_mean_insert" || $i eq "rf_median_insert"));
		my ($count,$sum,$mean,$min,$max) = (0,0,0,100000000000000000000000000000000000000000,0);
		foreach my $sample (keys %info) {
			my $val = $info{$sample}{$i};
			$val =~ s/\,//g;
			if ($val =~ /\d+/) {
				$count++;
				$sum += $val;
				$val  = sprintf("%.2f",$val) if ($i eq "pct_adapter");
				$min  = $val if ($val < $min);
				$max  = $val if ($val > $max);
			}
		}

		$mean = sprintf("%.1f",($sum / $count));
		$mean = sprintf("%.2f",($sum / $count)) if ($i eq "pct_adapter");
		$mean =~ s/\.\d+// if ($i eq "num_pairs" || $i eq "read_ln" || $i eq "est_lib_size" || $i eq "fr_mean_insert" || $i eq "fr_median_insert" || $i eq "rf_mean_insert" || $i eq "rf_median_insert" || $i eq "mean_cvg"); 
		$mean = commify($mean) if ($i eq "num_pairs" || $i eq "est_lib_size" || $i eq "rf_mean_insert" || $i eq "rf_median_insert");
		$min  = commify($min) if ($i eq "num_pairs" || $i eq "est_lib_size" || $i eq "rf_mean_insert" || $i eq "rf_median_insert");
		$max = commify($max) if ($i eq "num_pairs" || $i eq "est_lib_size" || $i eq "rf_mean_insert" || $i eq "rf_median_insert");
		print OUT "$i\t$mean\t$min\t$max\n";
	}
	close OUT;
}

sub writeOut {
	chdir $labset_dir;
	my $out = "$labset_dir/$labset_ssf.read_qc_metrics.tbl";
	open(OUT,">$out") || die "Cannot write to $out file:  $!\n\n";
	print OUT "id\tgenus\tgnumber\tspecies\tstrain\ttotal_pairs\tr1_mean_rd_ln\tr1_mean_qual\tr2_mean_qual\t";
	print OUT "r1_mean_gc\tr2_mean_gc\t";
	print OUT "pct_duplication\testimate_library_size\tpct_adapter\tpct_aligned\tfr_mean_insert\t";
	print OUT "fr_median_insert\trf_mean_insert\trf_median_insert\tmean_cvg\tpct_jump\n";
	foreach my $sample_id (keys %info) {
		print OUT "$sample_id\t$info{$sample_id}{'gnum'}\t$info{$sample_id}{'genus'}\t$info{$sample_id}{'species'}\t";
		print OUT "$info{$sample_id}{'strain'}\t";
		print OUT "$info{$sample_id}{'num_pairs'}\t$info{$sample_id}{'read_ln'}\t";
		print OUT "$info{$sample_id}{'r1_mean_qual'}\t$info{$sample_id}{'r2_mean_qual'}\t";
		print OUT "$info{$sample_id}{'r1_mean_gc'}\t$info{$sample_id}{'r2_mean_gc'}\t";
		print OUT "$info{$sample_id}{'pct_duplication'}\t$info{$sample_id}{'est_lib_size'}\t";
		print OUT "$info{$sample_id}{'pct_adapter'}\t$info{$sample_id}{'pct_aligned'}\t";
		print OUT "$info{$sample_id}{'fr_mean_insert'}\t$info{$sample_id}{'fr_median_insert'}\t";
		print OUT "$info{$sample_id}{'rf_mean_insert'}\t$info{$sample_id}{'rf_median_insert'}\t";
		print OUT "$info{$sample_id}{'mean_cvg'}\t$info{$sample_id}{'pct_jump'}\n";
	}
	close OUT;
}
