#!perl -w

use warnings;
use strict;

use Getopt::Long;
use Cwd 'abs_path';
use POSIX;

sub usage {
	print <<USAGE;

usage:
	perl $0 [options]
description:

options:
	-help: print help info
	-fq1    *   <str>   : input 1.fq, either .gz or .fq file.
	-fq2    *   <str>   : input 2.fq, either .gz or .fq file
	-out    *   <str>   : output prefix, out_1.fq out_2.fq out.dup.list   out.dup.stat
	-N      *   <float> : reads will be discarded if "N"s exceed the certain percentage of read length, default: 0.05
	-trim   *   <int>   : bases with quality lower than <int> at both side will be trimmed, default: 20
	-tech   *   <int>   : 1, amplicon sequencing; 2, metagenomic sequencing (default); 3, metatranscriptomic sequencing
	-phred  *   <int>   : input quals are: 1, phred+33 (default); 2, phred+64
e.g.:
	perl $0 -fq1 input_1.fq -fq2 input_2.fq -out output_QC -tech 2 -trim 20 -N 0.05 -phred 33 
USAGE
}

my ( $help, $fq1, $fq2, $out, $Ns, $trim, $tech, $phred )
  ;    ###added by dushuo
GetOptions(
	"help"    => \$help,
	"fq1=s"   => \$fq1,
	"fq2=s"   => \$fq2,
	"out=s"   => \$out,
	"N=f"     => \$Ns,
	"trim=i"  => \$trim,
	"tech=i"  => \$tech,
	"phred=i" => \$phred,
);

if ( defined $help || ( !defined $fq1 ) ) {
	&usage();
	exit 0;
}

if ( !defined($out) )     { $out     = "qc"; }
if ( !defined($Ns) )      { $Ns      = 0.05; }
if ( !defined($trim) )    { $trim    = 20; }
if ( !defined($tech) )    { $tech    = 2; }
if ( !defined($phred) )   { $phred   = 33; }

my $start_time = time();
open( LOG, ">$out.log" ) or die "can't create $out.log";

# STEP 1: Read input sequence file(s)
print LOG "STEP 1: Read input sequence file(s)\n";

my ( %fq1, %qual1, %fq2, %qual2 ) = ();
my $fq1_path = abs_path($fq1);
if($fq1_path =~ /.*\.gz/){
	open( FQ1, "gzip -dc $fq1_path |" ) or die "can't open file $fq1_path";
}
else{
	open( FQ1, "$fq1_path" ) or die "can't open file $fq1_path";
}
while (1) {
	my $line1 = <FQ1>;
	my $line2 = <FQ1>;
	my $line3 = <FQ1>;
	my $line4 = <FQ1>;
	unless ( defined($line3) && $line3 ne "" ) { last; }
	my $name = $line1;
	$name =~ s/^@|\s+.*$//ig;
	$line2 =~ s/\s+$//ig;
	$line4 =~ s/\s+$//ig;
	$fq1{$name}   = $line2;
	$qual1{$name} = $line4;
}
close FQ1;
my $num_seqs_1 = scalar( keys %fq1 );

if ( defined($fq2) ) {
	my $fq2_path = abs_path($fq2);
	if($fq2_path =~ /.*\.gz/){
		open( FQ2, "gzip -dc $fq2_path |" ) or die "can't open file $fq2_path";
	}
	else{
		open( FQ2, "$fq2_path" ) or die "can't open file $fq2_path";
	}
	while (1) {
		my $line1 = <FQ2>;
		my $line2 = <FQ2>;
		my $line3 = <FQ2>;
		my $line4 = <FQ2>;
		unless ( defined($line3) && $line3 ne "") { last; }
		my $name = $line1;
		$name =~ s/^@|\s+.*$//ig;
		$line2 =~ s/\s+$//ig;
		$line4 =~ s/\s+$//ig;
		$fq2{$name}   = $line2;
		$qual2{$name} = $line4;
	}
	close FQ1;
	my $num_seqs_2 = scalar( keys %fq2 );
	if ( $num_seqs_2 == $num_seqs_1 ) {
		print LOG
		  "	Inputs are pair-end reads, with $num_seqs_1 sequences detected\n";
	}
	else {
		print LOG
"	Inputs are pair-end reads, but two sides don't have equal number of sequences\n";
		last;
	}
}
else {
	print LOG
	  "	Inputs are single-end reads, with $num_seqs_1 sequences detected\n";
}

my $step1_time = time();
my $step1_duration = $step1_time - $start_time;
my ($hour, $min, $sec) = & time_format($step1_duration);

print LOG "	Running time: $hour day : $min min : $sec sec\n";

# STEP 2: Discard duplicated paired-end reads caused by PCR biases.
print LOG
"STEP 2: Discard duplicated paired-end reads caused by PCR biases (only for metagenomic sequencing, -t 2)\n";
my ($num_seqs_rm_dup_1, $dup_rate) = (0, 0);
if ( $tech == 2 ) {
	my %non_repeat_seqs = ();
	if ( defined($fq2) ) {
		foreach my $key ( sort keys %fq1 ) {
			if ( exists( $non_repeat_seqs{ $fq1{$key} . $fq2{$key} } ) ) {
				delete( $fq1{$key} );
				delete( $fq2{$key} );
				delete( $qual1{$key} );
				delete( $qual2{$key} );
			}
			else {
				$non_repeat_seqs{ $fq1{$key} . $fq2{$key} } = 1;
			}
		}
	}
	else {
		foreach my $key ( sort keys %fq1 ) {
			if ( exists( $non_repeat_seqs{ $fq1{$key} } ) ) {
				delete( $fq1{$key} );
				delete( $qual1{$key} );
			}
			else {
				$non_repeat_seqs{ $fq1{$key} } = 1;
			}
		}
	}
	$num_seqs_rm_dup_1 = scalar( keys %fq1 );
	$dup_rate = 1 - ( $num_seqs_rm_dup_1 / $num_seqs_1 );
	print LOG
"	Inputs are metagenomic data, duplicated reads caused by PCR biases need to be discarded\n";
	print LOG
"	$num_seqs_rm_dup_1 of reads are non-redundant, the duplication rate is $dup_rate\n";
}
else {
	$num_seqs_rm_dup_1 = $num_seqs_1;
	print LOG
	  "	Inputs are not metagenomic data, no need to remove duplicated reads\n";
}

my $step2_time = time();
my $step2_duration = $step2_time-$step1_time;
($hour, $min, $sec) = & time_format($step2_duration);

print LOG "	Running time: $hour day : $min min : $sec sec\n";

# STEP 3: Discard reads with excess number of "N"s
print LOG "STEP 3: Discard reads with excess number of \"N\"s\n";

if ( defined($fq2) ) {
	for my $key ( keys %fq1 ) {
		my $seq = $fq1{$key} . $fq2{$key};
		my $rt  = &rm_excess_Ns($seq);
		if ( $rt == 0 ) {
			delete( $fq1{$key});
			delete( $fq2{$key});
			delete( $qual1{$key});
			delete( $qual2{$key});
		}
	}
}
else {
	for my $key ( keys %fq1 ) {
		my $seq = $fq1{$key};
		my $rt  = &rm_excess_Ns($seq);
		if ( $rt == 1 ) {
			delete( $fq1{$key});
			delete( $qual1{$key});
		}
	}
}
my $num_seqs_rm_dup_N_1 = scalar( keys %fq1 );
my $dup_rate_N = 1 - ( $num_seqs_rm_dup_N_1 / $num_seqs_rm_dup_1 );
print LOG
"	Inputs are metagenomic data, duplicated reads caused by PCR biases need to be discarded\n";
print LOG
"	$num_seqs_rm_dup_N_1 of reads are non-redundant, the duplication rate is $dup_rate_N\n";

my $step3_time = time();
my $step3_duration = $step3_time-$step2_time;
($hour, $min, $sec) = & time_format($step3_duration);

print LOG "	Running time: $hour day : $min min : $sec sec\n";

# STEP 4: Trim bases with low quality at both ends
print LOG "STEP 4: Trim bases with low quality at both ends\n";

if ( defined($fq2) ) {
	for my $key ( keys %fq1 ) {
		my ($index_5_1, $index_3_1, $flag1) = &trim_reads_both_end($qual1{$key});
		my ($index_5_2, $index_3_2, $flag2) = &trim_reads_both_end($qual2{$key});
		if($flag1 == 1 && $flag2 == 1){
			$fq1{$key} = substr( $fq1{$key}, $index_5_1, length($fq1{$key}) - $index_5_1 - $index_3_1 );
			$fq2{$key} = substr( $fq2{$key}, $index_5_2, length($fq2{$key}) - $index_5_2 - $index_3_2 );
			$qual1{$key} = substr( $qual1{$key}, $index_5_1, length($qual1{$key}) - $index_5_1 - $index_3_1 );
			$qual2{$key} = substr( $qual2{$key}, $index_5_2, length($qual2{$key}) - $index_5_2 - $index_3_2 );
		}
		else{
			delete($fq1{$key});
			delete($qual1{$key});
			delete($fq2{$key});
			delete($qual2{$key});
		}
	}
}
else {
	for my $key ( keys %fq1 ) {
		my ($index_5_1, $index_3_1, $flag1) = &trim_reads_both_end($qual1{$key});
		if($flag1 == 1){
			$fq1{$key} = substr( $fq1{$key}, $index_5_1, length($fq1{$key}) - $index_5_1 - $index_3_1 );
			$qual1{$key} = substr( $qual1{$key}, $index_5_1, length($qual1{$key}) - $index_5_1 - $index_3_1 );
		}
		else{
			delete($fq1{$key});
			delete($qual1{$key});
		}
	}
}

my $num_seqs_rm_dup_N_trim_1 = scalar( keys %fq1 );
my $dup_rate_N_trim = 1 - ( $num_seqs_rm_dup_N_trim_1 / $num_seqs_rm_dup_N_1 );

print LOG
"	$num_seqs_rm_dup_N_trim_1 of reads are left, the discarded rate is $dup_rate_N_trim\n";

my $step4_time = time();
my $step4_duration = $step4_time-$step3_time;
($hour, $min, $sec) = & time_format($step4_duration);

print LOG "	Running time: $hour day : $min min : $sec sec\n";

# STEP 5: Output files
print LOG "STEP 5: export and save output files\n";
if ( defined($fq2) ) {
	open (OUT1, ">$out"."_1.fq") or die "can't creat file $out"."_1.fq";
	open (OUT2, ">$out"."_2.fq") or die "can't creat file $out"."_2.fq";
	foreach my $key (sort keys %fq1){
		print OUT1 "@".$key." 1:N:0\n$fq1{$key}\n"."+\n$qual1{$key}\n";
		print OUT2 "@".$key." 2:N:0\n$fq2{$key}\n"."+\n$qual2{$key}\n";
	}
	close OUT1;
	close OUT2;
}
else {
	open (OUT1, ">$out"."_1.fq") or die "can't creat file $out"."_1.fq";
	foreach my $key (sort keys %fq1){
		print OUT1 "@".$key." 1:N:0\n$fq1{$key}\n"."+\n$qual1{$key}\n";
	}
	close OUT1;
}

my $step5_time = time();
my $step5_duration = $step5_time-$step4_time;
($hour, $min, $sec) = & time_format($step5_duration);

print LOG "	Running time: $hour day : $min min : $sec sec\n";

# Done!
print LOG "Quality control has been done\n";
my $total_runningtime = $step5_time - $start_time;
($hour, $min, $sec) = & time_format($total_runningtime);

print LOG "The total running time: $hour day : $min min : $sec sec\n";

#######################################END#####################################

sub rm_excess_Ns {
	my $seq    = $_[0];
	my $length = length($seq);
	$seq =~ s/N|n//ig;
	if ( length($seq) / $length > 1-$Ns ) {
		return 1;
	}
	else {
		return 0;
	}
}

sub trim_reads_both_end {
	my $quality = $_[0];
	my @qualities = split //, $quality;
	my ( $index_5, $index_3, $counter ) = (0,0,0);
	for ( my $i = 0 ; $i < scalar(@qualities) ; $i++ ) {
		if ( ord( $qualities[$i] ) - $phred < $trim ) {
			$index_5++;
			$counter++;
		}
		else {
			last;
		}
	}
	for ( my $i = scalar(@qualities) - 1 ; $i >= 0 ; $i-- ) {
		if ( ord( $qualities[$i] ) - $phred < $trim ) {
			$index_3++;
			$counter++;
		}
		else {
			last;
		}
	}
	if (   scalar(@qualities) - $index_3 - $index_5 > 0 ){
		return ($index_5, $index_3, 1);
	}
	else{
		return ($index_5, $index_3, 0);
	}
}

sub time_format{
	my $duration = $_[0];
	my $hour = int($duration/3600);
	my $min = int(($duration%3600)/60);
	my $sec = $duration%60;
	return ($hour, $min, $sec);
}
