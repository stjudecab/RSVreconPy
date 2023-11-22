#!/usr/bin/perl

BEGIN {unshift @INC,"/research/groups/cab/projects/automapper/common/lli75/RSV"};

###################################################
##      load packages
###################################################

use strict;
use YAML::XS 'LoadFile';

# load RSV functions
use RSV_functions;


###################################################
##      load parameters from YAML file
###################################################
my $config_file = shift;
my $config = LoadFile($config_file) or die "Can not open $config_file as a YAML file!";

my $data_folder_name = $config->{'DATA_DIR'} or die "Please set up a DATA_DIR!";
my $reference_folder_name = $config->{'REFERENCE_DIR'} or die "Please set up a REFERENCE_DIR!";
my $working_folder_name = $config->{'OUTPUT_DIR'} or die "Please set up a OUTPUT_DIR!";

my $star_ThreadN = $config->{'STAR_THREAD_N'} // 16;
my $ref_use = $config->{'REF_USE'} // 'all';


###################################################
##      input data read
###################################################
$data_folder_name =~ s/\/\//\//g;

my %sample_hash;

opendir(my $dh, $data_folder_name) or die "Can not open $data_folder_name: $!";
while (readdir $dh) 
{
	if($_ =~ /\.fastq/)
	{
		my $file = $_;
		my $id = $file;
		$id =~ s/_R\d.+//;
		if($file =~ /_R1/)
		{
			$sample_hash{$id}{1} = $file;
		}
		elsif($file =~ /_R2/)
		{
			$sample_hash{$id}{2} = $file;
		}
	}
}
closedir $dh;


###################################################
##      setup working dir
###################################################
$working_folder_name =~ s/\/\//\//g;
if(!-e $working_folder_name)
{
	mkdir($working_folder_name) or die "Can not create the $working_folder_name!";
}


###################################################
##      start process
###################################################

foreach (sort {$a <=> $b} keys %sample_hash)
{
	print STDERR "Start process $_ ... \n";

	# skip Undetermined reads
	if ($_ =~ /Undetermined/) 
	{
		next;
	}

	my $cur_folder = $working_folder_name."/".$_;
	if(!-e $cur_folder)
	{
		mkdir($cur_folder) or die "Can not create the $cur_folder! Check the path and premission!";
	}

	my $original_read1 = $data_folder_name."/".$sample_hash{$_}{1};
	my $original_read2 = $data_folder_name."/".$sample_hash{$_}{2};

	
	##############################################################################
	### step 1    QC using fastp
	##############################################################################
	print STDERR "######################\tStep 1 Quality check and trimming using fastp ... \t######################\n";

	my $read1_trim = "reads_trimed_R1.fastq";
	my $read2_trim = "reads_trimed_R2.fastq";

	my $cmd = "cd $cur_folder; fastp -i $original_read1 -I $original_read2 -o $read1_trim -O $read2_trim";
	# print STDERR $cmd;
	`$cmd`;

	print STDERR "######################\nFinish Quality trimming . \n######################\n";

	
	##############################################################################
	## step 2    STAR mapping
	##############################################################################
	print STDERR "######################\tStep 2 STAR mapping ... \t######################\n";

	# get all reference DB
	my @all_ref_db = get_sub_folders($reference_folder_name);

	if ($ref_use ne "all") {
		my @ref_use_list = split(',',$ref_use);
		my @not_exist = elements_not_in_array(\@ref_use_list, \@all_ref_db);
		if (@not_exist) {
			my $string = join ',', @not_exist;
			die "These references can not be found under your reference folder: $string\n";
		} else {
			@all_ref_db = @ref_use_list;
		}
	}

	foreach my $ref_db (@all_ref_db) {
		my $refdb_A = "$reference_folder_name/$ref_db";

		my $map_A_folder = $cur_folder.'/$ref_db/';

		if(!-e $map_A_folder)
		{
			mkdir($map_A_folder) or die "Can not create the $map_A_folder! Check the path and premission!";
		}

		my $cmd = "cd $cur_folder; STAR --runThreadN $star_ThreadN --genomeDir $refdb_A --outSAMtype BAM SortedByCoordinate --readFilesIn $read1_trim $read2_trim --outFileNamePrefix $map_A_folder";
		# print STDERR $cmd;
		`$cmd`;
	}

	print STDERR "######################\nFinish STAR mapping . \n######################\n";

	
	##############################################################################
	## step 3    index and Count using Samtools and IGVtools 
	##############################################################################
	print STDERR "######################\tStep 3 Index and Count ... \t######################\n";

	
	foreach my $ref_db (@all_ref_db) {
		my $map_A_folder = $cur_folder.'/$ref_db/';
		my $refdb_A = "$reference_folder_name/$ref_db";

		my $cmd = "cd $map_A_folder; samtools index Aligned.sortedByCoord.out.bam";
		# print STDERR $cmd;
		`$cmd`;

		my $cmd = "cd $map_A_folder; igvtools count -z 5 -w 1 --bases  Aligned.sortedByCoord.out.bam  alignments.cov.wig  $refdb_A.fasta";
		# print STDERR $cmd;
		`$cmd`;
	}

	print STDERR "######################\nFinish Index and Count \n######################\n";

}

