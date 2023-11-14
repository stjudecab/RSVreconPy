#!/usr/bin/perl

BEGIN {unshift @INC,"/research/groups/cab/projects/automapper/common/lli75/RSV"};

###################################################
##      load packages
###################################################

use strict;
use JSON;
#use Data::Dumper;
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

# create output CSV
my $CSV_header_lv_1 = ',before_filtering, before_filtering, before_filtering, after_filtering, after_filtering, after_filtering, QC rate,';
foreach my $ref_db (@all_ref_db) {
	$CSV_header_lv_1 .= "Mapping $ref_db,Mapping $ref_db,Mapping $ref_db,Mapping $ref_db,";
}
$CSV_header_lv_1 .= "Subtype suggestion";

my $CSV_header_lv_2 = 'Sample name,total_reads, q20_rate, q30_rate, total_reads, q20_rate, q30_rate, QC rate,';
foreach my $ref_db (@all_ref_db) {
	$CSV_header_lv_2 .= 'Uniquely mapped reads %, MULTI-MAPPING READS %, UNMAPPED READS%, CHIMERIC READS%,';
}
$CSV_header_lv_2 .= "Subtype suggestion";

my $report = $working_folder_name."/Report.csv";
open(OUT, '>', $report) or die $!;
print OUT $CSV_header_lv_1."\n";
print OUT $CSV_header_lv_2."\n";

my $sequence_file = $working_folder_name."/Sequence.fasta";
open(FASTA, '>', $sequence_file) or die $!;

# collect information from each sample
my @Sample_folders = <"$working_folder_name/*">;
foreach my $cur_folder (@Sample_folders) {
	if(-f $cur_folder) {
		next;
	}

	print STDERR "Start process $cur_folder ... \n";

	# Step 1: collect QC info from JSON file
	my $QC_json_file = $cur_folder.'/fastp.json';

	my $json_str;
	open(TXT, '<', $QC_json_file) or die "Can not open $QC_json_file!\n";
	while(<TXT>) {
		$json_str .= $_;
	}
	close TXT;

	my $decoded_json = decode_json($json_str);

	my $before = $decoded_json->{'summary'}->{'before_filtering'};
	my $after = $decoded_json->{'summary'}->{'after_filtering'};
	######## fields: total_reads, total_bases, q20_bases, q30_bases, q20_rate, q30_rate, read1_mean_length, read2_mean_length, gc_cont

	# Step 2:  collect mapping info for all ref DB

	my @result_array = ();

	foreach my $ref_db (@all_ref_db) {
		my $cur_map_log_file = $cur_folder."/$ref_db/Log.final.out";

		my %cur_map_hash = {};
		open(FH, '<', $cur_map_log_file) or die "Can not open $cur_map_log_file!\n";
		while(<FH>){
			if ( $_ =~ /\|/ ){
				my $tmp_line = $_;
				$tmp_line =~ s/[\s\%\:]//g;
				my @tmp_arr = split(/\|/,$tmp_line);
				$cur_map_hash{$tmp_arr[0]} = $tmp_arr[1];
			}
		}

		push(@result_array, {%cur_map_hash});
	}

	# Step 3: output info to a CSV
	my $rate = $after->{'total_reads'}/$before->{'total_reads'}*100;
	my $QC_str = "$before->{'total_reads'},$before->{'q20_rate'},$before->{'q30_rate'},$after->{'total_reads'},$after->{'q20_rate'},$after->{'q30_rate'},$rate,";
	
	my $map_str = '';
	my %uniquely_mapped_reads_hash = ();
	for(my $i=0; $i<@result_array; $i++) { 
		#print $i."\t$all_ref_db[$i]\t$result_array[$i]{'Uniquelymappedreads'}\n";
		$map_str .= $result_array[$i]{'Uniquelymappedreads'}.",".pct_sum($result_array[$i]{'ofreadsmappedtomultipleloci'},$result_array[$i]{'ofreadsmappedtotoomanyloci'}).",".pct_sum($result_array[$i]{'ofreadsunmappedtoomanymismatches'},$result_array[$i]{'ofreadsunmappedtooshort'},$result_array[$i]{'ofreadsunmappedother'}).','.$result_array[$i]{'ofchimericreads'}.",";
		$uniquely_mapped_reads_hash{$all_ref_db[$i]} = $result_array[$i]{'Uniquelymappedreads'};
	} 

=pod
	foreach my $key (keys %uniquely_mapped_reads_hash) {
		print "$key => $uniquely_mapped_reads_hash{$key}\n";
	}
=cut

	my $subtype_str = determine_subtype(\%uniquely_mapped_reads_hash);

	my$sample = $cur_folder;
	$sample =~ s/.+\///;

	print OUT $sample.','.$QC_str.$map_str.$subtype_str."\n";


	# Step 4: assemble genome sequence
	my $wig_file = $cur_folder."/$subtype_str/alignments.cov.wig";

	my $genome_sequence = process($wig_file);
	print FASTA ">$sample\n";
	print FASTA $genome_sequence."\n";
}

close OUT;
close FASTA;


###################################################
##      program finished
###################################################

print "Report and sequences have been generated!\n";


