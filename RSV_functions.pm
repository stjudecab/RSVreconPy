#!/usr/bin/perl

package RSV_functions;

use strict;
use warnings;
use Exporter;

our @ISA= qw( Exporter );

# these CAN be exported.
our @EXPORT_OK = qw( get_sub_folders elements_not_in_array pct_sum process determine_subtype);

# these are exported by default.
our @EXPORT = qw( get_sub_folders elements_not_in_array pct_sum process determine_subtype );

###################################################
##      functions
###################################################

# get sub folders for a dir
sub get_sub_folders 
{
    my ($folder_path) = @_;
    my @sub_folders;

    opendir(my $dh, $folder_path) or die "Cannot open directory: $!";

    while (my $entry = readdir $dh) {
        next if $entry eq '.' || $entry eq '..';
        my $sub_folder_path = "$folder_path/$entry";
        if (-d $sub_folder_path) {
            push @sub_folders, $entry;
        }
    }

    closedir $dh;
    return @sub_folders;
}

# test if elements in array A in array B, if not, return elements that are not in B
sub elements_not_in_array 
{
    my ($array_a, $array_b) = @_;
    my @not_in_b;

    foreach my $element (@$array_a) {
        push @not_in_b, $element unless grep { $_ eq $element } @$array_b;
    }

    return @not_in_b;
}

# sum the pct values
sub pct_sum 
{
   my $sum = 0;
 
   foreach my $item (@_){
      $item =~ s/\%//;
      $sum += $item;
   }

   return $sum;
}

# assemble sequences from IGV counts
sub process
{
	my $file_name = shift;
	my $sequence = "";


	open(FILE,$file_name) or die ("can not open $file_name\n");
	while(my $line = <FILE>)
	{
		$line =~ s/\n//;
		if($line =~ /^track/)
		{
			next;
		}
		elsif($line =~ /^#/)
		{
			next;
		}
		elsif($line =~ /variableStep\schrom=([^\s]+)\s.+/)
		{
			next;
		}
		elsif($line =~ /^\d+/)
		{
			if($line =~ /^(\d+)\t(\d+)\.0\t(\d+)\.0\t(\d+)\.0\t(\d+)\.0.+/)
			{
				#print STDERR "$1\t";
				my $pos = $1;
				my $a = $2;
				my $c = $3;
				my $g = $4;
				my $t = $5;

				# make consensus sequence
				my $cov = $a + $c + $g + $t;
				if ($a / $cov > 0.5) {
					$sequence .= "A";
				} elsif ($c / $cov > 0.5) {
					$sequence .= "C";
				} elsif ($g / $cov > 0.5) {
					$sequence .= "G";
				} elsif ($t / $cov > 0.5) {
					$sequence .= "T";
				} else {
					if ($a > $c && $a > $g && $a > $t ) {
						$sequence .= "A";
					} elsif ($c > $a && $c > $g && $c > $t ){
						$sequence .= "C";
					}
					 elsif ($g > $a && $g > $c && $g > $t ){
						$sequence .= "G";
					}
					 elsif ($t > $a && $t > $g && $t > $c ){
						$sequence .= "T";
					}
				}		
			}
			else
			{
				print "error\n";
			}
		}
		else
		{
			print "error1\n";
		}
	}
	close FILE;

	return $sequence;
}

# determine the subtype
sub determine_subtype
{
	my %hash = %{$_[0]};

	foreach my $ref (sort {$hash{$b} <=> $hash{$a}} keys %hash) {
		return $ref;
	}
}

1;
