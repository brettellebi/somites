#!/usr/bin/perl

use strict;

my ($file1, $file2) = @ARGV;

my %dp4;
open IN, "$file1" or die $!;
while(<IN>) {
	chomp $_;
	my @items = split("\t", $_);
	my $key = $items[0] . "-" . $items[1];
	$dp4{$key} = \@items;
}
close IN;

open IN, "$file2" or die $!;
while(<IN>) {
        chomp $_;
        my @items = split("\t", $_);
        my $key = $items[0] . "-" . $items[1];
	if (exists $dp4{$key}) {
		if($items[4] eq "0/0" & $items[5] eq "1/1") {	
			my $dp_pos = get_dp_pos($items[2]);
			my $dp_pos_minor = get_dp_pos($items[3]);
			print $items[0] . "\t" . $items[1] . "\t" . $items[2] . "\t" . @{$dp4{$key}}[$dp_pos] . "\t" . $items[3] . "\t" . @{$dp4{$key}}[$dp_pos_minor] . "\n";
		}
		if($items[5] eq "0/0" & $items[4] eq "1/1") {
                        my $dp_pos = get_dp_pos($items[3]);
                        my $dp_pos_minor = get_dp_pos($items[2]);
                        print $items[0] . "\t" . $items[1] . "\t" . $items[3] . "\t" . @{$dp4{$key}}[$dp_pos] . "\t" . $items[2] . "\t" . @{$dp4{$key}}[$dp_pos_minor] . "\n";
                }
	}
}
close IN;

# allele order A_DP,C_DP,G_DP,T_DP,N_DP
sub get_dp_pos {
	my $allele = shift;
	if($allele eq "A") {
		return 4;
	}
	 if($allele eq "C") {
                return 5;
        }
	 if($allele eq "G") {
                return 6;
        }
 	if($allele eq "T") {
                return 7;
        }
}

