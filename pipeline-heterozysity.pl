#!/usr/bin/perl -w
use strict;


my $line="";
my @col=();
my @ele=();
my @cell=();
my $i="";
my $j="";
my $t1="";
my $t2="";
my %ind=();
my %count=();
my $countline=0;
open IN, "gunzip -c $ARGV[0]|";
while($line=<IN>){
	chomp($line);
	if($line=~/##/){next}
	elsif($line=~/#CHR/){
		@col=split(/\t/,$line);
		foreach $i (9 .. $#col){
			$ind{$i}=$col[$i];
		}
	}
	else{
		$countline++;
		@col=split(/\t/,$line);
		foreach $i (9 .. $#col){
			($t1,$t2)=($col[$i]=~/^(\d)\/(\d)/);
			if($t1 ne $t2){
				$count{$ind{$i}}++
			}
		}
	}
}

foreach $i (keys %count){
	print $i,"\t",$count{$i}/$countline,"\n";
}
