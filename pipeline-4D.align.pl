#!/usr/bin/perl -w
use strict;
use List::MoreUtils qw(uniq);

my %code=(
"GCT"=>"A","GCC"=>"A","GCA"=>"A","GCG"=>"A",
"CGT"=>"R","CGC"=>"R","CGA"=>"R","CGG"=>"R","AGA"=>"R","AGG"=>"R",
"AAT"=>"N","AAC"=>"N",
"GAT"=>"D","GAC"=>"D",
"TGT"=>"C","TGC"=>"C",
"GAA"=>"E","GAG"=>"E",
"CAA"=>"Q","CAG"=>"Q",
"GGT"=>"G","GGC"=>"G","GGA"=>"G","GGG"=>"G",
"CAT"=>"H","CAC"=>"H",
"ATT"=>"I","ATC"=>"I","ATA"=>"I",
"TTA"=>"L","TTG"=>"L","CTT"=>"L","CTC"=>"L","CTA"=>"L","CTG"=>"L",
"AAA"=>"K","AAG"=>"K",
"ATG" => "M",
"TTT" => "F","TTC" => "F",
"CCT" => "P","CCC" => "P","CCA" => "P","CCG" => "P",
"TCT" => "S","TCC" => "S","TCA" => "S","TCG" => "S","AGT" => "S","AGC" => "S",
"ACT" => "T","ACC" => "T","ACA" => "T","ACG" => "T",
"TGG" => "W",
"TAT" => "Y","TAC" => "Y",
"GTT"=>"V","GTC"=>"V","GTA"=>"V","GTG"=>"V",
"TAG"=>"-","TGA"=>"-","TAA"=>"-",
);

my $line="";
my $i="";
my $z="";
my $out="";
my $j="";
my $g="";
my $h="";
my $c=0;
my $spe="";
my $tmp="";
my $tmp1="";
my %site=();
my @pos=();
my @len=();
my %cds=();
my %ids=();
my $OG="";
my $file="";
my @emit=();
my $count=0;
my $id="";
my @column=();
my @invsplit=();
my %tag=();
my %tagog=();
my $inver="";
my %seq=();
my %inversion=(
        "inv1"=>"Tan_scaffold_1-1-1795769",
        "inv2"=>"Tan_scaffold_2-1-1714778",
        "inv3"=>"Tan_scaffold_5-28904837-30040300",
        "inv4"=>"Tan_scaffold_15-24077650-26558006",
        "inv5"=>"Tan_scaffold_16-1-2026078",
);
open IN,"Trachinotus_anak.gff3";
while($line=<IN>){
        chomp($line);
        if($line=~/mRNA/){
                ($tmp,$i)=($line=~/ID=FUN_(.+)-(T\d+)\;/);
		$id=$tmp.".".$i;
		#print STDERR $id,"\n";
		$tag{$id}="non";
                @column=split(/\t/,$line);
		#$region{$id}=$column[0]."-".$column[3]."-".$column[4];
		foreach $i (keys %inversion){
			@invsplit=split(/-/,$inversion{$i});
			if($column[0] eq $invsplit[0] and $column[3]>=$invsplit[1] and $column[4]<=$invsplit[2]){
				$tag{$id}=$i;
				#$print STDERR $id,"\t",$tag{$id},"\n";
				#		last;
			}
		}
		#print STDERR $id,"\t",$tag{$id},"\n";
        }

}
close IN;
my @alignment=glob("./Single_Copy_Orthologue_Sequences/*fa.nt.cleanup.aln");
foreach $i (@alignment){
	($OG)=($i=~/Sequences\/(.+)\.fa/);
	$c=0;
	open IN,"$i";
		while($line=<IN>){
			if($line=~/>/){
				$c++;
			}
			elsif($line=~/[RYMKSWHBVDN]/i){
					push(@emit,$OG);
			}
		}
	close IN;
	if($c<2){
		push(@emit,$OG);
	}
}

foreach $file (@alignment){
	$count=0;
	foreach $i (@emit){
                if($file=~/$i/){
                        $count++;
                        last;
                }
        }
	if($count==0){
		open IN,$file or die "$!";
		($OG)=($file=~/Sequences\/(.+)\.fa/);
#		print STDERR "open $file\nProcessing\t$i\n";
		while($line=<IN>){
			chomp $line;
			if($line=~/>/){
				#		($spe,$tmp)=($line=~/>(.+)_(.+)\./);
#				print STDERR $spe,"\t",$tmp,"\n";
				#	$ids{$spe}{$OG}=$tmp;
				if($line=~/>Tan/){
					$spe="Tan";
					($tmp)=($line=~/>Tan(.+)/);
                                	$tagog{$OG}=$tag{$tmp};
                        	}
				elsif($line=~/>Tbl/){
					$spe="Tbl";
				}
			}
			else{
				print STDERR $OG,"\t",$tagog{$OG},"\n";
				$cds{$tagog{$OG}}{$spe}{$OG}.=uc($line);
			}
		}
	}
	else{
		print STDERR "$file is omitted\n";
	}
}

#open OUT,">4D.all.fa";
#foreach $i (sort {$a cmp $b} keys %cds){
#	print OUT ">",$i,"\n";
#	foreach $j (sort {$a cmp $b} keys %{$cds{$i}}){
#		print OUT $cds{$i}{$j};
#	}
#	print OUT "\n";
#}
#close OUT;
foreach $z (keys %cds){
  $out="4D"."\.".$z.".fa";
  print STDERR "Next step: extract 4Dsite for $z\n";
  #open OUT,">$out";
  @len=();
  %site=();
  foreach $i (sort {$a cmp $b} keys %{$cds{$z}}){
	$tmp1="";
	#print OUT ">",$i,"\n";
	#	print STDERR "Processing species:",$i,"\n";
	foreach $j (sort {$a cmp $b} keys %{$cds{$z}{$i}}){
		#print STDERR $i,"\t",$j,"\n",$cds{$i}{$j},"\n";
#		print STDERR $id{$i}{$j},"\n";
		$tmp=$cds{$z}{$i}{$j};
		$tmp1.=&D4site($tmp);
	}
	print STDERR length($tmp1),"\n";
	push(@len,length($tmp1));
	@{$site{$i}}=split(//,$tmp1);
  }
  #close OUT;

  my @uniqlen=uniq @len;
  my @indN=keys %site;
  @pos=();
  print STDERR "The number of sequence length:",$#uniqlen+1,"\n";
  if($#uniqlen==0){
	foreach $i (0..$uniqlen[0]-1){
		$c=0;
		foreach $j (sort {$a cmp $b} keys %site){
			if(${$site{$j}}[$i]=~/[ATGC]/){
				$c++;
			}
		}
	#	print STDERR $c,"\t",$#indN+1,"\t",$i,"\n";
		if($c==2){
			push(@pos,$i);
		}
	}
  }
  open OUT, ">$out";
  %seq=();
  foreach $i (keys %site){
	print OUT ">",$i,"\n";
	foreach $j (sort @pos){
		print OUT ${$site{$i}}[$j];
		$seq{$i}.=${$site{$i}}[$j];
	}
	print OUT "\n";
  }
  $c=0;
  foreach $i (@pos){
  	if(${$site{"Tan"}}[$i] ne ${$site{"Tbl"}}[$i]){
		$c++;
	}
  }
  print STDERR $z,"\t",length($seq{"Tan"}),"\t",length($seq{"Tbl"}),"\t",$c,"\t",$c/length($seq{"Tan"}),"\n";
  close OUT;
}
sub D4site{
	my ($tmp)=@_;
	my @col=split(//,$tmp);
	my @tmparray=();
	my @newtmparray=();
	my $seq="";
#	print STDERR $tmp,"\n";
	for($j=0;$j<$#col-3;$j+=3){
		if($col[$j]=~/-/ and $col[$j+1]=~/-/ and $col[$j+2]=~/-/){
			$seq.="---";
		}
		elsif($col[$j]=~/[ATGC]/ and $col[$j+1]=~/[ATGC]/ and $col[$j+2]=~/[ATGC]/ and length($code{$col[$j].$col[$j+1].$col[$j+2]})>0){
			@tmparray=();
			for $h ("A","T","G","C"){
				if(length($code{$col[$j].$h.$col[$j+2]})<1){next}
				else{
					push(@tmparray,$code{$h.$col[$j+1].$col[$j+2]});
				}
			}
			@newtmparray = uniq(@tmparray);
			if($#newtmparray==0){
				$seq.=$col[$j]
			}
			else{
				$seq.="-";
			}
			###########
			@tmparray=();
			for $h ("A","T","G","C"){
				if(length($code{$col[$j].$h.$col[$j+2]})<1){next}
				else{
					push(@tmparray,$code{$col[$j].$h.$col[$j+2]});
				}
			}
			@newtmparray = uniq(@tmparray);
			if($#newtmparray==0){
				$seq.=$col[$j+1]
			}
			else{
				$seq.="-"
			}
			#############
			@tmparray=();
			for $h ("A","T","G","C"){
				if(length($code{$col[$j].$col[$j+1].$h})<1){next}
				else{
					push(@tmparray,$code{$col[$j+0].$col[$j+1].$h});
				}
			}
			@newtmparray =uniq(@tmparray);
			if($#newtmparray==0){
				$seq.=$col[$j+2]
			}
			else{
				$seq.="-"
			}
                }
		else{
			print STDERR $tmp,"\n";
			print STDERR $j,"\t",$col[$j],$col[$j+1],$col[$j+2],"\n";
			die "there is unaligned cds!\n";
		}
	}
	return($seq);
}









