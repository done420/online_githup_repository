#!/usr/bin/perl -w
use warnings;
use strict;

die "USAGE:perl $0 <GTF> <tophat.bed> <ensembfile> <output>\n" unless @ARGV==4;

#DATE:2014.12.15 Fix basespace bed promblem

my $gtf=$ARGV[0];
my $bed=$ARGV[1];
my $output=$ARGV[3];
my $ensemble = $ARGV[2];
my %gene_wide;
&main();

sub bed(){
	open (IN,"$bed") or die "$!";
	my $line;
	my %bed;
	#skip title
	$line=<IN>;
	while($line=<IN>){
		$line=~tr/\n\r//d;
		my @ele=split("\t",$line);
		my $chro=$ele[0];
			 $chro =~s/chr//g;
		my @info1=split("\,",$ele[10]);
		my @info2=split("\,",$ele[11]);
		my $mark=($ele[1]+$info1[0])."_".($ele[1]+$info2[1]+1)."_".$ele[5];
		if (not(exists($bed{$chro}{$mark}))){
			  $bed{$chro}{$mark}=$chro."\t".($ele[1]+$info1[0])."\t".($ele[1]+$info2[1]+1)."\t".$ele[4]."\t".$ele[5];
			}
			else{
				print $chro,"\t",$mark,"\n";
				}
		}
	close IN or die "$!";
	return %bed;
	}

sub gtf(){
	open (IN,"$gtf") or die "$!";
	my $line;
	my %gtf;
	while($line=<IN>){
		$line=~tr/\n\r//d;
		my @ele=split("\t",$line);
		if ($ele[2]){
		if ($ele[2]=~/exon/){
		    my $chro=$ele[0];
		    my $mark=$ele[3]."\t".$ele[4]."\t".$ele[6]."\t".$ele[0];
		    my @info=split("\;",$ele[8]);
		    $info[0]=~/(.+?)\"(.+?)\"/;
		    my $geneid=$2;
		    $info[1]=~/(.+?)\"(.+?)\"/;
		    my $transid=$2;
		    $info[2]=~/(.+?)\"(.+?)\"/;
		    my $exonnu=$2;
		    my $value=$mark."\t".$exonnu;
		   push @{$gtf{$geneid}{$transid}},$value;
		} 
	}
}
	close IN or die "$!";
	return %gtf;
	}
	
sub splice_site(){
	my %gtf=&gtf();
	my %site;
	foreach my $geneid (keys(%gtf)){
		foreach my $transid (keys(%{$gtf{$geneid}})){
			 my @array=();
			    @array=@{$gtf{$geneid}{$transid}};
			    if (@array==1){
			    	  next;
			    	}
			    my $i;
			    for ($i=0;$i<$#array;$i++){
			    	   my @ele1=split("\t",$array[$i]);
			    	   my @ele2=split("\t",$array[$i+1]);
			    	      my $chro=$ele1[3];
			    	      my $mark;
			    	      if ($ele1[2] eq "+"){
			    	          $mark=$ele1[1]."_".$ele2[0]."_".$ele1[2];
			    	      }
			    	      else{
			    	      	  $mark=$ele2[1]."_".$ele1[0]."_".$ele1[2];
			    	      	}
			    	   if (not(exists($site{$chro}{$mark}))){
			    	   	   $site{$chro}{$mark}=$geneid."_".$transid."#".$ele1[4].":".$ele2[4];
			    	   	}
			    	   	else{
			    	   		$site{$chro}{$mark}.="\t".$geneid."_".$transid."#".$ele1[4].":".$ele2[4];
			    	   		}
			    	}
			}
		}
	return %site;
	}

sub gene_wide(){
	open (IN,"$ensemble") or die "$!";
	my $line;
	$line=<IN>;
	while($line=<IN>){
		$line=~tr/\n\r//d;
		my @ele=split("\t",$line);
		   if (@ele == 6){
		   push @{$gene_wide{$ele[0]}},@ele[1..5];
		   }
		   else{
		   	push @{$gene_wide{$ele[0]}},@ele[1..4];
		   	push @{$gene_wide{$ele[0]}}, "NULL";
		   	}
		}
	close IN or die "$!";
	}

sub judge(){
	my $start=$_[0];
	my $end=$_[1];
	my $chro=$_[2];
	foreach my $geneid (keys(%gene_wide)){
		 my $gchro=$gene_wide{$geneid}[0];
		 my $ges=$gene_wide{$geneid}[2]-90;               #span 10 bp at least
		 my $gen=$gene_wide{$geneid}[3]+90;               #span 10 bp at least
		 my $genename=$gene_wide{$geneid}[4];
		 if ($gchro eq $chro){
		 	  if (($ges<=$start) and ($gen>=$end)){
		 	  	   my $mark=$geneid."\t".$genename;
		 	  	   return $mark;
		 	  	   last;
		 	  }
		 	}
		}
	}


sub main(){
 my %bed=&bed();
 my %site=&splice_site();
 
  &gene_wide();
	open (OUT1,">$output.new") or die "$!";
	open (OUT2,">$output.old") or die "$!";
	foreach my $key1 (sort(keys(%bed))){
		foreach my $key2 (sort(keys(%{$bed{$key1}}))){
			if (exists($site{$key1}{$key2})){
				  print OUT2 ($bed{$key1}{$key2},"\t",$site{$key1}{$key2},"\n");
				}
				else{
					my ($start,$end,$strand)=split("\_",$key2);
					my $ano=&judge($start,$end,$key1);
					print OUT1 ($bed{$key1}{$key2},"\t",$ano,"\n");
					}
			}
		}
	close OUT1 or die "$!";
	close OUT2 or die "$!";
	}