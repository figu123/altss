#!/usr/bin/perl
#Last changed Time-stamp: <2019-10-09 17:52:28 amrgalal>

#liest ein bed file fuer annotierte TSS ein und gibt aus wo die TSS liegen

#usage: perl ../scripts/TSS_classification.pl --ptt ../annotation/NC_000915/NC_000915.ptt --rnt ../annotation/NC_000915/NC_000915.rnt -i files/TSS_c.bed -o files/TSS_classified.table --workdir '.'

use lib '/mnt/storage/progs/TSSAR/local/lib/perl5';
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my $call = join(" ", $0, @ARGV);

my (%GENES_PLUS, %GENES_MINUS, %TSS_PLUS, %TSS_MINUS, %TRANSSTARTS_PLUS, %TRANSSTARTS_MINUS, %SCORE)=();
my $max_utr_length=250;

my ($ptt_file, $rnt_file, $door_file)=();
my ($i_fh, $o_fh)=('', 'TSS_classified.table');
my $workdir ='.';
my $tmpdir ='.';
my %DataFrameData=();
my $smear_fh='files/smear.table';
my %smear_data=();

my $logfile = "TSSAR_wv.log";   # TSSAR log file
my @log = ();                   # TSSAR log array

GetOptions(
  "ptt:s"       => \$ptt_file,
  "rnt:s"       => \$rnt_file,
  "door:s"      => \$door_file,
  "i:s"         => \$i_fh,
  "o:s"         => \$o_fh,
  "smear_i:s"     => \$smear_fh,
  "workdir:s"   => \$workdir
    );

# append / to workdir and tmpdir if necessary
$workdir .= "/" unless ($workdir =~ m/\/$/);
$tmpdir  .= "/" unless ($tmpdir =~ m/\/$/);
$smear_fh           = $workdir.$smear_fh;

# check log file
my $lf = "${workdir}$logfile";
unless (-e $lf ) {
  open(LOGFILE, ">${workdir}$logfile") or warn "Cannot create logfile ${workdir}$logfile";
}
else {
  open(LOGFILE, ">>${workdir}$logfile") or warn "Cannot open logfile ${workdir}$logfile for appending";
}

# redirect STDERR to $logfile
open(STDERR, ">>", "${workdir}$logfile");

print LOGFILE "TSS_classification.pl    ...  call: $call\n";

#einlesen der info wie verschmiert eine clustered TSS position ist
print LOGFILE "TSS_classification.pl    ...  start reading $smear_fh\n";
open IN, "< $smear_fh" or die "can t open $smear_fh for reading\n";
while(<IN>){
  chomp;
  my @F=split/\t/,$_;
  $smear_data{$F[1]}{$F[0]}=$F[2];
}
close IN;    

# einlesen des TSS.gff von STDIN
print LOGFILE "TSS_classification.pl    ...  start reading ${workdir}$i_fh\n";
open IN, "< ${workdir}$i_fh" or die "can t open file ${workdir}$i_fh\n";
while(<IN>){
  chomp;
  my @F=split"\t", $_;
  
  my ($diff,$pval,$id)=();
  if ($F[8]=~m/DIFF=(\d+);/){
    $diff=$1;
  }
  if ($F[8]=~m/PVAL=(.+)$/){
    $pval=$1;
    $pval=~s/NA/1/gi;
  }
  if ($F[8]=~m/ID=(TSS_\d+)/){
    $id=$1;
  }
  
  @{$DataFrameData{$F[6]}{$F[4]}}=($pval, $diff);
  $SCORE{$id}=$F[5];

  if ($F[6] eq '+'){
    $TSS_PLUS{$F[4]}=$id;
  }
  elsif ($F[6] eq '-'){
    $TSS_MINUS{$F[4]}=$id;
  }
  else {
    print STDERR "Fehler in gff-file line: $.\n";
  }
}
close IN;

#einlesen der door operon files
if ($door_file) {
  print LOGFILE "TSS_classification.pl    ...  start reading door annotation ($door_file)\n";
  open GFF, "< $door_file" or die "can t open  $door_file\n";
  while(<GFF>) {
    chomp;
    next if ( m/^$/);
    
    my @F =split /\t/, $_;
    my ($seq, $method, $domain, $Start, $Stop, $foo, $Strand, $bar, $comment)= @F;
    $comment=~m/id="([mpo_]{3}\d+)"/;
    my $Id=$1;

    if ($Strand eq '+'){
      $TRANSSTARTS_PLUS{$Start}=$Id;
      foreach my $pos ($Start+1 .. $Stop-1) {
	push @{$GENES_PLUS{$pos}}, $Id;
      }
    }
    if ($Strand eq '-') {
      $TRANSSTARTS_MINUS{$Stop}=$Id;
      foreach my $pos ($Start+1 .. $Stop-1) {
	push @{$GENES_MINUS{$pos}}, $Id;
      }
    }
  }
 close GFF; 
}

#einlesen des ncbi ptt file
if ($ptt_file) {
  my $los=0;
  print LOGFILE "TSS_classification.pl    ...  start reading ptt annotation ($ptt_file)\n";
  open PTT, "< $ptt_file" or die "can t open  $ptt_file\n";
  while(<PTT>) {
    chomp;
    if (!$los || m/^$/) {
      if (m/^Location/) {
	$los=1;
	next;
      }
      else {
	next;
      }
    }
    
    my @F =split /\t/, $_;
    my ($Location, $Strand, $Length, $PID, $Gene, $Synonym, $Code, $COG, $Product) = @F;
    $Location=~m/(\d+)\.\.(\d+)/;
    my ($Start, $Stop)=($1,$2);

    if ($Strand eq '+'){
      $TRANSSTARTS_PLUS{$Start}=$Synonym;
      foreach my $pos ($Start+1 .. $Stop-1) {
	push @{$GENES_PLUS{$pos}}, $Synonym;
      }
    }
    if ($Strand eq '-') {
      $TRANSSTARTS_MINUS{$Stop}=$Synonym;
      foreach my $pos ($Start+1 .. $Stop-1) {
	push @{$GENES_MINUS{$pos}}, $Synonym;
      }
    }
  }
 close PTT; 
}

#einlesen des ncbi rnt file
if ($rnt_file) {
  my $los=0;
  print LOGFILE "TSS_classification.pl    ...  start reading rnt annotation ($rnt_file)\n";
  open RNT, "< $rnt_file" or die "can t open  $rnt_file\n";
  while(<RNT>) {
    chomp;
    if (!$los || m/^$/) {
      if (m/^Location/) {
	$los=1;
	next;
      }
      else {
	next;
      }
    }
    
    my @F =split /\t/, $_;
    my ($Location, $Strand, $Length, $PID, $Gene, $Synonym, $Code, $COG, $Product) = @F;
    $Location=~m/(\d+)\.\.(\d+)/;
    my ($Start, $Stop)=($1,$2);

    if ($Strand eq '+'){
      $TRANSSTARTS_PLUS{$Start}=$Synonym;
      foreach my $pos ($Start+1 .. $Stop-1) {
	push @{$GENES_PLUS{$pos}}, $Synonym;
      }
    }
    if ($Strand eq '-') {
      $TRANSSTARTS_MINUS{$Stop}=$Synonym;
      foreach my $pos ($Start+1 .. $Stop-1) {
	push @{$GENES_MINUS{$pos}}, $Synonym;
      }
    }
  }
 close RNT; 
}

unless ($rnt_file || $ptt_file || $door_file) {
  print LOGFILE "TSS_classification.pl    ...  no annotation file specified (door, ptt nor rnt)\n";
}

print LOGFILE "TSS_classification.pl    ...  start writing output (${workdir}$o_fh)\n";
open OUT, "> ${workdir}$o_fh" or die "can t write to file ${workdir}$o_fh\n";

foreach my $pos (sort {$a<=>$b} keys %TSS_PLUS ) {

    my $status='';
    my $comment='';
    my @I_gene='';
    
    #checken ob es in einem annotierten gene liegt (I)
    if (defined $GENES_PLUS{$pos}) {
      $status.='I';
      $comment.='within gene(s) '.join(" ", @{$GENES_PLUS{$pos}}).'; ';
      @I_gene=@{$GENES_PLUS{$pos}};
    }
    
    #checken obs $max_utr_length upstream von translationstart liegt (P)
    foreach my $p ($pos .. $pos+$max_utr_length) {
      if ( defined $TRANSSTARTS_PLUS{$p} ){
	#next if (grep {m/$TRANSSTARTS_PLUS{$p}/} @I_gene);
	$status.='P';
	$comment.=$p-$pos."nt upstream of gene $TRANSSTARTS_PLUS{$p}; ";
	last;
      }
    }
    
    #check ob antisense
    if (defined $GENES_MINUS{$pos}) {
      $status.='Ai';
      $comment.='antisense to gene(s) '.join(" ", @{$GENES_MINUS{$pos}}).'; ';
    }
    elsif (defined $GENES_MINUS{$pos+30}) {
      my $au_dist_found=0;
      for my $au_dist_search (1 .. 30) {
	if ($au_dist_found==0 && defined $GENES_MINUS{$pos+$au_dist_search}) {
	  $au_dist_found=$au_dist_search;
	  $status.='Ad';
	  $comment.='antisense to gene(s) '.join(" ", @{$GENES_MINUS{$pos+$au_dist_found}})." (${au_dist_found}nt downstream); ";
	  last;
	}
      }
    }
    #else "O" zuweisen
    if ($status eq '') {
      $status='O';
      $comment='-';
    }
 
    print OUT "$pos\t+\t$TSS_PLUS{$pos}\t$SCORE{$TSS_PLUS{$pos}}\t${$DataFrameData{'+'}{$pos}}[1]\t${$DataFrameData{'+'}{$pos}}[0]\t$smear_data{'+'}{$pos}\t$status\t$comment\n";
}

foreach my $pos (sort {$a<=>$b} keys %TSS_MINUS ) {
    
    my $status='';
    my $comment='';
    my @I_gene='';
    
    #checken ob es in einem annotierten gene liegt (I)
    if (defined $GENES_MINUS{$pos}) {
      $status.='I';
      $comment.='within gene(s) '.join(" ", @{$GENES_MINUS{$pos}}).'; ';
      @I_gene=@{$GENES_MINUS{$pos}};
    }
    #checken obs $max_utr_length upstream von translationstart liegt (P)
    foreach my $p ( reverse($pos-$max_utr_length .. $pos) ) {
      if ( defined $TRANSSTARTS_MINUS{$p} ){
	#next if (grep {m/$TRANSSTARTS_MINUS{$p}/} @I_gene);
	$status.='P';
	$comment.=$pos-$p."nt upstream of gene $TRANSSTARTS_MINUS{$p}; ";
	last;
      }
    }
    
    #check ob antisense
    if (defined $GENES_PLUS{$pos}) {
      $status.='Ai';
      $comment.='antisense to gene(s) '.join(" ", @{$GENES_PLUS{$pos}}).'; ';
    }
    elsif (defined $GENES_PLUS{$pos-30}) {
      my $au_dist_found=0;
      for my $au_dist_search (1 .. 30) {
	if ($au_dist_found==0 && defined $GENES_PLUS{$pos-$au_dist_search}) {
	  $au_dist_found=$au_dist_search;
	  $status.='Ad';
	  $comment.='antisense to gene(s) '.join(" ", @{$GENES_PLUS{$pos-$au_dist_found}})." (${au_dist_found}nt downstream); ";
	}
      }
    }
    #else "O" zuweisen
    if ($status eq '') {
      $status='O';
      $comment='-';
    }
    
    print OUT "$pos\t-\t$TSS_MINUS{$pos}\t$SCORE{$TSS_MINUS{$pos}}\t${$DataFrameData{'-'}{$pos}}[1]\t${$DataFrameData{'-'}{$pos}}[0]\t$smear_data{'-'}{$pos}\t$status\t$comment\n";
}

print LOGFILE "TSS_classification.pl    ...  DONE\n";
