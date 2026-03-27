#!/usr/bin/perl
#Last changed Time-stamp: <2019-10-09 17:32:49 amrgalal>

# liest bed file und clustered consecutive TSS zusammen
#usage: perl bed_cluster.pl --score 'd' --i TSS.bed --o TSS_c.bed --workdir '.' --tmpdir '.' --range 2

use lib '/mnt/storage/progs/TSSAR/local/lib/perl5';
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Pod::Usage;
use File::Basename;
use File::Temp qw(tempfile);
use vars qw//;

my $call = join(" ", $0, @ARGV);

# set paths and i/o filenames
my $workdir  ='.';
my $tmpdir   ='.';
my $i_fh     = 'TSS.bed';
my $o_fh     = 'TSS_c.bed';
my $smear_fh = 'smear.table';
my (%smear_plus, %smear_minus)=();
my @out_data = ();

my $score_mode = 'd';
my $verbose = 0;

my $logfile = "TSSAR_wv.log";   # TSSAR log file
my @log = ();                   # TSSAR log array

## Parameter (how far ist consecutive?)
my $consecutiv_range =1;

# parse comand-line arguments
pod2usage(-verbose => 0)
    unless GetOptions(
      "verbose"     => \$verbose,
      "score:s"     => \$score_mode,
      "range:i"     => \$consecutiv_range,
      "workdir:s"   => \$workdir,
      "tmpdir:s"    => \$tmpdir,
      "i:s"         => \$i_fh,
      "o:s"         => \$o_fh,
      "smear_o:s"     => \$smear_fh,
      "man"         => sub{pod2usage(-verbose => 2)},
      "help|?"      => sub{pod2usage(-verbose => 1)},
    );

# append / to workdir and tmpdir if necessary
$workdir .= "/" unless ($workdir =~ m/\/$/);
$tmpdir  .= "/" unless ($tmpdir =~ m/\/$/);

$i_fh     = $workdir.$i_fh;
$o_fh     = $workdir.$o_fh;
$smear_fh = $workdir.$smear_fh;

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

print LOGFILE "gff_cluster.pl           ...  call: $call\n";

my (%TSS_PLUS, %TSS_MINUS);
my (@coverage_plus, @coverage_minus);
  
# read in TSS in bed formate
print LOGFILE "gff_cluster.pl           ...  start reading $i_fh\n";
open TSSIN, "< $i_fh" or die "can t read file $i_fh\n";
while(<TSSIN>){
  chomp;
  next if(m/^\s*$/);
  my @F=split "\t" ,$_;
  
  my $strand = $F[6]; #$ F 6
  my $pos = $F[4]; #F 4
  my $diff = '10';
  if ($F[8]=~m/DIFF=(\d+)/){
    $diff=$1;
  }

  if ($strand eq '+'){
    if ( exists($TSS_PLUS{$pos}) ){
      
	@{$TSS_PLUS{$pos}}   = @F if ($coverage_plus[$pos] < $diff);
	$coverage_plus[$pos] = $diff if ($coverage_plus[$pos] < $diff);
      
    }
    else {
      @{$TSS_PLUS{$pos}}   = @F;
      $coverage_plus[$pos] = $diff;
    }
  }
  
  elsif($strand eq '-') {
    if ( exists($TSS_MINUS{$pos}) ){

	@{$TSS_MINUS{$pos}}   = @F if ($coverage_minus[$pos] < $diff);
	$coverage_minus[$pos] = $diff if ($coverage_minus[$pos] < $diff);
      
    }
    else {
      @{$TSS_MINUS{$pos}}   = @F;
      $coverage_minus[$pos] = $diff;
    }
  }
}

close TSSIN;


print LOGFILE "gff_cluster.pl           ...  start analysing data (+ strand)\n";
# write clustered TSS
# cluster plus strand
my @peak=();
my @pos=();

foreach my $pos (0 .. $#coverage_plus) {
  if ( defined($coverage_plus[$pos]) ) {
    push @peak, $coverage_plus[$pos];
    push @pos, $pos;
  }
  elsif ( exists($peak[0]) ){
    next if ( ($pos-$pos[-1]) <= $consecutiv_range );
    my $min_pos = '';
    $min_pos = ($score_mode eq 'p') ? &smallestpval_pos(@peak) : &highestcoverage_pos(@peak);
    
    push @out_data,  join("\t", @{$TSS_PLUS{$pos[$min_pos]}});
    $smear_plus{$pos[$min_pos]}=abs($pos[-1]-$pos[0]);
    @peak=();
    @pos=();
  }
  else {
    next;
  }
}
if ( exists($peak[0]) ){
  my $min_pos = '';
  $min_pos = ($score_mode eq 'p') ? &smallestpval_pos(@peak) : &highestcoverage_pos(@peak);
    
  push @out_data,  join("\t", @{$TSS_PLUS{$pos[$min_pos]}});
  $smear_plus{$pos[$min_pos]}=abs($pos[-1]-$pos[0]);
  @peak=();
  @pos=();
}

print LOGFILE "gff_cluster.pl           ...  start analysing data (- strand)\n";  
#cluster minus strand
@peak=();
@pos=();
foreach my $pos (0 .. $#coverage_minus) {
  if ( defined($coverage_minus[$pos]) ) {
    push @peak, $coverage_minus[$pos];
    push @pos, $pos;
  }
  elsif ( exists($peak[0]) ){
    next if ( ($pos-$pos[-1]) <= $consecutiv_range );
    my $min_pos = '';
    $min_pos = ($score_mode eq 'p') ? &smallestpval_pos(@peak) : &highestcoverage_pos(@peak);
    
    push @out_data,  join("\t", @{$TSS_MINUS{$pos[$min_pos]}});
    $smear_minus{$pos[$min_pos]}=abs($pos[-1]-$pos[0]);
    @peak=();
    @pos=();
  }
  else {
    next;
  }
}
if ( exists($peak[0]) ){
  my $min_pos = '';
  $min_pos = ("$score_mode" eq "p") ? &smallestpval_pos(@peak) : &highestcoverage_pos(@peak);
  
  push @out_data, join("\t", @{$TSS_MINUS{$pos[$min_pos]}});
  $smear_minus{$pos[$min_pos]}=abs($pos[-1]-$pos[0]);
  @peak=();
  @pos=();
}

print LOGFILE "gff_cluster.pl           ...  start writing output to $o_fh\n";
open OUT, ">$o_fh" or die "can t write to $o_fh\n";
for my $l ( @out_data ) {
  print OUT "$l\n";
}
close OUT;

print LOGFILE "gff_cluster.pl           ...  start writing smear stat to $smear_fh\n";
open OUT, ">$smear_fh" or die "can t write to $smear_fh\n";
for my $l ( keys %smear_plus ) {
  print OUT "$l\t+\t$smear_plus{$l}\n";
}
for my $l ( keys %smear_minus ) {
  print OUT "$l\t-\t$smear_minus{$l}\n";
}
close OUT;


print LOGFILE "gff_cluster.pl           ...  DONE\n";

# subroutines
sub highestcoverage_pos {
  my @arr = @_;
  my $max_pos=0;
  my $max_dist=0;
  
  foreach my $pos (0 .. $#arr) {
    if ($arr[$pos]>$max_dist) {
      $max_dist=$arr[$pos];
      $max_pos=$pos;
    }
  }
  return $max_pos;
}

sub smallestpval_pos {
  my @arr = @_;
  my $min_pos=0;
  my $min_pval=9999;

  foreach my $pos (0 .. $#arr) {
    if ($arr[$pos]<$min_pval) {
      $min_pval=$arr[$pos];
      $min_pos=$pos;
    }
  }
  return $min_pos;
}
