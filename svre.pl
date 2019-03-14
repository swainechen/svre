#!/usr/bin/perl

# TODO
# DONE reference files seem to be only used for lengths
# DONE we should be able to get these from the bam file header directly
# we should use libraries to open the bam file
# DONE we need to do better at coordinate mapping for multiple chromosomes
# DONE we shouldn't generate unnecessary intermediate files (i.e. join file)
# DONE we need to clean up the output files and formats
# we need to process into predicted SVs better, including all significant hits
# if use single distances, then can't tell whether to add or subtract coords
#   for inversions
# DONE collapse into bins once we set ywin

#@@@@@@@@@@@
# Libraries
#@@@@@@@@@@@
use warnings;
use strict;
use sv;
use Cwd;
use File::Temp;
use Math::Round;
use Getopt::Long;
use List::Util qw(min max sum);

#@@@@@@@@@@@
# Variables 
#@@@@@@@@@@@
my $tempdir = File::Temp::tempdir(CLEANUP => 1);
my $command_line = join (" ", $0, @ARGV); chomp $command_line;
my $read;	# the main data file
my $r1 = "";	# sam/bam file 1
my $r2 = "";	# sam/bam file 2
my $read_length_f = 0;
my $read_length_r = 0;
my $ori = "FR";
my $ywin = 0;
my $optimise_cov = 1;
my $cov_bin = 100;
my $slide = 0;
my $pval = 1;
my $fdr = 0.05;
my $bootstrap = 1000000;
my $start = 0;
my $end = 0;
my $print_dist = 0;
my $continue = 0;
my $quiet = 0;
my $range_cluster = 0;

my $largest_genome = 0;
my $genome_size=0;
my $refh = {};
my $output = "svre-results";
my $output_ic = "";
my $output_detail = "";
my $output_list = "";
my $pval_string = '';
my $r_file = "";
my $now = time;
my $pos = {};
my $both = 0;
my $single = 0;
my $global_dist = {};
my $global_total = 0;
my $count = {};
my $ri = {};
my $loopstart;
my $loopend;
my $use_translocation = 1;
my $R_command = `which Rscript`;
chomp $R_command;

# generic variables we use all the time and others
my @f;
my @r;
my @m;
my $i;
my $j;
my @keys;
my $key;
my $ref;
my $sum;
my %index;
my $dist;
my ($pos1, $pos2);
my ($ref1, $ref2);
my $median_dist;
my @r1;
my $r1_freq;
my @rl;
my @dist;
my $up;
my $bin;
my $pvalue;
my $rcount;
my $bin_size;
my @pair;
my $dcount_for_thisbin;
my $ccount_for_thisbin;
my $success_rate_for_thisbin;
my $pair;
my $var;
my $output_png;
my $sv;
my $s;
my $d;
my @range;
my $rs;
my @dist_total;
my $read_length_total;
my $dist_f1;
my $dist_f2;
my $filter;
my $mapq_min = 30;

#@@@@@@@@@@@@@@@@@
# Program options 
#@@@@@@@@@@@@@@@@@
&Getopt::Long::Configure("pass_through");
GetOptions(
  "r1=s" => \$r1,
  "r2=s" => \$r2,
  "output=s" => \$output,
  "ori=s" => \$ori,
  "translocation!" => \$use_translocation,
  "ywindow=i" => \$ywin,
  "optimise!" => \$optimise_cov,
  "optimize!" => \$optimise_cov,
  "cov=i" => \$cov_bin,
  "slide!" => \$slide,
  "pval!" => \$pval,
  "fdr=f" => \$fdr,
  "bootstrap=i" => \$bootstrap,
  "start=i" => \$start,
  "end=i" => \$end,
  "print!" => \$print_dist,
  "continue!" => \$continue,
  "range=i" => \$range_cluster,
  "R_command=s" => \$R_command,
  "mapq=i" => \$mapq_min,
  "quiet" => \$quiet
);

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Main input file check - otherwise help screen
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
if ($r1 ne "" && -f $r1 && $r2 ne "" && -f $r2) {
  $r1 = Cwd::abs_path($r1);
  $r2 = Cwd::abs_path($r2);
  if ($output eq "") {
    $output = "svre-results";
  }
} else {
  print <<__USAGE__;
Usage: $0 -r1 <R1 bam> -r2 <R2 bam> -ori [FR|FF] [ options ]
Required parameters (default in parentheses):
  -r1 <file>        : single-end mapped forward reads (sam/bam format)	
  -r2 <file>        : single-end mapped reverse reads (sam/bam format)	

Optional parameters:
  -output <name>    : base file name for results ($output)
  -ori [FR|FF]      : FR is for Illumina, FF for Solid) ($ori)

  -mapq <int>       : minimum mapping quality to use (30)
  -translocation|notranslocation
                    : use reads mapping to different chromosomes (translocation)
  -ywindow <int>    : bin size for insert distances
                      (half of insert size rounded to nearest 50)
  -optimize|nooptimize
                    : Use signal-to-noise ratio to pick coverage bin (optimize)
  -cov <int>        : coverage bin, if -nooptimize ($cov_bin)
  -slide|noslide    : per base slide (noslide)

  -pval|nopval      : calculate p-values by resampling (pval)
  -bootstrap <int>  : # of resamplings to do (1000000)
  -fdr <float>      : FDR cut-off, using Benjamini-Yekutieli ($fdr)

  -range <int>      : cluster ranges in output if within this many bp
                      (median read distance)
  -print|noprint    : print the insert sizes between the pairs and exit ($print_dist)
  -start <int>
  -end <int>        : restrict analysis region
                      (defaults: -1 * genome size to genome size)
  -quiet            : no screen output
__USAGE__
  exit;
}

#@@@@@@@@@@@@@@@@@@@@@@@
# Get chromosomes, size
#@@@@@@@@@@@@@@@@@@@@@@@
# if reference files given, parse reference files to get reference names and length
# adjust the length to be divisible by the y-window binning length
my @samheader1 = ();
my @samheader2 = ();
my $mapper = "";
my $mapper_version = "";
my $mapper_command_line = "";

if ($r1 =~ /\.sam$/i) {
  @samheader1 = `samtools view -SH $r1`;
} elsif ($r1 =~ /\.bam$/i) {
  @samheader1 = `samtools view -H $r1`;
} else {
  die "Can't figure out $r1 file type (sam/bam)\n";
}
if ($r1 =~ /\.sam$/i) {
  @samheader2 = `samtools view -SH $r2`;
} elsif ($r1 =~ /\.bam$/i) {
  @samheader2 = `samtools view -H $r2`;
} else {
  die "Can't figure out $r2 file type (sam/bam)\n";
}

foreach $i (@samheader1) {
  if ($i =~ /^\@SQ/) {
    $ref = "";
    $j = 0;
    @f = split /\t/, $i;
    if ($f[1] =~ /SN:(\S+)/) {
      $ref = $1;
    }
    if ($f[2] =~ /LN:(\d+)/) {
      $j = $1;	# i.e. length
    }
    if ($ref ne "" && $j) {
      $refh->{$ref} = $j;
    }
  }
  if ($i =~ /^\@PG/) {
    @f = split /\t/, $i;
    if ($f[1] =~ /ID:(.*?)/) {
      $mapper = $1;
    }
    if ($i =~ /VN:(.*?)\s+/) {
      $mapper_version = $1;
    }
    if ($i =~ /CL:(.*?)/) {
      $mapper_command_line = $1;
    }
  }
}
foreach $i (@samheader2) {
  if ($i =~ /^\@SQ/) {
    $ref = "";
    $j = 0;
    @f = split /\t/, $i;
    if ($f[1] =~ /SN:(\S+)/) {
      $ref = $1;
    }
    if ($f[2] =~ /LN:(\d+)/) {
      $j = $1;	# i.e. length
    }
    if ($ref ne "" && $j) {
     die "SAM headers don't match on line $i\n" if !defined $refh->{$ref} || $refh->{$ref} != $j;
    }
  }
}

#@@@@@@@@@@@@@@@@@@@@@
# Parse sam/bam files
#@@@@@@@@@@@@@@@@@@@@@ 

# For bwa 0.7.5a:
# in bwa single end mapping, the flags just give the information of strand
# - whether the mapping is on the forward (0) or reverse (16) strand
# however, there are tags that give us the information whether the mapping
# is unique, if it has mismatches etc.
# XT:A - type of mapping (unique U, repeat R, etc)
# NM:i - number of mismatches to the reference (edit distance - so will not tell indels)
# X0:i - number of best hits for the read - if multiple, this will be >1
# X1:i - number of suboptimal hits for the read - this is ok even if >1
# XM:i - number of mismatches in the exact alignment of the read - including indels
# -- so test for XT:A:U, NM:i:0, XM:i:0, X0:i:1
#
# For bwa 0.7.10:
# for bwa mem:
# we seem to only get this:
# NM:i - edit distance to reference
# MD:Z - mismatching positions/bases - though if NM:i is 0, this seems to be the same as alignment score
# AS:i - Alignment score
# XS:i - suboptimal alignment score
# -- so test for NM:i:0 and XS:i:0
# for bwa aln (bwa samse):
# still have XT:A, NM, X0, X1, XM - so same tests as above
#
# for bowtie, we only get:
# NM:i - edit distance to reference (mismatches)
# CM:i - edit distance in colorspace
# XA:i - stratum, or # mismatches in seed region
# XM:i - number of alternative alignments
# -- so test for XA:i:0 and NM:i:0
#
# for bowtie2, we get:
# NM:i - edit distance to reference (mismatches)
# XN:i - number of ambiguous bases
# XM:i - mismatches in alignment
# XO:i - number of gap opens
# XG:i - number of gap extensions
# -- so test for XM:i:0, XG:i:0, XO:i:0, NM:i:0
# 

$read = {};	# keyed with read name, then r1/r2, ref/length/pos
		# this will go away when we make $precount
@dist_total = ();
$read_length_total = 0;
$filter = 0;
if (-B $r1) {
  open R, "samtools view $r1|";
} else {
  open R, "<$r1";
}
while (<R>) {
#  if ($mapper eq "bwa" && $mapper_version =~ /^0\.7\.10/) {
  if ($mapper eq "bwa" && $mapper_command_line =~ /bwa\s+mem/) {
    next if $_ !~ /NM\:i\:0/;
    next if $_ !~ /XS\:i\:0/;
#  } elsif ($mapper eq "bwa") {
  } elsif ($mapper eq "bwa" && $mapper_command_line =~ /bwa\s+sam[ps]e/) {
    next if $_ !~ /XT\:A\:U/;
    next if $_ !~ /NM\:i\:0/;
    next if $_ !~ /XM\:i\:0/;
    next if $_ !~ /X0\:i\:1/;
  } elsif ($mapper eq "Bowtie" || $mapper eq "bowtie") {
    next if $_ !~ /XA\:i\:0/;
    next if $_ !~ /NM\:i\:0/;
  } elsif ($mapper eq "Bowtie2" || $mapper eq "bowtie2") {
    next if $_ !~ /XM\:i\:0/;
    next if $_ !~ /XG\:i\:0/;
    next if $_ !~ /XO\:i\:0/;
    next if $_ !~ /NM\:i\:0/;
  } else {
    next if $_ !~ /NM\:i\:0/;
  }
  @f = split /\t/, $_;
  # sam flag 4 is unmapped, 16 is reverse strand map for the read
  next if $f[1] & 4;
  if ($mapq_min > 0) {
    next if $f[4] < $mapq_min;
  }
  $read->{$f[0]}->{r1}->{ref} = $f[2];
  $read->{$f[0]}->{r1}->{length} = length($f[9]);
  if ($f[1] & 16) {
    $read->{$f[0]}->{r1}->{pos} = "-$f[3]";
  } else {
    $read->{$f[0]}->{r1}->{pos} = $f[3];
  }
}
close R;
  
if (-B $r2) {
  open R, "samtools view $r2|";
} else {
  open R, "<$r2";
}
while (<R>){
#  if ($mapper eq "bwa" && $mapper_version =~ /^0\.7\.10/) {
  if ($mapper eq "bwa" && $mapper_command_line =~ /bwa\s+mem/) {
    next if $_ !~ /NM\:i\:0/;
    next if $_ !~ /XS\:i\:0/;
#  } elsif ($mapper eq "bwa") {
  } elsif ($mapper eq "bwa" && $mapper_command_line =~ /bwa\s+sam[ps]e/) {
    next if $_ !~ /XT\:A\:U/;
    next if $_ !~ /NM\:i\:0/;
    next if $_ !~ /XM\:i\:0/;
    next if $_ !~ /X0\:i\:1/;
  } elsif ($mapper eq "Bowtie" || $mapper eq "bowtie") {
    next if $_ !~ /XA\:i\:0/;
    next if $_ !~ /NM\:i\:0/;
  } elsif ($mapper eq "Bowtie2" || $mapper eq "bowtie2") {
    next if $_ !~ /XM\:i\:0/;
    next if $_ !~ /XG\:i\:0/;
    next if $_ !~ /XO\:i\:0/;
    next if $_ !~ /NM\:i\:0/;
  } else {
    next if $_ !~ /NM\:i\:0/;
  }
  @f = split /\t/, $_;
  next if $f[1] & 4;
  if ($mapq_min > 0) {
    next if $f[4] < $mapq_min;
  }
  $read->{$f[0]}->{r2}->{ref} = $f[2];
  $read->{$f[0]}->{r2}->{length} = length($f[9]);
  if ($f[1] & 16) {
    $read->{$f[0]}->{r2}->{pos} = "-$f[3]";
  } else {
    $read->{$f[0]}->{r2}->{pos} = $f[3];
  }
}
close R;

foreach $key (keys %$read) {
  if (defined $read->{$key}->{r1}->{ref} && defined $read->{$key}->{r2}->{ref}) {
    # these should both be well mapped, and all %read fields defined
    $pos1 = $read->{$key}->{r1}->{pos};
    $pos2 = $read->{$key}->{r2}->{pos};
    $ref1 = $read->{$key}->{r1}->{ref};
    $ref2 = $read->{$key}->{r2}->{ref};

    # to filter by reference
    next if (!defined $refh->{$ref1});
    next if (!defined $refh->{$ref2});

    if ($ref1 ne $ref2) {
      next if !$use_translocation;
      # these distances don't contribute to dist_total for median calculation
      $pos->{$ref1}->{$pos1}->{$pos2."___$ref2"}->{count}++;
      $pos->{$ref2}->{$pos2}->{$pos1."___$ref1"}->{count}++;
      push @{$pos->{$read->{$key}->{r1}->{ref}}->{$pos1}->{pair}}, $pos2."___$ref2";
      push @{$pos->{$read->{$key}->{r2}->{ref}}->{$pos2}->{pair}}, $pos1."___$ref1";
      $both++;
      print STDERR $pos1."___$ref1 -- $pos2"."___$ref2\n" if $print_dist;
    } else {
      $read_length_f = $read->{$key}->{r1}->{length};
      $read_length_r = $read->{$key}->{r2}->{length};
      $read_length_total += $read_length_f;
      $read_length_total += $read_length_r;
      $both++;

      # bowtie always reports 0-based (so must add 1) leftmost coordinate based
      # on forward strand, so when a reverse coordinate is reported it is
      # actually the 3' end and not 5' end
      # so must adjust accordingly - we always take the distance b/w the 5'
      # coordinates of the mates.
      if ($pos1 < 0) {
        $pos1 = $pos1 - $read_length_f;
      } else {
        $pos1 = $pos1 + 1;
      }
      if ($pos2 < 0) {
        $pos2 = $pos2 - $read_length_r;
      } else {
        $pos2 = $pos2 + 1;
      }

      $dist = abs(abs($pos1) - abs($pos2));
      # wrap around situations
      # 1. test F < abs(R) if F > 0 - no wrap around
      # pos strand, could have -F>.....-R> or -F>.....<R-
      # 2. test F > abs(R) if F > 0 - need wrap around
      # pos strand -R>.....-F> or <R-.....-F>
      # 3. test abs(F) > abs(R) if F < 0 - no wrap around
      # neg strand, could have <R-.....<F- or -R>.....<F-
      # 4. test abs(F) < abs(R) if F < 0 - need wrap around
      # neg strand <F-.....<R- or <F-.....-R>
      if (($pos1 > 0 && $pos2 > 0) || ($pos1 < 0 && $pos2 < 0)) {
        # we need 2 distances
#        if ($pos1 > 0) {
#          if ($pos1 <= $pos2) {
#            $dist_f1 = $dist;
#            $dist_f2 = $refh->{$read->{$key}->{r1}->{ref}} - $dist;
#          } else {
#            $dist_f1 = $refh->{$read->{$key}->{r1}->{ref}} - $dist;
#            $dist_f2 = $dist;
#          }
#        } else {
#          if (abs($pos1) >= abs($pos2)) {
#            $dist_f1 = $dist;
#            $dist_f2 = $refh->{$read->{$key}->{r1}->{ref}} - $dist;
#          } else {
#            $dist_f1 = $refh->{$read->{$key}->{r1}->{ref}} - $dist;
#            $dist_f2 = $dist;
#          }
#        }
        # to be fair, even though we lose half the coordinates, we limit
        # to no more than half the chromsome - so that this read counts
        # twice, just like the distances for pairs that map to opposite strands
        if ($dist > $refh->{$read->{$key}->{r1}->{ref}}/2) {
          $dist = $refh->{$read->{$key}->{r1}->{ref}} - $dist;
        }
        $dist_f1 = $dist;
        $dist_f2 = $dist;
      } else {
        if ( ($pos1 > 0 && $pos1 > abs($pos2)) ||	# <R-.....-F>
             ($pos1 < 0 && abs($pos1) < abs($pos2)) ) {	# <F-.....-R>
          $dist = $refh->{$read->{$key}->{r1}->{ref}} - $dist;
        }
        # distance should be symmetric
        $dist_f1 = $dist;
        $dist_f2 = $dist;
      }

      if ($ori eq "FF") {
        if (($pos1 > 0 && $pos2 < 0) || ($pos1 < 0 && $pos2 > 0)) {
          # opposite strands - wrong orientation, so negative distance
          $dist_f1 = -$dist_f1;
          $dist_f2 = -$dist_f2;
        }
      } elsif ($ori eq "FR") {
        if (($pos1 > 0 && $pos2 > 0) || ($pos1 < 0 && $pos2 < 0)) {
          # map to same strands - wrong orientation, negative distance
          $dist_f1 = -$dist_f1;
          $dist_f2 = -$dist_f2;
        }
      }

      $pos->{$read->{$key}->{r1}->{ref}}->{$pos1}->{$dist_f1}->{count} += 1;
      push @{$pos->{$read->{$key}->{r1}->{ref}}->{$pos1}->{pair}}, $pos2;

      $pos->{$read->{$key}->{r1}->{ref}}->{$pos2}->{$dist_f2}->{count} += 1;
      push @{$pos->{$read->{$key}->{r1}->{ref}}->{$pos2}->{pair}}, $pos1;

      print STDERR "$dist_f1\n" if $print_dist;
      print STDERR "$dist_f2\n" if $print_dist;
      push @dist_total, $dist_f1;
      push @dist_total, $dist_f2;
    }
  } else {
    # these singles also don't contribute to the total distance
    if (defined $read->{$key}->{r1}->{ref}) {
      $ref1 = $read->{$key}->{r1}->{ref};
      $pos1 = $read->{$key}->{r1}->{pos};
# These could be repeat sequences where mapping doesn't happen due to the
# requirement for unique mapping. Unclear what information they hold.
#      $pos->{$ref1}->{$pos1}->{"0___unmapped"}->{count}++;
#      push @{$pos->{$ref1}->{$pos1}->{pair}}, "0___unmapped";
      print STDERR $pos1."___$ref1 -- 0___unmapped\n" if $print_dist;
      $single++;
    } elsif (defined $read->{$key}->{r2}->{ref}) {
      $ref2 = $read->{$key}->{r2}->{ref};
      $pos2 = $read->{$key}->{r2}->{pos};
#      $pos->{$ref2}->{$pos2}->{"0___unmapped"}->{count}++;
#      push @{$pos->{$ref2}->{$pos2}->{pair}}, "0___unmapped";
      print STDERR $pos1."___$ref2 -- 0___unmapped\n" if $print_dist;
      $single++;
    }
  }
}
undef %$read;

$median_dist = sv::array_median(@dist_total); # coz the distance takes direction both ways, but total count is only once
if ($ywin == 0) {
  $ywin = int($median_dist + (-$median_dist % 50))/2 if $median_dist > 0;
  $ywin = int(-$median_dist + ($median_dist % 50))/2 if $median_dist < 0;
}

# can set global distribution and counts now
# should collapse into bins instead of distances also
my $precount;	# keyed with ref, position, then y-binned distance
		# distance can also be "pair" - an array ref for where the
		#   other read of the pair maps
		# $precount is the data before it gets binned by coverage
foreach $ref (keys %$pos) {
  foreach $i (keys %{$pos->{$ref}}) {
    foreach $dist (keys %{$pos->{$ref}->{$i}}) {
      if ($dist eq "pair") {
        push @{$precount->{$ref}->{$i}->{pair}}, @{$pos->{$ref}->{$i}->{pair}};
        next;
      }
      if ($dist =~ /^(-?\d+)___(.*)/) {
        my $location = $1;
        my $name = $2;
        $j = int ($location/$ywin) * $ywin;
        $j = $j - $ywin if ($location < 0 || $location == -0);
        $precount->{$ref}->{$i}->{$j."___".$name} = $pos->{$ref}->{$i}->{$dist}->{count};
        $global_dist->{$j."___".$name} += $pos->{$ref}->{$i}->{$dist}->{count};
        $global_total += $pos->{$ref}->{$i}->{$dist}->{count};
      } else {
        if ($pos->{$ref}->{$i}->{$dist}->{count}) {
          $j = int($dist/$ywin) * $ywin;
          $j = $j - $ywin if ($dist < 0 || $dist == -0);
          $precount->{$ref}->{$i}->{$j} = $pos->{$ref}->{$i}->{$dist}->{count};
          $global_dist->{$j} += $pos->{$ref}->{$i}->{$dist}->{count};
          $global_total += $pos->{$ref}->{$i}->{$dist}->{count};
        }
      }
    }
  }
}
undef %$pos;

my $global_cdf = {};
foreach $bin (sort keysort keys %$global_dist) {
  $global_dist->{$bin} = $global_dist->{$bin} / $global_total;
  if ($print_dist) {
    print "$dist\t";
    print $global_dist->{$dist}, "\n";
  }
}
exit if $print_dist;

if (!$quiet) {
  print STDERR "Paired-end reads with one / both uniquely mapped: $single / $both\n";
  print STDERR "Median paired-end distance: $median_dist; Y-window bin size: $ywin\n";
#  print STDERR "Average read length: ", int($read_length_total/(2*$both)),"\n";
}
if ($range_cluster == 0) {
  $range_cluster = $median_dist;
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Reference length adjustment 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# adjust the length to be divisible by the y-window binning length
foreach $key (keys %$refh) {
  $i = $refh->{$key};

  # to make length fit y-win intervals
  $i = int(($i / $ywin) + 0.5) * $ywin;
  if ($i < $refh->{$key}){
    $i += $ywin;
  }
  $refh->{$key} = $i;
  if ($largest_genome < $i) {
    $largest_genome = $i;
  }
  $genome_size += $i;
}

printf STDERR ("Total genome size (rounded to bin size): %d (%d Y-bins)\n", $genome_size, 2*$genome_size/$ywin) if !$quiet;

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# THRESHOLDING FOR HIGHEST SNR
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# from fimx and strep dataset looks like optimal snr is 2x genome coverage
$ri = {};	# main relative entropy hash - keyed on reference,
		# bin (cov_bin), then one of rcount, binsize, entropy
$count = {};	# more detail - keyed on reference, bin (cov_bin), distance
		# distance can also be "pair" (array ref again)
my $max_ri = ();
my $max_count = ();
my $max_cov = 0;
my $snr = 0;
my $max_snr = 0;
if ($optimise_cov) {
  $cov_bin = 2*$read_length_total/$genome_size;
  $cov_bin = int(($cov_bin/100) + 0.5) * 100;
  if ($cov_bin == 0) {
    $cov_bin = 100;
  }
  
  print STDERR "Calculating coverage bin with highest SNR\n" if !$quiet;
  print STDERR "Cov\tSNR\n" if !$quiet;
  my $current = $cov_bin;
  my $inc = 10**(int(log($cov_bin)/log(10)));
  my $snh = {};
  ($ri, $count) = sv::ric($refh, $precount, $current, $global_dist);
  $snr = sv::snr($ri);
  if ($snr > $max_snr) {
    $max_snr = $snr;
    $max_ri = $ri;
    $max_count = $count;
    $max_cov = $current;
  }
  $snh->{$current} = $snr;
  printf STDERR "%d\t%.2f\n", $current, $snr if !$quiet;
  
  my $current_minus = 0;
  my $inc_minus = 10**(int(log($cov_bin)/log(10)) - 1) * 2;
  if ($current > $inc) {
    $current_minus = $current - $inc;
  } elsif ($current > $inc_minus) {
    $current_minus = $current - $inc_minus;
  }
  my $snr_minus = 0;
  if ($current_minus > 0) {
    ($ri, $count) = sv::ric($refh, $precount, $current_minus, $global_dist);
    $snr_minus = sv::snr($ri);
    if ($snr_minus > $max_snr) {
      $max_snr = $snr_minus;
      $max_ri = $ri;
      $max_count = $count;
      $max_cov = $current_minus;
    }
    $snh->{$current_minus} = $snr_minus;
    printf STDERR "%d\t%.2f\n", $current_minus, $snr_minus if !$quiet;
  }
  
  my $current_plus = $current + $inc;
  my $snr_plus = 0;
  ($ri, $count) = sv::ric($refh, $precount, $current_plus, $global_dist);
  $snr_plus = sv::snr($ri);
  $snh->{$current_plus} = $snr_plus;
  printf STDERR "%d\t%.2f\n", $current_plus, $snr_plus if !$quiet;
  
  if ($snr_plus > $snr) {
    $up = 0;  
    while ($current_plus <= $inc * 10 && $up < 2){
      $current_minus = $current;
      $snr_minus = $snr;
      $current = $current_plus;
      $snr = $snr_plus;
      $current_plus = $current + $inc;
      $snr_plus = 0;
      ($ri, $count) = sv::ric($refh, $precount, $current_plus, $global_dist);
      $snr_plus = sv::snr($ri);
      if ($snr_plus > $max_snr) {
        $max_snr = $snr_plus;
        $max_ri = $ri;
        $max_count = $count;
        $max_cov = $current_plus;
      }
      $snh->{$current_plus} = $snr_plus;
      printf STDERR "%d\t%.2f\n", $current_plus, $snr_plus if !$quiet;
	  
      if ($snr_plus && $snr_plus < $snr) {
        $up++;
      } else {
        $up = 0;
      }
    }
  }

  if ($snr_minus > $snr){
    while ($snr_minus > $snr && $current_minus > 0){
      $current_plus = $current;
      $snr_plus = $snr;
      $current = $current_minus;
      $snr = $snr_minus;
      if ($current > $inc) {
        $current_minus = $current - $inc;
      } elsif ($current > $inc_minus) {
        $current_minus = $current - $inc_minus;
      }
      $snr_minus = 0;
      if (!defined $snh->{$current_minus}) {
        ($ri, $count) = sv::ric($refh, $precount, $current_minus, $global_dist);
        $snr_minus = sv::snr($ri);
        if ($snr_minus > $max_snr) {
          $max_snr = $snr_minus;
          $max_ri = $ri;
          $max_count = $count;
          $max_cov = $current_minus;
        }
        $snh->{$current_minus} = $snr_minus;
        printf STDERR "%d\t%.2f\n", $current_minus, $snr_minus if !$quiet;
      }
    }
  }
  
  $ri = $max_ri;
  $count = $max_count;
  $cov_bin = $max_cov;
#  printf STDERR "Optimal coverage bin: %d (~%d coverage bins)\n", $cov_bin, ($both + $single)/$cov_bin if !$quiet;
  printf STDERR "Optimal coverage bin: %d (~%d coverage bins, average %d bp)\n", $cov_bin, $both/$cov_bin, $genome_size/$both*$cov_bin if !$quiet;
} else {
  ($ri, $count) = sv::ric($refh, $precount, $cov_bin, $global_dist);
  $snr = sv::snr($ri);
  printf STDERR "Coverage bin: %d (~%d coverage bins, average %d bp); SNR: %.2f\n", $cov_bin, $both/$cov_bin, $genome_size/$both*$cov_bin, $snr if !$quiet;
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Other file paths and variable check
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
$output_ic = $output."_re.txt";
$output_detail = $output."_detail.txt";
$output_list = $output."_SV.txt";
$output_png = $output."_graph.png";

#@@@@@@@@@@@@@@@@@@@@
# Process the counts
#@@@@@@@@@@@@@@@@@@@@
my $no_of_bins = 0;
foreach $ref (keys %$count) {
  $no_of_bins += scalar(keys %{$count->{$ref}});
}
my $no_of_ybins = round($genome_size*2 / $ywin); 

my @refsize = sort {$b<=>$a} values %$refh;
my @reforder = ();
my $refnames = "";
my $refsizes = "";
foreach my $i (@refsize) {
  push @reforder, grep { $refh->{$_} == $i } keys %$refh; 
}
my @refnames = keys %$count;
foreach my $name (@reforder) {
  push @refnames, grep { $count->{$_} eq $name } keys %$count; 
}

#@@@@@@@@@@@@@@@@@@@@@@
# P-value calculations
#@@@@@@@@@@@@@@@@@@@@@@
my $info_freq = {};
# P-value simulation by bootstrapping
if ($pval) {
  $info_freq = sv::bootstrap($global_dist, $bootstrap, $cov_bin, $quiet);
}

open I, ">$output_ic" || die "Cannot open file $output_ic: $!\n";
print I join ("\t", "# Chromosome", "Bin coord", "Bin size (reads)", "Bin size (bp)", "Relative entropy", "Adjusted P-value"), "\n";
open D, ">$output_detail" || die "Cannot open file $output_detail: $!\n";
print D "### $command_line\n";
print D "### Coverage bin: $cov_bin; y-window: $ywin\n";
print D join ("\t", "## Chromosome", "Bin coord", "Bin size (reads)", "Bin size (bp)"), "\n";
print D "## Mate pair location ranges are clustered if within $range_cluster bp\n";
print D join ("\t", "# Distance", "Read Pairs", "Expected", "Observed", "Relative Entropy"), "\n";
$pvalue = {};
my @pv = ();
my @ent = ();
my $corrected_pval;
my $significant = {};

# adjust the p-values by Benjamini-Yekutieli
if ($pval) {
  foreach $ref (@reforder) {
    foreach $bin (keys %{$ri->{$ref}}) {
      push @ent, $ri->{$ref}->{$bin}->{entropy};
    }
  }
  @ent = sort {$b <=> $a} @ent;
  $j = 0;
  $pos1 = 1/$bootstrap;
  while ($j <= $#ent && $ent[$j] > max(keys %$info_freq)) {
    $pvalue->{$ent[$j]} = 1/$bootstrap;
    push @pv, 1/$bootstrap;
    $j++;
  }
  $pos1 = $info_freq->{max(keys %$info_freq)};
  foreach $i (sort {$b <=> $a} keys %$info_freq) {
    while ($j <= $#ent && $ent[$j] > $i) {
      $pvalue->{$ent[$j]} = $pos1;
      push @pv, $pos1;
      $j++;
    }
    $pos1 += $info_freq->{$i};
    last if $j > $#ent;
  }
  while ($j <= $#ent) {
    $pvalue->{$ent[$j]} = 1;
    push @pv, 1;
    $j++;
  }
}
if ($pval) {
  @pv = sort {$a <=> $b} @pv;
  $corrected_pval = sv::pval_BY(\@pv);
  foreach $i (keys %$pvalue) {
    $pvalue->{$i} = $corrected_pval->{$pvalue->{$i}};
  }
}

print STDERR "Generating output:\n";
foreach $ref (@reforder) {
  print STDERR "  $ref ... " if !$quiet;
  foreach $bin (sort {$a <=> $b} keys %{$ri->{$ref}}) {
    print I join ("\t",
      $ref,
      $bin,
      $ri->{$ref}->{$bin}->{rcount},
      $ri->{$ref}->{$bin}->{binsize},
      $ri->{$ref}->{$bin}->{entropy});
    if ($pval) {
      print I "\t", $pvalue->{$ri->{$ref}->{$bin}->{entropy}}, "\n";
    } else {
      print I "\t0\n";
    }
    if ($pval && $pvalue->{$ri->{$ref}->{$bin}->{entropy}} < $fdr) {
      $significant->{$ref}->{$bin} = $ri->{$ref}->{$bin};
      print D join ("\t", "## $ref", $bin, $ri->{$ref}->{$bin}->{rcount}, $ri->{$ref}->{$bin}->{binsize}), "\n";
      @pair = ();
      my $out = {};
      foreach $dist (sort keysort keys %{$count->{$ref}->{$bin}}) {
        next if $dist eq "pair";
        next if $dist eq "rcount";
        next if $dist eq "pos";
        print D join ("\t",
          $dist,
          $count->{$ref}->{$bin}->{$dist},
          $global_dist->{$dist},
          $count->{$ref}->{$bin}->{$dist}/$ri->{$ref}->{$bin}->{rcount},
          $count->{$ref}->{$bin}->{$dist}/$ri->{$ref}->{$bin}->{rcount}*log($count->{$ref}->{$bin}->{$dist}/$ri->{$ref}->{$bin}->{rcount}/$global_dist->{$dist})), "\n";
      }
      if (defined $count->{$ref}->{$bin}->{pair} && ref($count->{$ref}->{$bin}->{pair}) eq "ARRAY"){
        @pair = @{$count->{$ref}->{$bin}->{pair}};
      }
      $pair = "";
      if (scalar(@pair) != 0) {
        $pair = sv::range(\@pair, $range_cluster); #collapse range within 50bp
        print D "# Pairs: $pair\n";
      } else {
        print D "# Pairs: none\n";
      }
    }
  }
  print STDERR "done\n" if !$quiet;
}
close I;

#@@@@@@@@@@@@@@@@@@
# Graphical output
#@@@@@@@@@@@@@@@@@@
$refnames = join("\",\"", @reforder);
$refnames = "\"".$refnames."\"";
$refsizes = join(",", @refsize);

#=begin GHOSTCODE
$r_file = $tempdir."/r";
open R, ">$r_file";
if ($pval) {
print R <<__END__;
ref <- c($refnames)
refsize <- c($refsizes)
x <- read.table(\"$output_ic\", sep="\\t", comment.char="", header=T)
np <- length(ref) * 2
s <- round(log10(refsize))
s[1] <- s[1] + 1.5
qvalue <- $fdr
ifelse(length(ref) > 1, linepos <- rep(c(3,2,2,3), length(ref)/2), linepos <- 2)
refpos <- rep(c(max(range(x\$Relative.entropy)), max(range(x\$Relative.entropy))-0.1), length(ref))
png(\"$output_png\", width=1000, height=700)
layout(matrix(1:np,2,np/2,byrow=TRUE),s); par(oma=c(3,2,3,4));
for (i in 1:length(ref)) {
  ifelse(i == 1, par(mar=c(0,2,4,0)), par(mar=c(0,0,4,0)))
  if (nrow(x[x\$Bin.coord >=0 & x\$X..Chromosome == ref[i],]) > 0) {
    plot(x[x\$Bin.coord >=0 & x\$X..Chromosome == ref[i], c(\"Bin.coord\",\"Relative.entropy\")], type="o", cex=0.5, xlim=c(0, refsize[i]), ylim=range(x\$Relative.entropy), xlab="", ylab="", xaxt="n", yaxt="n", axes=FALSE)
    points(x[x\$X..Chromosome == ref[i] & x\$Adjusted.P.value <= qvalue & x\$Bin.coord >= 0, c(\"Bin.coord\",\"Relative.entropy\")], cex=0.5, pch=8, col="red")
  } else {
    plot(0, 0, cex=0, xlim=c(0, refsize[i]), ylim=range(x\$Relative.entropy), xlab="", ylab="", xaxt="n", yaxt="n", axes=FALSE)
  }
  if (i %% 2 == 1) {
    mtext(ref[i], side=3, cex=0.75, line=linepos[i])
  } else {
    axis(3, xlim=c(0, refsize[i]), cex.axis=1.25)
  }
  ifelse(i == 1, axis(2, ylim=range(x\$Relative.entropy), cex.axis=1.25), ifelse(i == length(ref), axis(4, ylim=range(x\$Relative.entropy), cex.axis=1.25), NA))
} 
for (i in 1:length(ref)) {
  ifelse(i == 1, par(mar=c(4,2,0,0)), par(mar=c(4,0,0,0)))
  if (nrow(x[x\$Bin.coord < 0 & x\$X..Chromosome == ref[i],]) > 0) {
    plot(-x\$Bin.coord[x\$Bin.coord < 0 & x\$X..Chromosome == ref[i]], x\$Relative.entropy[x\$Bin.coord < 0 & x\$X..Chromosome == ref[i]], type="o", cex=0.5, xlim=c(0, refsize[i]), ylim=rev(range(x\$Relative.entropy)), xlab="", ylab="", xaxt="n", yaxt="n", axes=FALSE)
    points(-x\$Bin.coord[x\$X..Chromosome == ref[i] & x\$Adjusted.P.value <= qvalue & x\$Bin.coord < 0], x\$Relative.entropy[x\$X..Chromosome == ref[i] & x\$Adjusted.P.value <= qvalue & x\$Bin.coord < 0], cex=0.5, pch=8, col="red")
  } else {
    plot(0, 0, cex=0, xlim=c(0, refsize[i]), ylim=range(x\$Relative.entropy), xlab="", ylab="", xaxt="n", yaxt="n", axes=FALSE)
  }
  if (i %% 2 == 0) {
    mtext(ref[i], side=1, cex=0.75, line=linepos[i])
  } else {
    axis(1, xlim=c(0, refsize[i]), cex.axis=1.25)
  }
  ifelse(i == 1, axis(2, ylim=rev(range(x\$Relative.entropy)), cex.axis=1.25), ifelse(i == length(ref), axis(4, ylim=rev(range(x\$Relative.entropy)), cex.axis=1.25), NA))
} 
mtext(\"Forward strand position\", side=3, line=1, outer=TRUE)
mtext(\"Reverse strand position\", side=1, line=1, outer=TRUE)
mtext(paste(\"RIC score (Q-value = \",qvalue,\")\"), side=2, line=0, cex=1.5, outer=TRUE)
dev.off()
__END__
} else {
print R <<__END__;
ref <- c($refnames)
refsize <- c($refsizes)
x <- read.table(\"$output_ic\")
np <- length(ref) * 2
s <- round(log10(refsize))
s[1] <- s[1] + 1.5
ifelse(length(ref) > 1, linepos <- rep(c(3,2,2,3), length(ref)/2), linepos <- 2)
refpos <- rep(c(max(range(x\$Relative.entropy)), max(range(x\$Relative.entropy))-0.1), length(ref)) 
png(\"$output_png\", width=1000, height=700)
layout(matrix(1:np,2,np/2,byrow=TRUE),s); par(oma=c(3,2,3,1));
for (i in 1:length(ref)) {
  ifelse(i == 1, par(mar=c(0,2,4,0)), par(mar=c(0,0,4,0)))
  plot(x[x\$Bin.coord >=0 & x\$X..Chromosome == ref[i], c(\"Bin.coord\",\"Relative.entropy\")], type="o", cex=0.5, ylim=range(x\$Relative.entropy), xlab="", ylab="", xaxt="n", yaxt="n", axes=FALSE)
  if (i %% 2 == 1) {
    mtext(ref[i], side=3, cex=0.75, line=linepos[i])
  } else {
    axis(3, xlim=range(x\$Bin.coord[x\$X..Chromosome==ref[i] & x\$Bin.coord >=0]), cex.axis=1.25)
  }
  ifelse(i == 1, axis(2, ylim=range(x\$Relative.entropy), cex.axis=1.25), ifelse(i == length(ref), axis(4, ylim=range(x\$Relative.entropy), cex.axis=1.25), NA))
} 
for (i in 1:length(ref)) {
  ifelse(i == 1, par(mar=c(4,2,0,0)), par(mar=c(4,0,0,0)))
  plot(-x\$Bin.coord[x\$Bin.coord < 0 & x\$X..Chromosome == ref[i]], x\$Relative.entropy[x\$Bin.coord < 0 & x\$X..Chromosome == ref[i]], type="o", cex=0.5, ylim=rev(range(x\$Relative.entropy)), xlab="", ylab="", xaxt="n", yaxt="n", axes=FALSE)
  if (i %% 2 == 0) {
    mtext(ref[i], side=1, cex=0.75, line=linepos[i])
  } else {
    axis(1, xlim=range(x\$Bin.coord[x\$X..Chromosome==ref[i] & x\$Bin.coord >=0]), cex.axis=1.25)
  }
  ifelse(i == 1, axis(2, ylim=rev(range(x\$Relative.entropy)), cex.axis=1.25), ifelse(i == length(ref), axis(4, ylim=rev(range(x\$Relative.entropy)), cex.axis=1.25), NA))
} 
mtext(\"Forward strand position\", side=3, line=1, outer=TRUE)
mtext(\"Reverse strand position\", side=1, line=1, outer=TRUE)
mtext(\"RIC score\", side=2, line=0, cex=1.5, outer=TRUE)
dev.off()
__END__
}
close R;
#system "$R_command CMD BATCH --slave $r_file $tempdir/Rout";
system "$R_command $r_file";
#=end GHOSTCODE

#=cut

#@@@@@@@@@@@@@@@
# typing the SV
#@@@@@@@@@@@@@@@
# we can look at
# circularity at ends (should be at two of the four ends)
# deletions (inserts too large)
# duplication (distance near the size of the chromosome)
# insertions (lots of singles) - but this is confounded by being near repeats
# inversions (wrong mapping orientation)
# translocations (mapping to different chromosomes)
# false positives (too many spread out)
# widen all the ranges out by bin size at this point
if ($pval) {
  open S, ">$output_list";
  print S "## Individual SV calls\n";
  print S join ("\t", "# Chromosome", "Bin", "Entropy", "Adjusted P", "SV:target region(prob)"), "\n";
  my $merge_i = 0;
  my $merge_sv = {};	# keyed on number, then type, confidence, bp1, bp2
			# each bp has ref and coordinates expanded to ranges
  foreach $ref (keys %$significant) {
    foreach $bin (sort {$a <=> $b} keys %{$significant->{$ref}}) {
      my $sv = {};	# only for this (coverage) bin - keyed on y-binned
			# distance, then entropy, type, target
      my $b_entropy;
      my $b_sign;
      foreach $dist (keys %{$count->{$ref}->{$bin}}) {
        next if $dist eq "pair";
        next if $dist eq "rcount";
        next if $dist eq "pos";
        $b_entropy = $count->{$ref}->{$bin}->{$dist}/$ri->{$ref}->{$bin}->{rcount}*log($count->{$ref}->{$bin}->{$dist}/$ri->{$ref}->{$bin}->{rcount}/$global_dist->{$dist});
        next if $b_entropy <= 0;

        # translocation
        if (!sv::isfloat($dist)) {
          $sv->{$dist}->{entropy} = $b_entropy;
          $sv->{$dist}->{type} = "Translocation";
          $sv->{$dist}->{target} = $dist;	# should be format \d+___ref
          $sv->{$dist}->{target} =~ s/___\S+$//;
          $sv->{$dist}->{target} = sv::refadd($sv->{$dist}->{target}, -($ri->{$ref}->{$bin}->{binsize}), $refh->{$ref}) . ".." . sv::refadd($sv->{$dist}->{target}, $ri->{$ref}->{$bin}->{binsize}, $refh->{$ref});
          next;
        }

        # circularity
        # for FR - should be positive max coordinate, negative low coordinates
        # we should see positive distances near the chromosome length
        # for FF - should be all the extremes
        # we should see positive distances that are near normal if we are doing
        # only one distance per pair - so we shouldn't see it here
        if ($ori eq "FR" && ( ($bin > $refh->{$ref} - $median_dist) ||
                              ($bin < 0 && $bin > -$median_dist   ) )) {
          if ($dist > $refh->{$ref} - 1.5*$median_dist) {
            $sv->{$dist}->{entropy} = $b_entropy;
            $sv->{$dist}->{type} = "Circularity";
            $sv->{$dist}->{target} = $refh->{$ref} - $dist;
            $sv->{$dist}->{target} = sv::refadd($sv->{$dist}->{target}, -($ri->{$ref}->{$bin}->{binsize}), $refh->{$ref}) . ".." . sv::refadd($sv->{$dist}->{target}, $ri->{$ref}->{$bin}->{binsize}, $refh->{$ref});
            next;
          }
        }

        # deletion/duplication
        if ($dist > 1.5*$median_dist) {
          $sv->{$dist}->{entropy} = $b_entropy;
          $sv->{$dist}->{type} = "Deletion";
          $sv->{$dist}->{target} = $bin + abs($dist) - $median_dist;
          if ($dist > 0.5 * $refh->{$ref}) {
            $sv->{$dist}->{type} = "Duplication";
            $sv->{$dist}->{target} = $bin + ($refh->{$ref} - abs($dist)) - $median_dist;
          }
          if ($ori eq "FF") {
            $i = sv::unwind_distance($ref, $bin, $dist, $precount, $ri);
            if ($i eq "second") {
              $sv->{$dist}->{target} = $bin - abs($dist) + $median_dist;
            }
          }
          $sv->{$dist}->{target} = sv::refmod($sv->{$dist}->{target}, $refh->{$ref});
          $sv->{$dist}->{target} = sv::refadd($sv->{$dist}->{target}, -($ri->{$ref}->{$bin}->{binsize}), $refh->{$ref}) . ".." . sv::refadd($sv->{$dist}->{target}, $ri->{$ref}->{$bin}->{binsize}, $refh->{$ref});
          next;
        }

        # insertion
        if ($dist >= 0 && $dist < $median_dist) {
          $sv->{$dist}->{entropy} = $b_entropy;
          $sv->{$dist}->{type} = "Insertion";
          $sv->{$dist}->{target} = $bin + $dist;
          if ($sv->{$dist}->{target} > $refh->{$ref}) {
            $sv->{$dist}->{target} -= $refh->{$ref};
          } elsif ($sv->{$dist}->{target} > 0 && $bin < 0) {
            $sv->{$dist}->{target} -= $refh->{$ref};
          }
          $sv->{$dist}->{target} = sv::refadd($sv->{$dist}->{target}, -($ri->{$ref}->{$bin}->{binsize}), $refh->{$ref}) . ".." . sv::refadd($sv->{$dist}->{target}, $ri->{$ref}->{$bin}->{binsize}, $refh->{$ref});
          next;
        }

        # inversion
        if ($dist < 0) {
          $sv->{$dist}->{entropy} = $b_entropy;
          $sv->{$dist}->{type} = "Inversion";
          $sv->{$dist}->{target} = $bin + abs($dist);
          if ($ori eq "FR") {
            $i = sv::unwind_distance($ref, $bin, $dist, $precount, $ri);
            if ($i eq "second") {
              $sv->{$dist}->{target} = $bin - abs($dist);
            }
          }
          $sv->{$dist}->{target} = sv::refmod($sv->{$dist}->{target}, $refh->{$ref});
          $sv->{$dist}->{target} = sv::refadd($sv->{$dist}->{target}, -($ri->{$ref}->{$bin}->{binsize}), $refh->{$ref}) . ".." . sv::refadd($sv->{$dist}->{target}, $ri->{$ref}->{$bin}->{binsize}, $refh->{$ref});
          next;
        }
      }
      my $pos_e = 0;
      foreach $i (keys %$sv) {
        $pos_e += $sv->{$i}->{entropy};
      }
      my @out = ();
      @f = sort keysort keys %$sv;
      @r = ();
      $j = 0;
      $key = "";
      # targets at this point can only be ints or \d+___ref
      # actually change targets to ranges
      foreach $i (0..$#f) {
        if ($j == 0) {
          $j = $sv->{$f[$i]}->{entropy}/$pos_e;
          $key = $sv->{$f[$i]}->{type};
          @r = ($sv->{$f[$i]}->{target});
          next;
        }
        if ($sv->{$f[$i]}->{type} eq "Translocation" &&
            $key eq "Translocation" &&
            $sv->{$f[$i]}->{target} eq $sv->{$f[$i-1]}->{target}) {
          $j += $sv->{$f[$i]}->{entropy}/$pos_e;
          @r = ($sv->{$f[$i]}->{target});
        } elsif (sv::isfloat($f[$i]) && sv::isfloat($f[$i-1]) &&
                 $f[$i] - $f[$i-1] <= $ywin &&
                 $sv->{$f[$i]}->{type} eq $key) {
          $j += $sv->{$f[$i]}->{entropy}/$pos_e;
          push @r, $sv->{$f[$i]}->{target};
        } else {
          if ($j > $fdr) {
            $merge_i++;
            $merge_sv->{$merge_i}->{type} = $key;
            @{$merge_sv->{$merge_i}->{confidence}} = ($j);
            if ($key ne "Translocation") {
              push @out, sprintf("%s:%s(%.2f)", $key, sv::range(\@r, $range_cluster), $j);
              $merge_sv->{$merge_i}->{bp1}->{ref} = $ref;
              $merge_sv->{$merge_i}->{bp2}->{ref} = $ref;
              $merge_sv->{$merge_i}->{bp1}->{coord} = abs($bin) . ".." . (abs($bin) + $ri->{$ref}->{$bin}->{binsize} + 1);	# to make sure we overlap later with the next bin
              $merge_sv->{$merge_i}->{bp2}->{coord} = sv::absrange(\@r, $range_cluster);
            } else {
              push @out, sprintf("%s:%s(%.2f)", $key, $r[0], $j);
              $merge_sv->{$merge_i}->{bp1}->{ref} = $ref;
              $merge_sv->{$merge_i}->{bp2}->{ref} = $f[$i];	# this should be \d+___ref
              $merge_sv->{$merge_i}->{bp2}->{ref} =~ s/\d+___//;
              $merge_sv->{$merge_i}->{bp1}->{coord} = abs($bin) . ".." . (abs($bin) + $ri->{$ref}->{$bin}->{binsize} + 1);	# to make sure we overlap later with the next bin
              $merge_sv->{$merge_i}->{bp2}->{coord} = sv::absrange(\@r, $range_cluster);
            }
          }
          $j = 0;
        }
      }
      if ($j > $fdr) {
        $merge_i++;
        $merge_sv->{$merge_i}->{type} = $key;
        @{$merge_sv->{$merge_i}->{confidence}} = ($j);
        if ($key ne "Translocation") {	# not a translocation
          push @out, sprintf("%s:%s(%.2f)", $key, sv::range(\@r, $ri->{$ref}->{$bin}->{binsize}), $j);
          $merge_sv->{$merge_i}->{bp1}->{ref} = $ref;
          $merge_sv->{$merge_i}->{bp2}->{ref} = $ref;
          $merge_sv->{$merge_i}->{bp1}->{coord} = abs($bin) . ".." . (abs($bin) + $ri->{$ref}->{$bin}->{binsize} + 1);	# to make sure we overlap later with the next bin
          $merge_sv->{$merge_i}->{bp2}->{coord} = sv::absrange(\@r, $range_cluster);
        } else {		# a translocation
          push @out, sprintf("%s:%s(%.2f)", $key, $r[0], $j);
          $merge_sv->{$merge_i}->{bp1}->{ref} = $ref;
          $merge_sv->{$merge_i}->{bp2}->{ref} = $f[$#f];	# this should be \d+___ref
          $merge_sv->{$merge_i}->{bp2}->{ref} =~ s/\d+___//;
          $merge_sv->{$merge_i}->{bp1}->{coord} = abs($bin) . ".." . (abs($bin) + $ri->{$ref}->{$bin}->{binsize} + 1);	# to make sure we overlap later with the next bin
          $merge_sv->{$merge_i}->{bp2}->{coord} = sv::absrange(\@r, $range_cluster);
        }
      }
      if ($bin >= 0) {
        $b_sign = 1;
      } else {
        $b_sign = -1;
      }
      print S join ("\t", $ref,
                          "$bin.." . ($bin + $b_sign * $ri->{$ref}->{$bin}->{binsize}),
                          $ri->{$ref}->{$bin}->{entropy},
                          $pvalue->{$ri->{$ref}->{$bin}->{entropy}},
                          join (";", @out)), "\n";
    }
  }
  print S "## Combined SV calls\n";
  print S join ("\t", "# Breakpoint 1", "Breakpoint 2", "Type", "Confidence"), "\n";
  my ($k, $l);
  my $collapsed = {};
  my $merged;
  my $changed = 1;
  while ($changed) {
    $changed = 0;
    foreach $k (keys %$merge_sv) {
      next if !defined $merge_sv->{$k};	# we delete as we combine
      foreach $l (keys %$merge_sv) {
        next if $k == $l;
        next if !defined $merge_sv->{$l};
        $merged = 0;
        if ($merge_sv->{$k}->{type} eq $merge_sv->{$l}->{type}) {
          if ($merge_sv->{$k}->{type} eq "Translocation") {
            if ($merge_sv->{$k}->{bp1}->{ref} eq $merge_sv->{$l}->{bp1}->{ref} && $merge_sv->{$k}->{bp2}->{ref} eq $merge_sv->{$l}->{bp2}->{ref}) {
              if (sv::overlap($merge_sv->{$k}->{bp1}->{coord}, $merge_sv->{$l}->{bp1}->{coord}, $range_cluster) && sv::overlap($merge_sv->{$k}->{bp2}->{coord}, $merge_sv->{$l}->{bp2}->{coord}), $range_cluster) {
                $merged = 1;
                $merge_sv->{$k}->{bp1}->{coord} = sv::range([ $merge_sv->{$k}->{bp1}->{coord}, $merge_sv->{$l}->{bp1}->{coord} ], $range_cluster);
                $merge_sv->{$k}->{bp2}->{coord} = sv::range([ $merge_sv->{$k}->{bp2}->{coord}, $merge_sv->{$l}->{bp2}->{coord} ], $range_cluster);
              }
            } elsif ($merge_sv->{$k}->{bp1}->{ref} eq $merge_sv->{$l}->{bp2}->{ref} && $merge_sv->{$k}->{bp2}->{ref} eq $merge_sv->{$l}->{bp1}->{ref}) {
              if (sv::overlap($merge_sv->{$k}->{bp1}->{coord}, $merge_sv->{$l}->{bp2}->{coord}, $range_cluster) && sv::overlap($merge_sv->{$k}->{bp2}->{coord}, $merge_sv->{$l}->{bp1}->{coord}, $range_cluster)) {
                $merged = 1;
                $merge_sv->{$k}->{bp1}->{coord} = sv::range([ $merge_sv->{$k}->{bp1}->{coord}, $merge_sv->{$l}->{bp2}->{coord} ], $range_cluster);
                $merge_sv->{$k}->{bp2}->{coord} = sv::range([ $merge_sv->{$k}->{bp2}->{coord}, $merge_sv->{$l}->{bp1}->{coord} ], $range_cluster);
              }
            }
          } else {	# not translocation so only need to check first bp refs
            if ($merge_sv->{$k}->{bp1}->{ref} eq $merge_sv->{$l}->{bp1}->{ref}) {
              if (sv::overlap($merge_sv->{$k}->{bp1}->{coord}, $merge_sv->{$l}->{bp1}->{coord}, $range_cluster) && sv::overlap($merge_sv->{$k}->{bp2}->{coord}, $merge_sv->{$l}->{bp2}->{coord}, $range_cluster)) {
                $merged = 1;
                $merge_sv->{$k}->{bp1}->{coord} = sv::range([ $merge_sv->{$k}->{bp1}->{coord}, $merge_sv->{$l}->{bp1}->{coord} ], $range_cluster);
                $merge_sv->{$k}->{bp2}->{coord} = sv::range([ $merge_sv->{$k}->{bp2}->{coord}, $merge_sv->{$l}->{bp2}->{coord} ], $range_cluster);
              } elsif (sv::overlap($merge_sv->{$k}->{bp1}->{coord}, $merge_sv->{$l}->{bp2}->{coord}, $range_cluster) && sv::overlap($merge_sv->{$k}->{bp2}->{coord}, $merge_sv->{$l}->{bp1}->{coord}, $range_cluster)) {
                $merged = 1;
                $merge_sv->{$k}->{bp1}->{coord} = sv::range([ $merge_sv->{$k}->{bp1}->{coord}, $merge_sv->{$l}->{bp2}->{coord} ], $range_cluster);
                $merge_sv->{$k}->{bp2}->{coord} = sv::range([ $merge_sv->{$k}->{bp2}->{coord}, $merge_sv->{$l}->{bp1}->{coord} ], $range_cluster);
              }
            }
          }
          if ($merged) {
  #          $merge_sv->{$k}->{confidence} = 1 - (1-$merge_sv->{$k}->{confidence}) * (1-$merge_sv->{$l}->{confidence});
            push @{$merge_sv->{$k}->{confidence}}, @{$merge_sv->{$l}->{confidence}};
            delete($merge_sv->{$l});
            $changed = 1;
          }
        }
      }
      $i = join (",", $merge_sv->{$k}->{bp1}->{coord}, $merge_sv->{$k}->{bp2}->{coord});
      my @sort = split /,/, $i;
      $merge_sv->{$k}->{sort} = 0;
      foreach $i (@sort) {
        @f = split /\.\./, $i;
        foreach $j (@f) {
          $merge_sv->{$k}->{sort} = $j if $merge_sv->{$k}->{sort} == 0;
          $merge_sv->{$k}->{sort} = $j if $merge_sv->{$k}->{sort} > $j;
        }
      }
      # these are all supposed to be ranged and positive now
      $merge_sv->{$k}->{bp1}->{coord} =~ /(^\d+)/;
      $i = $1;
      $merge_sv->{$k}->{bp2}->{coord} =~ /(^\d+)/;
      $j = $1;
      if ($i > $j) {
        my $temp = $merge_sv->{$k}->{bp1};
        $merge_sv->{$k}->{bp1} = $merge_sv->{$k}->{bp2};
        $merge_sv->{$k}->{bp2} = $temp;
      }
    }
  }
  foreach $k (sort {$merge_sv->{$a}->{sort} <=> $merge_sv->{$b}->{sort}} keys %{$merge_sv}) {
    print S join ("\t", "$merge_sv->{$k}->{bp1}->{ref}:$merge_sv->{$k}->{bp1}->{coord}", "$merge_sv->{$k}->{bp2}->{ref}:$merge_sv->{$k}->{bp2}->{coord}", $merge_sv->{$k}->{type}, sv::geom_mean(@{$merge_sv->{$k}->{confidence}})), "\n";
  }
}
$now = time - $now;
printf D ("\n### Total running time: %02d:%02d:%02d\n", int($now/3600), int(($now % 3600)/60), int($now % 60));
close D;
close S;
exit;

sub keysort {
  if (sv::isfloat($a) && sv::isfloat($b)) {
    return $a <=> $b;
  } elsif ($a =~ /^-?\d+___\S+$/ && $b =~ /^-?\d+___\S+$/) {
    $a =~ /^(-?\d+)___(\S+$)/;
    my $p1 = $1;
    my $r1 = $2;
    $b =~ /^(-?\d+)___(\S+$)/;
    my $p2 = $1;
    my $r2 = $2;
    if ($r1 eq $r2) {
      return $p1 <=> $p2;
    } else {
      return $r1 cmp $r2;
    }
  } else {
    return $a cmp $b;
  }
}
