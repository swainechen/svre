#!/usr/bin/perl

package sv; # using Perl v5.10.1

use 5; # requires atleast Perl 5 interpreter
use strict; # restricted - global vars referred by full package name
use warnings; # warnings for potential errors

require Exporter;
# allows to export names into caller's namespace

our(@ISA, @EXPORT, @EXPORT_OK);
# predeclares as global variables, to be used under 'strict'. Another way is to use vars, but is somewhat deprecated.

@ISA = qw(Exporter);
# if a method is not found in this package, Perl searches for it in the packages in the @ISA array

# Symbols, functions, variables, etc. to be exported by default
@EXPORT = qw();

# Symbols, functions, variables, etc. to be exported on request
@EXPORT_OK = qw();

# Other modules
use Math::Round;
use Math::CDF;
use Math::Trig;
use Math::BigFloat;
use List::Util qw(min max sum);
#--------------------------------------------------------------
# Sub-routines
# for arrays and hashes, pass them to the subroutine as references

# get co-ordinator ranges
# similar to slchen make_list, by increasing the range for +- 50 positions
# accept a range interval
sub range {
  my ($a, $r) = @_; # $a is a reference to an array, $r is range
  my $ranger = "..";
  my $ranger_re = $ranger;
  $ranger_re =~ s/\./\\./g;
  my $sep = ",";
  my @return;
  my $arrays;
  my $suffix;
  my ($i, $j, $k);
  my (@f, @g);
  my ($start, $end);
  my ($p0, $p1);
  my ($r0, $r1);
  my $inc = $r-1;
  $inc = 1 if $inc <= 0;

  foreach $i (@$a) {
    @f = split /,/, $i;
    foreach $j (0..$#f) {
      @g = split /$ranger_re/, $f[$j];
      if ($#g == 0) {
        if (isfloat($g[0])) {
          push @{$arrays->{__NUMBERS__}}, $g[0];
        } elsif ($g[0] =~ /^(-?\d+)___(\S+)$/) {
          push @{$arrays->{$2}}, $1;
        }
      } elsif ($#g == 1) {
        if (isfloat($g[0])) {	# assume $g[1] is a number to catch any errors
          if ($g[0] > $g[1]) {
            $k = $g[0];
            $g[0] = $g[1];
            $g[1] = $k;
          }
          push @{$arrays->{__NUMBERS__}}, $g[0];
          $k = 1;
          while ($g[0] + $k*($inc) < $g[1]) {
            push @{$arrays->{__NUMBERS__}}, $g[0] + $k*($inc);
            $k++;
          }
          push @{$arrays->{__NUMBERS__}}, $g[1];
        } else {	# again assume format to catch any errors
          $g[0] =~ /^(-?\d+)___(\S+)$/;
          $p0 = $1;
          $r0 = $2;
          $g[1] =~ /^(-?\d+)___(\S+)$/;
          $p1 = $1;
          $r1 = $2;
          print STDERR "Internal range format error: $f[$j]\n" if $r0 ne $r1;
          push @{$arrays->{$r0}}, $p0;
          if ($p0 > $p1) {
            $k = $p0;
            $p0 = $p1;
            $p1 = $k;
          }
          $k = 1;
          while ($p0 + $k*($inc) < $p1) {
            push @{$arrays->{$r0}}, $p0 + $k*($inc);
            $k++;
          }
          push @{$arrays->{$r0}}, $p1;
        }
      }
    }
  }

  foreach $i (keys %$arrays) {
    @{$arrays->{$i}} = sortu(@{$arrays->{$i}});
    undef $start;
    undef $end;
    if ($i eq "__NUMBERS__") {
      $suffix = "";
    } else {
      $suffix = "___".$i;
    }
    foreach $j (@{$arrays->{$i}}) {
      if (!defined $start) {
        $start = $j;
        $end = $j;
        next;
      }
      if ($j <= $end + $r) {
        $end = $j;
      } else {
        if ($end == $start) {
          push @return, $start.$suffix;
        } else {
          push @return, (join($ranger, $start, $end).$suffix);
        }
        $start = $j;
        $end = $j;
      }
    }
    if (defined $start && defined $end) {
      if ($end == $start) {
        push @return, $start.$suffix;
      } else {
        push @return, (join($ranger, $start, $end).$suffix);
      }
    }
  }
  return join($sep, @return);
}

sub absrange {	# same as above but return only absolute values
		# we will assume these are only numbers here
  my ($a, $r) = @_;
  my $b = [];
  my $i;
  my @f;
  foreach $i (@$a) {
    if (isfloat($i)) {
      push @$b, abs($i);
    } else {
      $i =~ /(-?\d+)\.\.(-?\d+)/;
      push @$b, abs($1) . ".." . abs($2);
    }
  }
  return (range($b, abs($r)));
}

sub overlap {
  my ($r, $s, $range) = @_;	# format 1,3,6..10,15 for both - can only take
				# numbers also we're only going to look at
				# absolute values
  $range = 0 if !defined $range || $range <= 0;
  my ($i, $j);
  my ($f1, $f2);
  my ($l1, $l2);
  my $t;
  my @r = split /,/, $r;
  my @s = split /,/, $s;
  foreach $i (@r) {
    if (isfloat($i)) {
      $f1 = $i;
      $l1 = $i;
    } else {
      ($f1, $l1) = split /\.\./, $i;
    }
    $f1 = abs($f1);
    $l1 = abs($l1);
    if ($f1 > $l1) {
      $t = $f1;
      $f1 = $l1;
      $l1 = $t;
    }
    foreach $j (@s) {
      if (isfloat($j)) {
        $f2 = $j;
        $l2 = $j;
      } else {
        ($f2, $l2) = split /\.\./, $j;
      }
      $f2 = abs($f2);
      $l2 = abs($l2);
      if ($f2 > $l2) {
        $t = $f2;
        $f2 = $l2;
        $l2 = $t;
      }
      if ( ( ($f1 >= $f2 - $range) && ($f1 <= $l2 + $range) ) ||
           ( ($l1 >= $f2 - $range) && ($l1 <= $l2 + $range) ) ||
           ( ($f2 >= $f1 - $range) && ($f2 <= $l1 + $range) ) ||
           ( ($l2 >= $f1 - $range) && ($l2 <= $l1 + $range) )    ) {
        return 1;
      }
    }
  }
  return 0;	# didn't get any match
}

# check if given number is within an array of numbers
sub pair_check {
  my @pairs = @{$_[0]}; # always pass the array as a reference to the sub when calling
  my $a_string = $_[1];  # joined (with ,) string of numbers or ranges 
  my $pos = $_[2];
  my @response = ();
  my $response = 0;

  $pos = $pos + 0.1*$pos;

  my @a = split(",", $a_string);
  foreach my $pair ( @pairs){
    foreach my $a (@a){
      if( $a !~ /\.\./ ){
        if( abs($pair) == abs($a)){
          $response = 1;
          last;
        } elsif( (abs($a) - $pos <= abs($pair) && abs($pair) <= abs($a) + $pos) || (abs($a) + $pos <= abs($pair) && abs($pair) <= abs($a) - $pos) ){
          $response = 1;
          last;
        }
      } else {
        my @as = split(/\.\./, $a);
        if( abs($pair) == abs($as[0]) || abs($pair) == abs($as[1]) ){
          $response = 1;
          last;
        } elsif( (abs($as[0]) - $pos <= abs($pair) && abs($pair) <= abs($as[1]) + $pos) || (abs($as[1]) - $pos <= abs($pair) && abs($pair) <= abs($as[0]) + $pos) ){
          $response = 1;
          last;
        }
      }
    }
    push @response, $response;
  }

  return @response;
}

# distributions

# Normal distribution
# Usage: normal(value, mean, sd, interval)
sub normal {
  my ($x, $mean, $sd, $interval) = @_;
  my $a = $x + $interval;
  my $prob_x = 0;
  my $prob_a = 0;

  $prob_x = Math::CDF::pnorm(($x-$mean)/$sd);
  $prob_a = Math::CDF::pnorm(($a-$mean)/$sd);

  return $prob_a - $prob_x;
}

# Poisson distribution
# Usage: poisson(rate, value)
sub poisson {
  my ($success_rate, $x) = @_;
  my $y = Math::BigFloat->new(round($x));
  my $z = Math::BigFloat->new($success_rate);

  my $numerator = $z->bpow($y);

  if( $numerator == 0 ){
    return 0;
  } else {
    my $fac = $y->bfac();
    my $e = exp($success_rate);
    $e = Math::BigFloat->new($e);
    my $denominator = $fac->bmul($e);
    my $result = $numerator->bdiv($denominator);
    return $result;
  }
}

# approximation of poisson probability
# return ln(poisson)
sub poissonApprox {
  my ($lambda, $k) = @_;
  my $a = $k * log($lambda);
  my $ln_poisson = $a - $lambda - factorialApprox($k);
  return $ln_poisson;
}

# Ramanujam ln(n!) factorial estimation
# For large values of n
# return ln(n!)
sub factorialApprox {
  my $x = shift @_;
  if( $x < 2 ){
    return 0;
  } else {
    my $a = $x * log($x) - $x;
    my $b = log($x*(1+4*$x*(1+2*$x)))/6;
    return $a+$b+log(pi)/2;
  }
}

# generate null distribution hashes

# null for y-bins
# Usage: null_distribution(genome_size, mean, sd, ywin)
# ywin can be 1 for convolution - per base distribution
sub null_distribution {
  my $genome_size = $_[0];
  my $mean = $_[1];
  my $sd = $_[2];
  my $ywin = $_[3];
  my $null = {};

  for( my $i = $ywin-$genome_size; $i <= $genome_size; $i+=$ywin ){
    $null->{$i} = normal($i, $mean, $sd, $ywin);
  } 

  return $null;
}

# null from file
# Usage: null_distribution_from_file(file, ywin)
# ywin can be 1 for convolution - per base distribution
sub null_distribution_from_file {
  my $file = $_[0];
  my $ywin = $_[1];
  my $null = {};
  my $total = 0;
  my $j = 0;

  open N, $file;
  while(<N>){
    chomp;
    $j = int($_/$ywin) * $ywin;
    $j = $j - $ywin if $_ < 0;
    $null->{$j} += 1;
    $total++;
  }

  foreach my $dist ( keys %$null ){
    $null->{$dist} = $null->{$dist} / $total;
  }

  return $null;
}

# rms variance
# Usage: rms_variance(info content, info content elements array, no. of bins)
sub rms_variance {
  my $info = $_[0];
  my @info_parts = $_[1];
  my $no_of_ybins = $_[2];
  my $var = 0;

  my $info_mean = $info / round($no_of_ybins);
  @info_parts = map{ ($_ - $info_mean) ** 2} @info_parts;
  my $sum = sum(0, @info_parts); 
  # for those bins that don't map
  $sum += ((0-$info_mean)**2) * (round($no_of_ybins) - scalar(@info_parts));

  $var = sqrt($sum/ round($no_of_ybins));
  return $var;
}

# signal-to-noise ratio
# usage: snr(ref to relative entropy hash)
# the input hash should be keyed on reference then bin
sub snr {
  my ($ri) = @_;

  my $min = 0; my $max = 0;
  my $lt = {}; my $mt = {};
  my $w1 = 0; my $w2 = 0;
  my $sd1 = 0; my $sd2 = 0;
  my $v = 0; my $vmin = 0;
  my $t = 0; my $snr = 0;
  my $noise = 0; my $noise_n = 0;
  my $signal = 0; my $signal_n = 0;
  my @ric = (); my $i;
  my $resolution = 0.01;

  foreach $w1 (keys %$ri) {	# ref
    foreach $w2 (keys %{$ri->{$w1}}) {	# bin
      push @ric, $ri->{$w1}->{$w2}->{entropy};
    }
  }
  @ric = sort {$a <=> $b} @ric;

  if ($resolution > ($ric[$#ric] - $ric[0])/1000) {
    $resolution = ($ric[$#ric] - $ric[0])/1000;
  }
  for ($i = $ric[0]; $i <= $ric[$#ric]; $i += $resolution) {
    ($min, $max) = sv::binary_search(\@ric, $i);
    last if $min == $#ric;
    push @{$lt->{$i}}, @ric[0..$min];
    push @{$mt->{$i}}, @ric[$min+1..$#ric];
    $w1 = scalar(@{$lt->{$i}})/scalar(@ric);
    $w2 = scalar(@{$mt->{$i}})/scalar(@ric);
    $sd1 = stdev(\@{$lt->{$i}});
    $sd2 = stdev(\@{$mt->{$i}});
    $v = ($w1*$sd1*$sd1) + ($w2*$sd2*$sd2);
    if ($v != 0) {
      if ($t == 0) {
        $t = $i;
        $vmin = $v;
      } elsif ($v < $vmin) {
        $t = $i;
        $vmin = $v;
      }
    }
  }

  ($min, $max) = sv::binary_search(\@ric, $t); 
  $noise_n = $min+1;
  $noise = sum(@ric[0..$min]);
  $signal_n = $#ric - $min;
  $signal = sum(@ric[$min+1..$#ric]);

  if ($noise * $signal_n) {
    $snr = ($signal * $noise_n)/($noise * $signal_n);
    return $snr;
  } else{
    return 0;
  }
}

# bootstrap 
# Usage: bootstrap(probability distribution hash, n, coverage to bootstrap for, whether to be quiet)
sub bootstrap {
  my ($prob_dist, $bootstrap, $cov, $quiet) = @_;
  my $info_freq = {};
  my $times = 1;
  my $cdf;
  my @cdf;
  my ($i, $c, $random);
  my $bins;
  my @random;
  my $progress;

#  foreach $i (sort {$a<=>$b} keys %$prob_dist) {
  foreach $i (sort {$prob_dist->{$b}<=>$prob_dist->{$a}} keys %$prob_dist) {
    $c += $prob_dist->{$i};
    $cdf->{$c} = $i;
    push @cdf, $c;
  }

  if (!$quiet) {
    print STDERR "Bootstrapping for P-value\n";
    print STDERR "[0.......50......100%]\n[";
  }

  while ($times <= $bootstrap) {
    $bins = {};
    @random = ();
    $random = 0;
    $i = 1;

    $progress = $times % int($bootstrap/20);
    print STDERR "*" if !$quiet && $progress == 0;
    while ($i <= $cov) {
      push @random, rand();
      $i++;
    }
    @random = sort(@random); 
    $random = shift @random;

    my ($min, $max) = (0, $#cdf);
    my ($start, $end) = (0, $#cdf); 

    ($min, $max) = binary_search(\@cdf, $random);
    $start = $min;

    ($min, $max) = binary_search(\@cdf, $random[$#random]);
    $end = $max;

    CDF: foreach $c (@cdf[$start..$end]) {
      while ($random <= $c) {
        $bins->{$cdf->{$c}} += 1;

        if (scalar(@random) != 0) {
          $random = shift @random;
        } else {
          $random = 2;
          last CDF;
        }
      }
    }

    my $info = 0;
    my $info_bin = 0;
    my $freq_bin = 0;
#    foreach $c (sort {$a<=>$b} keys %$bins) {}
    foreach $c (keys %$bins) {
      $freq_bin = $bins->{$c}/$cov;
      $info_bin = $freq_bin * log($freq_bin/$prob_dist->{$c}); 
      $info += $info_bin;
    }
    $info_freq->{$info} += 1;
    $times++;
  }
  print STDERR "]\n" if !$quiet;
  foreach $c (keys %$info_freq) {
    $info_freq->{$c} /= $bootstrap;
  }

  return $info_freq;
}

# binary search
# takes in an sorted array, a value to be searched
# returns the range of index to be searched
# adapted from the book "Mastering Algotrithms with Perl pg 3"
sub binary_search {
  my ($array, $value) = @_;
  my ($low, $high) = (0, $#$array);

  while ($low < $high) {
    my $try = int(($low+$high)/2);
    if ($array->[$try] == $value) {
      $low = $try;
      $high = $try;
      last;
    }
    if ($array->[$try] < $value) {
      $low = $try + 1; 
    }
    if ($array->[$try] > $value) {
      $high = $try - 1; 
    }
  }

  if ($low <= 0) {
    $low = 0;
  }
  if ($high <= 0) {
    $high = 0;
  }
  if ($low >= $#$array) {
    $low = $#$array;
  }
  if ($high >= $#$array) {
    $high = $#$array;
  }

  if ($low < $high) {
    return ($low, $high);
  } else {
    return ($high, $low);
  }
}

# creating structural variations
# ones marked adapted are from the book "Beginers Perl for Bioinfomatics"

# A subroutine to randomly select a position in a string.
# WARNING: make sure you call srand to seed the
# random number generator before you call this function.
# adapted*
sub randomposition {

    my $length = shift @_;

    # Notice the "nested" arguments:
    # $string is the argument to length
    # $length is the argument to rand
    # rand($length) is the argument to int
    # int(rand($length)) is the argument to return
    # But we write it without parentheses, as permitted.
    #
    # rand returns a decimal number between 0 and its argument.
    # int returns the integer portion of a decimal number.
    #
    # The whole expression returns a random number between 0 and length-1,
    # which is how the positions in a string are numbered in Perl.

    return int rand $length;
}

# A subroutine to randomly generate a size
sub randomsize {
  my $length = shift @_; 
  my $read_length = shift @_;
  my $random = 0;

  $length = $length * 0.01;
  while ($random < $read_length * 2) {
    $random = int(rand($length));
  } 
  return $random;
}

# A subroutine to randomly select an element from an array
# WARNING: make sure you call srand to seed the
# random number generator before you call this function.
# adapted *
sub randomelement {
    my(@array) = @_;
    return $array[rand @array];
}


# A subroutine to select at random one of the four nucleotides
# WARNING: make sure you call srand to seed the
# random number generator before you call this function.
# adapted *
sub randomnucleotide {
    my(@nucleotides) = ('A', 'C', 'G', 'T');
    return randomelement(@nucleotides);
}

# a subroutine to generate a random sequence of the given size
sub randomsequence {
  my $size = shift @_;
  my $i = 0;
  my $seq = "";
  while ($i < $size) {
    $seq .= randomnucleotide;
    $i++;
  }
  return $seq;
}

# a subroutine to create an insertion in a given sequence
sub insertion {
  my $parent = shift @_;
  my $read_length = shift @_;
  my $position = randomposition(length($parent));
  my $size = randomsize(length($parent), $read_length);
  my $insert = randomsequence($size);

  substr($parent, $position, 0) = $insert;

  print STDERR "Insertion\n";
  print STDERR "At position $position, of size $size\n";
  print STDERR "Inserted sequence: $insert\n";

  return ($parent, $position, $size, $insert);
}

# a subroutine to create a deletion in a given sequence
sub deletion {
  my $parent = shift @_;
  my $read_length = shift @_;
  my $position = randomposition(length($parent));
  my $size = randomsize(length($parent), $read_length);

  while ( $size > length($parent) - $position ){
    $position = randomposition(length($parent));
  }

  my $deleted = substr($parent, $position, $size);
  substr($parent, $position, $size) = "";

  print STDERR "Deletion\n";
  print STDERR "At position $position, of size $size\n";
  print STDERR "Deleted sequence: $deleted\n";

  return ($parent, $position, $size, $deleted);
}

# a subroutine to create a inversion in a given sequence 
sub inversion {
  my $parent = shift @_;
  my $read_length = shift @_;
  my $position = randomposition(length($parent));
  my $size = randomsize(length($parent), $read_length);

  while ( $size > length($parent) - $position ){
    $position = randomposition(length($parent));
  }

  my $inverted = substr($parent, $position, $size);
  $inverted = revcomp($inverted);
  substr($parent, $position, $size) = "$inverted";

  print STDERR "Inversion\n";
  print STDERR "At position $position, of size $size\n";
  print STDERR "Inverted sequence: $inverted\n";

  return ($parent, $position, $size, $inverted);
}

# a subroutine to create a tandem duplication
sub tandem {
  my $parent = shift @_;
  my $read_length = shift @_;
  my $position = randomposition(length($parent));
  my $size = randomsize(length($parent), $read_length);
  
  while ( $size > length($parent) - $position ){
    $position = randomposition(length($parent));
  }

  print STDERR "Tandem duplication\n";
  print STDERR "Original sequence from position $position\n";

  my $tandem = substr($parent, $position, $size);
  $position = $position + $size;
  substr($parent, $position, 0) = "$tandem";

  print STDERR "At position $position, of size $size\n";
  print STDERR "Duplicated sequence: $tandem\n";

  return ($parent, $position-$size, $position, $size, $tandem);
}


# a subroutine to create a duplication (non-tandem)
sub duplication {
  my $parent = shift @_;
  my $read_length = shift @_;
  my $position = randomposition(length($parent));
  my $size = randomsize(length($parent), $read_length);

  while ( $size > length($parent) - $position ){
    $position = randomposition(length($parent));
  }

  my $duplicate = substr($parent, $position, $size);
  my $position2 = $position;
  while( $position2 == $position ){
    $position2 = randomposition(length($parent));
  }
  substr($parent, $position2, 0) = "$duplicate";

  print STDERR "Non-Tandem duplication\n";
  print STDERR "Original sequence from position $position\n";
  print STDERR "At position $position2, of size $size\n";
  print STDERR "Duplicated sequence: $duplicate\n";

  return ($parent, $position, $position2, $size, $duplicate);
}

# a subroutine to create all sv in a given sequence
sub createsv{
  my $parent = shift @_;
  my $read_length = shift @_;
  my $ref = {};
  my $seq = "";
  my $pos = 0; 
  my $pos2 = 0;
  my $size = 0;
  my $mutant = "";
  my $position = {};
  my $check = 0;
  my $p = 0;

  ($seq, $pos, $pos2, $size, $mutant) = duplication($parent, $read_length);
  $ref->{duplication}->{mutant} = $mutant;
  $ref->{duplication}->{size} = $size;
  $ref->{duplication}->{pos} = $pos2; # this is so we can sort acc to pos later
  $ref->{duplication}->{pos2} = $pos;
  $position->{$pos} = $pos + $size;
  $position->{$pos2} = $pos2 + $size;

  CHECK: while( $check == 0 ){
    ($seq, $pos, $size, $mutant) = deletion($parent, $read_length);
    foreach $p (keys %$position){
      if( $pos >= $p && $pos <= $position->{$p} ){
        redo CHECK;
      }
    }
    $check = 1;
  }
  $ref->{deletion}->{mutant} = $mutant;
  $ref->{deletion}->{size} = $size;
  $ref->{deletion}->{pos} = $pos;
  $position->{$pos} = $pos + $size;

  $check = 0;
  CHECKa: while( $check == 0 ){
    ($seq, $pos, $size, $mutant) = inversion($parent, $read_length);
    foreach $p (keys %$position){
      if( $pos >= $p && $pos <= $position->{$p} ){
        redo CHECKa;
      }
    }
    $check = 1;
  }
  $ref->{inversion}->{mutant} = $mutant;
  $ref->{inversion}->{size} = $size;
  $ref->{inversion}->{pos} = $pos;
  $position->{$pos} = $pos + $size;

  $check = 0;
  CHECKb: while( $check == 0 ){
    ($seq, $pos, $pos2, $size, $mutant) = tandem($parent, $read_length);
    foreach $p (keys %$position){
      if( $pos >= $p && $pos <= $position->{$p} ){
        redo CHECKb;
      }
      if( $pos2 >= $p && $pos2 <= $position->{$p} ){
        redo CHECKb;
      }
    }
    $check = 1;
  }
  $ref->{tandem}->{mutant} = $mutant;
  $ref->{tandem}->{size} = $size;
  $ref->{tandem}->{pos} = $pos2; #this is so that we can sort by positions later
  $ref->{tandem}->{pos2} = $pos;
  $position->{$pos} = $pos + $size;

  $check = 0;
  CHECKc: while( $check == 0 ){
    ($seq, $pos, $size, $mutant) = insertion($parent, $read_length);
    foreach $p (keys %$position){
      if( $pos >= $p && $pos <= $position->{$p} ){
        redo CHECKc;
      }
    }
    $check = 1;
  }
  $ref->{insertion}->{mutant} = $mutant;
  $ref->{insertion}->{size} = $size;
  $ref->{insertion}->{pos} = $pos;
  $position->{$pos} = $pos + $size;

#print "done all work\n";
#print "now implementing\n";

  foreach my $sv (sort {$ref->{$b}->{pos} <=> $ref->{$a}->{pos}} keys %$ref ){
#print "$sv\t";
#print $ref->{$sv}->{pos}, "\n";
    if( $sv eq "duplication" ){
      substr($parent, $ref->{$sv}->{pos}, 0) = $ref->{$sv}->{mutant}; 
      print STDERR "Non-Tandem duplication\n";
      print STDERR "Original sequence from position ",$ref->{$sv}->{pos2},"\n";
      print STDERR "At position ", $ref->{$sv}->{pos}," of size ", $ref->{$sv}->{size}, "\n";
      print STDERR "Duplicated sequence: ", $ref->{$sv}->{mutant}, "\n";
    } elsif( $sv eq "tandem" ){
      substr($parent, $ref->{$sv}->{pos}, 0) = $ref->{$sv}->{mutant};
      print STDERR "Tandem duplication\n";
      print STDERR "Original sequence from position ",$ref->{$sv}->{pos2},"\n";
      print STDERR "At position ", $ref->{$sv}->{pos}," of size ", $ref->{$sv}->{size}, "\n";
      print STDERR "Duplicated sequence: ", $ref->{$sv}->{mutant}, "\n";
    } elsif( $sv eq "insertion" ){
      substr($parent, $ref->{$sv}->{pos}, 0) = $ref->{$sv}->{mutant};
      print STDERR "Insertion\n";
      print STDERR "At position ", $ref->{$sv}->{pos}," of size ", $ref->{$sv}->{size}, "\n";
      print STDERR "Inserted sequence: ", $ref->{$sv}->{mutant}, "\n";
    } elsif( $sv eq "inversion" ){
      substr($parent, $ref->{$sv}->{pos}, $ref->{$sv}->{size}) = $ref->{$sv}->{mutant};
      print STDERR "Inversion\n";
      print STDERR "At position ", $ref->{$sv}->{pos}," of size ", $ref->{$sv}->{size}, "\n";
      print STDERR "Inverted sequence: ", $ref->{$sv}->{mutant}, "\n";
    } elsif( $sv eq "deletion" ){
      substr($parent, $ref->{$sv}->{pos}, $ref->{$sv}->{size}) = "";
      print STDERR "Deletion\n";
      print STDERR "At position ", $ref->{$sv}->{pos}," of size ", $ref->{$sv}->{size}, "\n";
      print STDERR "Deleted sequence: ", $ref->{$sv}->{mutant}, "\n";
    } 
  }

  return $parent;
}

# average and stdev from http://edwards.sdsu.edu/research/index.php/kate/302-calculating-the-average-and-standard-deviation
sub average {
  my($data) = @_;
  if (not @$data) {
    die("Empty array\n");
  }
  my $total = 0;
  foreach (@$data) {
    $total += $_;
  }
  my $average = $total / @$data;
  return $average;
}

sub stdev {
  my ($data) = @_;
  if (@$data == 1) {
    return 0;
  }
  my $average = &average($data);
  my $sqtotal = 0;
  foreach(@$data) {
    $sqtotal += ($average-$_) ** 2;
  }
  my $std = ($sqtotal / (@$data-1)) ** 0.5;
  return $std;
}

sub array_median {
  # again, throwing out non number values
  my @a = @_;
  my ($t, $i, $j, @b);
  foreach $t (@a) {
    if (isfloat($t)) {
      push @b, $t;
    }
  }
  @b = sort { $a <=> $b } @b;
  return $b[0] if scalar @b == 1;
  return undef if !scalar @b;
  $i = int $#b/2;
  $j = int ($#b+1)/2;
  return ($b[$i] + $b[$j])/2;
}

sub revcomp {
  my ($inline) = $_[0];
  my $outline = reverse ($inline);
  $outline =~ tr/ABCDGHKMNRSTVWXYabcdghkmnrstvwxy/TVGHCDMKNYSABWXRtvghcdmknysabwxr/;
  return $outline;
}

# sort unique, be smart about numbers and non numbers
sub sortu {
  my @a = @_;
  my @defined = ();
  my @undefined = ();
  my $number = 1;
  foreach my $i (0..$#a) {
    if (!defined $a[$i]) {
      push @undefined, $i;
    }
    push @defined, $i;
    if (!isfloat($a[$i])) {
      $number = 0;
    }
  }
  if ($number) {
    @a = sort { $a <=> $b } @a[@defined];
  } else {
    @a = sort { $a cmp $b } @a[@defined];
  }
  foreach my $i (reverse 1..$#a) {
    if ($number) {
      if ($a[$i] == $a[$i-1]) { splice @a, $i, 1; }
    } else {
      if ($a[$i] eq $a[$i-1]) { splice @a, $i, 1; }
    }
  }
  foreach my $i (@undefined) {
    push @a, undef;
  }
  return @a;
}

sub isfloat {
  my ($string) = @_;
  if (!defined $string) {
    return 0;
  }
  chomp $string;
  if ($string =~ /^([+-])?(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/) {
    return 1;
  } else {
    return 0;
  }
}

sub ric {
  my ($refh, $precount, $cov_bin, $global_dist) = @_;
  my $ri = {};
  my $count = {};
  my $rcount;
  my $rcount_p;
  my $rcount_n;
  my $pos_count_p;
  my $pos_count_n;
  my ($p1, $p2);
  my (@pp, @pn);
  my ($i, $j);
  my ($info, $exp_prob, $prob);
  my $dist;
  my @sortbin;
  my $bin_size;

  foreach my $ref (keys %$refh) {
    $rcount_p = 0;
    $rcount_n = 0;
    $pos_count_p = 0;
    $pos_count_n = 0;
    $p1 = 0;
    $p2 = 0;
    @pp = ();
    @pn = ();
    
    foreach $i (1..$refh->{$ref}) {
      $pos_count_p++;
      $pos_count_n++;
      if (defined $precount->{$ref}->{$i}) {
        if ($rcount_p >= $cov_bin) {
          $rcount_p = 0;
        }
        if ($rcount_p == 0) {
          $p1 = $i;
          $count->{$ref}->{$p1}->{pos} = 0;
          $pos_count_p = 0;
        }
        if (defined $precount->{$ref}->{$i}) {
          if (defined $precount->{$ref}->{$i}->{pair} && ref($precount->{$ref}->{$i}->{pair}) eq "ARRAY"){
            push @{$count->{$ref}->{$p1}->{pair}}, @{$precount->{$ref}->{$i}->{pair}};
          }
          foreach $j (keys %{$precount->{$ref}->{$i}}){
            next if $j eq "pair";
            $count->{$ref}->{$p1}->{$j} += $precount->{$ref}->{$i}->{$j};
            $rcount_p += $precount->{$ref}->{$i}->{$j};
            $count->{$ref}->{$p1}->{rcount} += $precount->{$ref}->{$i}->{$j};
          }
        }
        $count->{$ref}->{$p1}->{pos} = $pos_count_p;
      }

      if (defined $precount->{$ref}->{-1*$i}) {
        if ($rcount_n >= $cov_bin) {
          $rcount_n = 0;
        }
        if ($rcount_n == 0) {
          $p2 = -1 * $i;
          $count->{$ref}->{$p2}->{pos} = 0;
          $pos_count_n = 0;
        }
        if (defined $precount->{$ref}->{-1*$i}) {
          if (defined $precount->{$ref}->{-1*$i}->{pair} && ref($precount->{$ref}->{-1*$i}->{pair}) eq "ARRAY") {
            push @{$count->{$ref}->{$p2}->{pair}}, @{$precount->{$ref}->{-1*$i}->{pair}};
          }
          foreach $j (keys %{$precount->{$ref}->{-1*$i}}){
            next if $j eq "pair";
            $count->{$ref}->{$p2}->{$j} += $precount->{$ref}->{-1*$i}->{$j};
            $rcount_n += $precount->{$ref}->{-1*$i}->{$j};
            $count->{$ref}->{$p2}->{rcount} += $precount->{$ref}->{-1*$i}->{$j};
          }
        }
        $count->{$ref}->{$p2}->{pos} = $pos_count_n;
      }
    }

    @sortbin = sort {$a<=>$b} keys %{$count->{$ref}};
    BIN: foreach my $bin (@sortbin) {
      next if !defined $count->{$ref}->{$bin};	# in case we deleted the last one
      next if !defined $count->{$ref}->{$bin}->{rcount};
      $rcount = $count->{$ref}->{$bin}->{rcount};
      $bin_size = $count->{$ref}->{$bin}->{pos};

      if ($bin == $sortbin[0] && defined($sortbin[1]) && $rcount < $cov_bin/2) {
        $count->{$ref}->{$sortbin[1]}->{rcount} += $rcount;
        $count->{$ref}->{$sortbin[1]}->{pos} += $bin_size;
        push @{$count->{$ref}->{$sortbin[1]}->{pair}}, @{$count->{$ref}->{$bin}->{pair}};
        foreach $dist (keys %{$count->{$ref}->{$bin}}) {
          next if $dist eq "pair";
          next if $dist eq "rcount";
          next if $dist eq "pos";
          $count->{$ref}->{$sortbin[1]}->{$dist} += $count->{$ref}->{$bin}->{$dist};
        }
        delete($count->{$ref}->{$sortbin[0]});
        next BIN;
      }
	  
      if (defined($sortbin[$#sortbin - 1]) && $bin == $sortbin[$#sortbin - 1] && $count->{$ref}->{$sortbin[$#sortbin]}->{rcount} < $cov_bin/2) {
        $count->{$ref}->{$bin}->{rcount} += $count->{$ref}->{$sortbin[$#sortbin]}->{rcount};
        $count->{$ref}->{$bin}->{pos} += $count->{$ref}->{$sortbin[$#sortbin]}->{pos};
        push @{$count->{$ref}->{$bin}->{pair}}, @{$count->{$ref}->{$sortbin[$#sortbin]}->{pair}};
        foreach $dist (keys %{$count->{$ref}->{$sortbin[$#sortbin]}}) {
          next if $dist eq "pair";
          next if $dist eq "rcount";
          next if $dist eq "pos";
          $count->{$ref}->{$bin}->{$dist} += $count->{$ref}->{$sortbin[$#sortbin]}->{$dist};
        }
        $rcount = $count->{$ref}->{$bin}->{rcount};
        $bin_size = $count->{$ref}->{$bin}->{pos};
        delete($count->{$ref}->{$sortbin[$#sortbin]});
      }

      $info = 0;
      $exp_prob = 0;
      $prob = 0;
#      foreach my $dist (sort{$a<=>$b} keys %{$count->{$ref}->{$bin}}){ }
      foreach $dist (keys %{$count->{$ref}->{$bin}}) { 
        next if $dist eq "pair";
        next if $dist eq "rcount";
        next if $dist eq "pos";
        my $info_w = 0;
        $prob = 0;
        $prob = $count->{$ref}->{$bin}->{$dist}/$rcount;
        $exp_prob = $global_dist->{$dist};

        $info_w = $prob*(log($prob/$exp_prob));
        $info = $info + $info_w;
      }
      $ri->{$ref}->{$bin}->{rcount} = $rcount;
      $ri->{$ref}->{$bin}->{binsize} = $bin_size;
      $ri->{$ref}->{$bin}->{entropy} = $info;
    }
  }
  return ($ri, $count);
}

sub pval_BY {
  # translated from R p.adjust
  my ($p) = @_;	# $p should be a ref to an array, sorted in increasing order
  my $return = {};
  my $n = scalar(@$p);
  my $i;
  my @sofar = ();
  my $q = 0;
  my $min = 1;
  my $new;

  foreach $i (1..$n) {
    $q += 1/$i;
  }
  foreach $i (reverse (0..$#$p)) {
    $new = $q * $n/($i+1) * $p->[$i];
    if ($new < $min) {
      $min = $new;
    }
    $return->{$p->[$i]} = $min;
  }
  return $return;
}

sub unwind_distance {
  my ($ref, $bin, $dist, $precount, $ri) = @_;
  # for FR, negative distances are unclear if we use only one distance per pair
  # for FF, it's positive distances
  my $i;
  my $j;
  my $r = { "first" => 0, "second" => 0 };
  # first means:
  #   i.e. -->.....--> $bin contains the read on the left - positive coords
  #   i.e. <--.....<-- $bin contains the read on the right - negative coords
  #   here the pair position is greater than the read we're looking at
  # second means you should subtract the distance to the first read
  #   i.e. -->.....--> $bin contains the read on the right
  #   i.e. <--.....<-- $bin contains the read on the left
  #   here the pair position is less than the read we're looking at
  # both means we have pairs on both sides
  # we should be able to work this out, this probably happens if we have
  # multiple reads mapping to the same position, and their pairs map on both
  # sides of where we are - but hopefully this doesn't happen too much given
  # the amount of coverage that would be needed
  foreach $i ($bin .. $bin + $ri->{$ref}->{$bin}->{binsize}) {
    next if !defined $precount->{$ref}->{$i}->{$dist};
    next if !defined $precount->{$ref}->{$i}->{pair};	# should never happen
    foreach $j (@{$precount->{$ref}->{$i}->{pair}}) {
      next if $j !~ /^[+-]?\d+$/;
      if ($j > $i) {
        $r->{first} = 1;
      } else {
        $r->{second} = 1;
      }
    }
  }
  if ($r->{first} && $r->{second}) {
    return("both");
  }
  if ($r->{first}) {
    return("first");
  }
  if ($r->{second}) {
    return("second");
  }
  return("error - no pairs");
}

sub geom_mean {
  my (@a) = @_;
  my $a;
  my $count = 0;
  my $total = 1;
  foreach $a (@a) {
    if (isfloat($a) && $a >= 0) {
      $total *= $a;
      $count++;
    } else {
      print STDERR "sub geom_mean input $a is not a nonnegative number\n";
    }
  }
  if ($total == 0 || $count == 0) {
    return 0;
  } else {
    return ($total ** (1/$count));
  }
}

sub refmod {
  # put the coordinate back in a reasonable range for the reference
  my ($coord, $ref) = @_;
  if ($coord == 0) {
    return $ref;
  } elsif ($coord > 0) {
    while ($coord > $ref) {
      $coord -= $ref;
    }
    return $coord;
  } else {
    while ($coord < -$ref) {
      $coord += $ref;
    }
    return $coord;
  }
}

sub refadd {
  # do addition of distances, taking into account negative coordinates and
  # circularity - $ref is total length, has to be positive
  my ($a, $b, $ref) = @_;
  my $return;
  return $a if $b == 0;
  $ref = 10000000 if $ref <= 0;
  $a = $ref if $a == 0;
  $return = $a + $b;
  while ($return > $ref) {
    $return -= $ref;
  }
  while ($return <= -$ref) {
    $return += $ref;
  }
  if ($return == 0) {
    return $ref;
  } else {
    return $return;
  }
}

#-------------------------------------------------------------
# Perl uses 'eval' to process the code. Thus for 'eval' to
# evaluate to TRUE and not fail, we provide a 1
1;

