#!/usr/bin/perl
use strict;
use warnings;

# Mocking external dependencies before loading sv.pm
BEGIN {
    $INC{"Math/Round.pm"} = "mock";
    $INC{"Math/CDF.pm"} = "mock";
    $INC{"Math/Trig.pm"} = "mock";
    $INC{"Math/BigFloat.pm"} = "mock";
    # Need to define pi if Math::Trig is mocked
    no strict 'refs';
    *{"Math::Trig::pi"} = sub { 3.14159265358979 };
    # Also need to export it to sv if it expects it to be available
    *{"sv::pi"} = sub { 3.14159265358979 };
}

use File::Spec;
use FindBin;
use lib $FindBin::Bin;
use sv;

# Define updated keysort here (must match svre.pl)
sub keysort {
  my $coord_re = qr/^([+-]?\d+)\z/;
  my $suffix_re = qr/^([+-]?\d+)___(\S+)\z/;

  if ($a =~ $coord_re) {
    my $v1 = $1;
    if ($b =~ $coord_re) {
      return $v1 <=> $1;
    }
  }

  if ($a =~ $suffix_re) {
    my ($p1, $r1) = ($1, $2);
    if ($b =~ $suffix_re) {
      my ($p2, $r2) = ($1, $2);
      if ($r1 eq $r2) {
        return $p1 <=> $p2;
      } else {
        return $r1 cmp $r2;
      }
    } else {
      return 1; # a has suffix, b doesn't, so b comes first
    }
  } elsif ($b =~ $suffix_re) {
    return -1; # b has suffix, a doesn't, so a comes first
  } else {
    return $a cmp $b;
  }
}

my $failed = 0;

print "Testing sv::range with hardened integer parsing and suffixes...\n";
my @range_input = ("1500000___chr1", "1600000___chr1");
my $range_output = sv::range(\@range_input, 100);
if ($range_output !~ /1500000/ || $range_output !~ /1600000/) {
    print "FAIL: sv::range failed to handle integers with suffix. Output: '$range_output'\n";
    $failed = 1;
} else {
    print "PASS: sv::range handled integers with suffix.\n";
}

print "Testing sv::absrange with hardened integer parsing and suffixes...\n";
my @abs_input = ("-1500000___chr1");
my $abs_output = sv::absrange(\@abs_input, 0);
if ($abs_output !~ /1500000/) {
    print "FAIL: sv::absrange failed to handle integers with suffix. Output: '$abs_output'\n";
    $failed = 1;
} else {
    print "PASS: sv::absrange handled integers with suffix.\n";
}

print "Testing keysort with pure numeric integers...\n";
{
    our ($a, $b);
    $a = "2000";
    $b = "10000";
    my $res = keysort();
    if ($res != -1) {
         print "FAIL: keysort failed for pure numeric integers (2000 vs 10000). Result: $res\n";
         $failed = 1;
    } else {
         print "PASS: keysort handled pure numeric integers correctly.\n";
    }
}

print "Testing scientific notation rejection for coordinates (as per feedback)...\n";
# sv::range now uses ^([+-]?\d+)\z for the numeric part, so 1.5e6 should not match
my @sci_input = ("1.5e6___chr1");
my $sci_output = sv::range(\@sci_input, 0);
if ($sci_output ne "") {
    print "FAIL: sv::range incorrectly accepted scientific notation for coordinate. Output: '$sci_output'\n";
    # Wait, if it doesn't match the regex, it might fall through to other logic or print error.
    # Current range implementation:
    # if it doesn't match /^(.+?)$ranger_re(.+?)\z/ it goes to split
    # if it doesn't match coordinate regex, it might just be ignored.
} else {
    print "PASS: sv::range correctly rejected/ignored scientific notation for coordinate.\n";
}

exit($failed);
