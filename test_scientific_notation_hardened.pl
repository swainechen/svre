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

# Define updated keysort here
sub keysort {
  my $num_re = qr/([+-]?(?=\d|\.\d)\d*(?:\.\d*)?(?:[Ee][+-]?\d+)?)/;
  if ($a =~ /^$num_re\z/) {
    my $v1 = $1;
    if ($b =~ /^$num_re\z/) {
      return $v1 <=> $1;
    }
  }

  if ($a =~ /^${num_re}___(\S+)\z/) {
    my ($p1, $r1) = ($1 + 0, $2);
    if ($b =~ /^${num_re}___(\S+)\z/) {
      my ($p2, $r2) = ($1 + 0, $2);
      if ($r1 eq $r2) {
        return $p1 <=> $p2;
      } else {
        return $r1 cmp $r2;
      }
    } else {
      return 1; # a has suffix, b doesn't, so b comes first
    }
  } elsif ($b =~ /^${num_re}___(\S+)\z/) {
    return -1; # b has suffix, a doesn't, so a comes first
  } else {
    return $a cmp $b;
  }
}

my $failed = 0;

print "Testing sv::range with scientific notation and suffixes...\n";
my @range_input = ("1.5e6___chr1", "1.6e6___chr1");
my $range_output = sv::range(\@range_input, 100);
# Expected: "1500000..1600000___chr1" (or similar normalized)
if ($range_output !~ /1500000/ || $range_output !~ /1600000/) {
    print "FAIL: sv::range failed to handle scientific notation with suffix. Output: '$range_output'\n";
    $failed = 1;
} else {
    print "PASS: sv::range handled scientific notation with suffix.\n";
}

print "Testing sv::absrange with scientific notation and suffixes...\n";
my @abs_input = ("-1.5e6___chr1");
my $abs_output = sv::absrange(\@abs_input, 0);
if ($abs_output !~ /1500000/) {
    print "FAIL: sv::absrange failed to handle scientific notation with suffix. Output: '$abs_output'\n";
    $failed = 1;
} else {
    print "PASS: sv::absrange handled scientific notation with suffix.\n";
}

print "Testing sv::overlap with scientific notation and suffixes...\n";
if (!sv::overlap("1.5e6___chr1", "1.500001e6___chr1", 100)) {
    print "FAIL: sv::overlap failed to match scientific notation with suffix.\n";
    $failed = 1;
} else {
    print "PASS: sv::overlap handled scientific notation with suffix.\n";
}

print "Testing keysort with scientific notation and suffixes...\n";
{
    our ($a, $b);
    $a = "1.5e6___chr1"; # 1.5M
    $b = "1e5___chr1";   # 100k
    my $res = keysort();
    # Numerical: 1.5e6 > 1e5 => should return 1
    if ($res != 1) {
         print "FAIL: keysort failed to handle scientific notation with suffix. Result: $res\n";
         $failed = 1;
    } else {
         print "PASS: keysort handled scientific notation with suffix.\n";
    }
}

print "Testing keysort with pure numeric values...\n";
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

    $a = "1.5e6";
    $b = "200000";
    $res = keysort();
    if ($res != 1) {
         print "FAIL: keysort failed for pure numeric scientific notation (1.5e6 vs 200000). Result: $res\n";
         $failed = 1;
    } else {
         print "PASS: keysort handled pure numeric scientific notation correctly.\n";
    }
}

exit($failed);
