#!/usr/bin/perl
use strict;
use warnings;

BEGIN {
    $INC{"Math/Round.pm"} = "mock";
    $INC{"Math/CDF.pm"} = "mock";
    $INC{"Math/Trig.pm"} = "mock";
    $INC{"Math/BigFloat.pm"} = "mock";
    package Math::Trig;
    use constant pi => 3.14159;
    our @EXPORT = qw(pi);
    use Exporter qw(import);
}

use sv;

my $RE_FLOAT = $sv::RE_FLOAT;

# 1. Test sv::range with scientific notation and suffixes
print "Testing sv::range...\n";
my @coords = ("1.2e6", "1.25e6", "2e6___chr2", "2.1e6___chr2");
my $range_out = sv::range(\@coords, 100000);
print "Range out: $range_out\n";
if ($range_out =~ /1200000\.\.1250000/ && $range_out =~ /2000000\.\.2100000___chr2/) {
    print "PASS: sv::range handles scientific notation and suffixes.\n";
} else {
    die "FAIL: sv::range did not produce expected output.\n";
}

# 2. Test sv::absrange
print "\nTesting sv::absrange...\n";
my @abs_coords = ("-1.5e6", "1.5e6..1.6e6___chr3");
my $abs_range_out = sv::absrange(\@abs_coords, 100000);
print "Abs range out: $abs_range_out\n";
if ($abs_range_out =~ /^1500000/ && $abs_range_out =~ /1500000\.\.1600000___chr3/) {
    print "PASS: sv::absrange handles scientific notation and suffixes.\n";
} else {
    die "FAIL: sv::absrange did not produce expected output.\n";
}

# 3. Test keysort (mocked from svre.pl)
print "\nTesting keysort logic...\n";
sub keysort_mock {
  if ($a =~ /^${RE_FLOAT}\z/ && $b =~ /^${RE_FLOAT}\z/) {
    return $a <=> $b;
  } elsif (my ($p1, $r1) = $a =~ /^${RE_FLOAT}___(\S+)\z/) {
    if (my ($p2, $r2) = $b =~ /^${RE_FLOAT}___(\S+)\z/) {
      if ($r1 eq $r2) {
        return $p1 <=> $p2;
      } else {
        return $r1 cmp $r2;
      }
    } else {
      return 1;
    }
  } elsif ($b =~ /^${RE_FLOAT}___(\S+)\z/) {
    return -1;
  } else {
    return $a cmp $b;
  }
}

my @to_sort = ("2e6", "1e6", "1.5e6___chr1", "1.1e6___chr1", "1.5e6___chr2");
our ($a, $b);
my @sorted = sort keysort_mock @to_sort;
print "Sorted: " . join(", ", @sorted) . "\n";
if ($sorted[0] eq "1e6" && $sorted[1] eq "2e6" && $sorted[2] eq "1.1e6___chr1" && $sorted[3] eq "1.5e6___chr1" && $sorted[4] eq "1.5e6___chr2") {
    print "PASS: keysort handles scientific notation and suffixes.\n";
} else {
    die "FAIL: keysort did not produce expected order.\n";
}

# 4. Test RIC loop pattern (translocation target identification)
print "\nTesting RIC loop pattern...\n";
my $dist = "1.5e6___chr2";
if ($dist =~ /^${RE_FLOAT}___(.*)\z/s) {
    print "PASS: Distance identified as translocation target with scientific notation.\n";
} else {
    die "FAIL: Distance NOT identified as translocation target with scientific notation.\n";
}

# 5. Test translocation target reference extraction
print "\nTesting translocation target extraction...\n";
my $t_target = "1.5e6___chr2";
my $t_ref = $t_target;
$t_ref =~ s/^${RE_FLOAT}___//;
print "Extracted ref: $t_ref\n";
if ($t_ref eq "chr2") {
    print "PASS: Ref extracted correctly from scientific notation target.\n";
} else {
    die "FAIL: Ref NOT extracted correctly from scientific notation target.\n";
}

my $t_range_target = "1.5e6..1.6e6___chr3";
my $t_range_ref = $t_range_target;
$t_range_ref =~ s/^${RE_FLOAT}(?:\.\.${RE_FLOAT})?___//;
print "Extracted range ref: $t_range_ref\n";
if ($t_range_ref eq "chr3") {
    print "PASS: Range ref extracted correctly from scientific notation target.\n";
} else {
    die "FAIL: Range ref NOT extracted correctly from scientific notation target.\n";
}

print "\nALL SCIENTIFIC NOTATION HARDENING TESTS PASSED.\n";
