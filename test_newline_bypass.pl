#!/usr/bin/perl
use strict;
use warnings;
use File::Spec;
use lib '.';
use sv;

# Mock $a and $b for keysort testing
our ($a, $b);

# 1. Test isfloat with newlines
print "Testing sv::isfloat with newlines...\n";
if (sv::isfloat("100\n")) {
    print "FAIL: sv::isfloat matched '100\\n'\n";
    exit 1;
} else {
    print "PASS: sv::isfloat rejected '100\\n'\n";
}

# 2. Test sv::overlap suffix stripping
print "Testing sv::overlap suffix stripping with newlines...\n";
# If overlap handles newlines correctly, it should strip '___ref\n' and match '100'
if (sv::overlap("100___ref\n", "100", 0)) {
    print "PASS: sv::overlap stripped '___ref\\n' correctly\n";
} else {
    print "FAIL: sv::overlap failed to strip '___ref\\n' or failed to match\n";
    exit 1;
}

# 3. Test keysort in svre.pl (logic verification)
sub keysort_test {
    $a = $_[0];
    $b = $_[1];
    if (sv::isfloat($a) && sv::isfloat($b)) {
        return $a <=> $b;
    } elsif ($a =~ /^-?\d+___\S+\z/ && $b =~ /^-?\d+___\S+\z/) {
        $a =~ /^(-?\d+)___(\S+\z)/;
        my $p1 = $1;
        my $r1 = $2;
        $b =~ /^(-?\d+)___(\S+\z)/;
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

print "Testing keysort validation with newlines...\n";
my $res = keysort_test("100___ref\n", "100___ref");
if ($res != 0) {
    print "PASS: keysort handled newline in input correctly (did not match complex regex)\n";
} else {
    print "FAIL: keysort incorrectly matched '100___ref\\n' as valid coordinate format\n";
    exit 1;
}

print "ALL NEWLINE BYPASS TESTS PASSED.\n";
exit 0;
