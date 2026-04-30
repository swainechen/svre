#!/usr/bin/perl
use strict;
use warnings;

# Mocking modules before loading sv.pm
BEGIN {
    $INC{'Math/Round.pm'} = 1;
    $INC{'Math/CDF.pm'} = 1;
    $INC{'Math/Trig.pm'} = 1;
    $INC{'Math/BigFloat.pm'} = 1;
}

package Math::Trig;
sub pi { return 3.14159265358979; }

package sv;
sub isfloat { return 1; } # simplified
use constant pi => 3.14159265358979;

package main;
use lib '.';
require sv;

my @array = (1, 2, "nan", 4, 5);
# Note: "nan" in Perl might be treated as string or number depending on context.
# Let's use an actual NaN if possible.
my $nan = 0 + "nan";

@array = (1, 2, $nan, 4, 5);
# The array must be sorted for binary_search, but NaN doesn't sort well.
# sv::snr sorts the array: @ric = sort {$a <=> $b} @ric;
@array = sort {$a <=> $b} @array;

print "Array: @array\n";

print "Calling sv::binary_search with NaN in array...\n";
eval {
    local $SIG{ALRM} = sub { die "TIMEOUT\n" };
    alarm(3);
    sv::binary_search(\@array, 3);
    alarm(0);
};
if ($@) {
    print "Result: $@";
} else {
    print "Completed without timeout.\n";
}
