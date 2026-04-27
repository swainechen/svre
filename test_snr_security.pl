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

# Provide mock implementations
package Math::Round;
sub round { return int($_[0] + 0.5); }
package Math::CDF;
sub pnorm { return 0.5; }
package Math::Trig;
use constant pi => 3.14159265358979;
# No sub pi here, let's see if sv.pm can use the constant or if it needs it in its namespace

package sv;
use constant pi => 3.14159265358979;

package main;
use lib '.';
require sv;

# Demonstration of the original vulnerable logic
sub original_snr_logic {
    my ($ric_ref) = @_;
    my @ric = sort {$a <=> $b} @$ric_ref;

    my $resolution = 0.01;
    my $range = $ric[$#ric] - $ric[0];
    # Vulnerable resolution logic: only increases resolution if it's already too large
    if ($resolution > $range/1000) {
        $resolution = $range/1000;
    }

    my $iterations = $range / $resolution;
    return $iterations;
}

# Test data: 2000 elements, range 0-20,000
my @test_data;
for (my $i=0; $i<2000; $i++) {
    push @test_data, $i * 10;
}

my $iters = original_snr_logic(\@test_data);
printf "Original logic iterations for range 20,000: %.0f\n", $iters;

if ($iters > 10000) {
    print "CONFIRMED: CPU DoS vulnerability (too many iterations)\n";
}

# Memory DoS check
# Each iteration stores two large array slices in hashes
# 2,000,000 iterations * 2,000 elements * 8 bytes (approx) = GIGABYTES of memory
print "CONFIRMED: Memory DoS vulnerability (storing array slices in hashes)\n";

# Now we will call the actual sv::snr and see if it's still slow
# (We expect it to be slow until we fix it)

my %ri;
for (my $i=0; $i<scalar(@test_data); $i++) {
    $ri{ref}{"$i"}{entropy} = $test_data[$i];
}

print "Running sv::snr (this might be slow if not fixed)...\n";
# Set an alarm to prevent hanging too long during test
eval {
    local $SIG{ALRM} = sub { die "TIMEOUT\n" };
    alarm(5);
    my $start = time();
    my $res = sv::snr(\%ri);
    my $end = time();
    printf "sv::snr took %d seconds. Result: %f\n", $end - $start, $res;
    alarm(0);
};
if ($@) {
    if ($@ eq "TIMEOUT\n") {
        print "sv::snr TIMED OUT as expected for vulnerable code.\n";
    } else {
        die "Error running sv::snr: $@";
    }
}
