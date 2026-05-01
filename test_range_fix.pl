use strict;
use warnings;

BEGIN {
    package Math::Round;
    $INC{'Math/Round.pm'} = 'mock';
    our @EXPORT = qw(round);
    use base 'Exporter';
    sub round { return int($_[0] + 0.5); }

    package Math::Trig;
    $INC{'Math/Trig.pm'} = 'mock';
    our @EXPORT = qw(pi);
    use base 'Exporter';
    sub pi { return 3.14159265358979; }

    package Math::CDF;
    $INC{'Math/CDF.pm'} = 'mock';

    package Math::BigFloat;
    $INC{'Math/BigFloat.pm'} = 'mock';
}

use lib '.';
use sv;

# Test case 1: Single point
my $r1 = sv::range(["100"], 1);
die "Test 1 failed: expected 100, got $r1" unless $r1 eq "100";

# Test case 2: Simple range
my $r2 = sv::range(["100..200"], 1);
die "Test 2 failed: expected 100..200, got $r2" unless $r2 eq "100..200";

# Test case 3: Overlapping ranges
my $r3 = sv::range(["100..200", "150..250"], 1);
die "Test 3 failed: expected 100..250, got $r3" unless $r3 eq "100..250";

# Test case 4: Disjoint ranges
my $r4 = sv::range(["100..200", "300..400"], 1);
die "Test 4 failed: expected 100..200,300..400, got $r4" unless $r4 eq "100..200,300..400";

# Test case 5: Ranges with different suffixes, merging with cluster distance
my $r5 = sv::range(["100___ref1", "150..250___ref1", "300___ref2"], 60);
# 100 and 150..250 should merge because 150 <= 100 + 60
die "Test 5 failed: got $r5" unless $r5 =~ /100\.\.250___ref1/ && $r5 =~ /300___ref2/;

# Test case 6: Large range (The DoS test)
my $large_binsize = 10_000_000;
my $large_range = ["0..$large_binsize"];
print "Testing large range (10M bp)...\n";
my $start_time = time();
my $r6 = sv::range($large_range, 1);
my $end_time = time();
my $duration = $end_time - $start_time;
print "Large range took $duration seconds.\n";
die "Test 6 failed: large range took too long ($duration s)" if $duration > 5;
die "Test 6 failed: expected 0..$large_binsize, got $r6" unless $r6 eq "0..$large_binsize";

# Test case 7: Merging with cluster distance
my $r7 = sv::range(["100", "110"], 15);
die "Test 7 failed: expected 100..110, got $r7" unless $r7 eq "100..110";

print "All sv::range tests passed!\n";
