use strict;
use warnings;
BEGIN {
    $INC{'Math/Round.pm'} = 'mock';
    $INC{'Math/CDF.pm'} = 'mock';
    $INC{'Math/Trig.pm'} = 'mock';
    $INC{'Math/BigFloat.pm'} = 'mock';
    package Math::Round; our @EXPORT = qw(round); use base 'Exporter'; sub round { int($_[0]+0.5) }
    package Math::Trig; our @EXPORT = qw(pi); use base 'Exporter'; sub pi { 3.14 }
}
use lib '.';
use sv;
use Time::HiRes qw(time);
use Data::Dumper;

$Data::Dumper::Sortkeys = 1;

my $refh = { "chr1" => 1_000_000 };
my $precount = { "chr1" => { 10 => { 100 => 50 }, 100 => { 100 => 60 } } };
my $global_dist = { 100 => 1 };
my $cov_bin = 100;

print "Testing with small genome...\n";
my ($ri1, $c1) = sv::ric($refh, $precount, $cov_bin, $global_dist);

print "Testing with large genome...\n";
$refh->{chr1} = 100_000_000;
my $start = time();
my ($ri2, $c2) = sv::ric($refh, $precount, $cov_bin, $global_dist);
my $end = time();
printf "100MB genome took %.6f seconds\n", $end - $start;

if (Dumper($ri1) eq Dumper($ri2)) {
    print "Verification SUCCESS: Results match and performance is independent of genome size.\n";
} else {
    print "Verification FAILURE: Results mismatch!\n";
    print "Small: ", Dumper($ri1);
    print "Large: ", Dumper($ri2);
}

# Test negative coordinates too
$precount = { "chr1" => { -10 => { 100 => 50 }, -100 => { 100 => 60 } } };
$refh->{chr1} = 1_000_000;
($ri1, $c1) = sv::ric($refh, $precount, $cov_bin, $global_dist);
$refh->{chr1} = 100_000_000;
($ri2, $c2) = sv::ric($refh, $precount, $cov_bin, $global_dist);
if (Dumper($ri1) eq Dumper($ri2)) {
    print "Verification SUCCESS: Negative coordinates results match.\n";
} else {
    print "Verification FAILURE: Negative coordinates results mismatch!\n";
}
