#!/usr/bin/perl
use strict;
use warnings;

# Mocking modules before loading sv.pm to handle missing environment dependencies
BEGIN {
    $INC{"Math/Round.pm"} = 1;
    $INC{"Math/CDF.pm"} = 1;
    $INC{"Math/Trig.pm"} = 1;
    $INC{"Math/BigFloat.pm"} = 1;
    package Math::Round; sub round { $_[0] };
    package Math::CDF; sub pnorm { 0.5 };
    package Math::Trig; sub pi () { 3.14159 };
    *sv::pi = \&Math::Trig::pi;
    package Math::BigFloat; sub new { bless {}, shift };
}

use lib '.';
require sv;
use Test::More tests => 1;

# Test Case: Chromosome present in reference but missing in precount
# This previously caused a crash due to undefined hash dereference in sv::ric
my $refh = { 'chr1' => 1000, 'chr2' => 1000 };
my $precount = {
    'chr1' => {
        100 => { 50 => 10, 'rcount' => 10 }
    }
};
my $global_dist = { 50 => 1.0 };
my $cov_bin = 100;

eval {
    sv::ric($refh, $precount, $cov_bin, $global_dist);
};

ok(!$@, "sv::ric handles chromosomes with missing data gracefully") or diag("sv::ric crashed: $@");

done_testing();
