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

package sv;
use constant pi => 3.14159265358979;

package main;
use lib '.';
require sv;
use Test::More tests => 5;

# Test sv::ric with zero counts
my $refh = { chr1 => 1000 };
my $precount = { chr1 => { 100 => { 50 => 10, "pair" => [200] } } };
my $global_dist = { 50 => 0 }; # Potential exp_prob = 0

eval {
    my ($ri, $count) = sv::ric($refh, $precount, 10, $global_dist);
    ok(1, "sv::ric did not crash with exp_prob = 0");
};
if ($@) {
    fail("sv::ric crashed: $@");
}

# Test sv::bootstrap with zero prob_dist
my $prob_dist = { 50 => 0 };
eval {
    # We need to make sure bootstrap doesn't hang or crash.
    # With prob_dist 0, CDF might be weird, but let's see if our guard in the loop works.
    # Note: bootstrap might have other issues if prob_dist sums to 0, but we specifically fixed the log.
    my $res = sv::bootstrap({50 => 0.1, 100 => 0}, 10, 5, 1);
    ok(1, "sv::bootstrap did not crash with zero prob in distribution");
};
if ($@) {
    fail("sv::bootstrap crashed: $@");
}

# Test sv::poissonApprox
eval {
    my $res = sv::poissonApprox(0, 10);
    is($res, -1000, "sv::poissonApprox handles lambda=0");

    $res = sv::poissonApprox(-1, 10);
    is($res, -1000, "sv::poissonApprox handles lambda < 0");
};

# Test sv::factorialApprox
eval {
    my $res = sv::factorialApprox(undef);
    is($res, 0, "sv::factorialApprox handles undef");
};

done_testing();
