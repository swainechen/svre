use strict;
use warnings;
use Test::More tests => 2;

# Mock Math::CDF::pnorm to avoid dependency issues in the test environment
BEGIN {
    $INC{'Math/Round.pm'} = 1;
    $INC{'Math/CDF.pm'} = 1;
    $INC{'Math/Trig.pm'} = 1;
    $INC{'Math/BigFloat.pm'} = 1;

    # We need to define pi in the main namespace or sv namespace if it's used as a bareword
    package Math::Trig;
    use base 'Exporter';
    our @EXPORT = qw(pi);
    sub pi { return 3.14159; }

    package Math::CDF;
    sub pnorm { return 0.5; }

    package sv;
    # If sv.pm uses pi as a bareword, it might expect it to be imported.
    # But sv.pm has 'use Math::Trig;' which should import it if it's exported.
}

use lib '.';
require 'sv.pm';

subtest 'sv::normal hardening' => sub {
    plan tests => 3;

    eval {
        my $res = sv::normal(10, 5, 0, 1);
        is($res, 0, "sv::normal returns 0 when sd is 0");
    };
    ok(!$@, "sv::normal doesn't crash when sd is 0") or diag("Crashed: $@");

    eval {
        my $res = sv::normal(10, 5, -1, 1);
        is($res, 0, "sv::normal returns 0 when sd is negative");
    };
};

subtest 'sv::normal input validation' => sub {
    plan tests => 2;
    eval {
        my $res = sv::normal(10, 5, "not a number", 1);
        is($res, 0, "sv::normal returns 0 when sd is non-numeric");
    };
    ok(!$@, "sv::normal doesn't crash with non-numeric input") or diag("Crashed: $@");
};
