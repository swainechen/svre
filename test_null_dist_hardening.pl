use strict;
use warnings;
use Test::More tests => 1;

BEGIN {
    $INC{'Math/Round.pm'} = 1;
    $INC{'Math/CDF.pm'} = 1;
    $INC{'Math/Trig.pm'} = 1;
    $INC{'Math/BigFloat.pm'} = 1;

    package Math::Trig;
    use base 'Exporter';
    our @EXPORT = qw(pi);
    sub pi { return 3.14159; }

    package Math::CDF;
    sub pnorm { return 0.5; }
}

use lib '.';
require 'sv.pm';

subtest 'sv::null_distribution hardening' => sub {
    plan tests => 3;

    eval {
        my $res = sv::null_distribution(100, 50, 0, 10);
        is_deeply($res, {}, "sv::null_distribution returns empty hash when sd is 0");
    };
    ok(!$@, "sv::null_distribution doesn't crash when sd is 0") or diag("Crashed: $@");

    eval {
        my $res = sv::null_distribution(100, 50, -1, 10);
        is_deeply($res, {}, "sv::null_distribution returns empty hash when sd is negative");
    };
};
