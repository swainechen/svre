use strict;
use warnings;
use Test::More tests => 2;

BEGIN {
    $INC{'Math/Round.pm'} = 1;
    $INC{'Math/CDF.pm'} = 1;
    $INC{'Math/Trig.pm'} = 1;
    $INC{'Math/BigFloat.pm'} = 1;

    package Math::Round;
    sub round { return int($_[0] + 0.5); }

    package Math::CDF;
    sub pnorm { return 0.5; }

    package Math::Trig;
    use base 'Exporter';
    our @EXPORT = qw(pi);
    sub pi { return 3.14159; }
}

use lib '.';
use sv;

subtest 'sv::pval_BY empty input' => sub {
    my $res = eval { sv::pval_BY([]) };
    is($@, '', "sv::pval_BY doesn't crash on empty input");
    is_deeply($res, {}, "sv::pval_BY returns empty hash on empty input");
};

subtest 'sv::pval_BY normal input' => sub {
    # p-values must be sorted in increasing order for sv::pval_BY
    my $pvals = [0.01, 0.05, 0.5];
    my $res = sv::pval_BY($pvals);
    ok(exists $res->{'0.01'}, "Contains key 0.01");
    ok(exists $res->{'0.05'}, "Contains key 0.05");
    ok(exists $res->{'0.5'}, "Contains key 0.5");
};
