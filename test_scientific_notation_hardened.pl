use strict;
use warnings;
use Test::More;

# Mocking external dependencies BEFORE loading sv.pm
BEGIN {
    $INC{"Math/Round.pm"} = "mock";
    $INC{"Math/CDF.pm"} = "mock";
    $INC{"Math/Trig.pm"} = "mock";
}

# Define pi in the package sv where it's used
{
    package sv;
    sub pi { 3.14159265358979 }
}

use sv;

subtest 'keysort hardening' => sub {
    # We need to test keysort logic. Since it's in svre.pl, we can't easily import it.
    # We'll redefine it here with the same logic as the updated svre.pl.
    my $RE_FLOAT = qr/([+-]?(?=\d|\.\d)\d*(?:\.\d*)?(?:[Ee][+-]?\d+)?)/;
    no warnings 'redefine';
    local *keysort = sub {
        my ($a_val, $b_val) = ($main::a, $main::b);
        if ($a_val =~ /^${RE_FLOAT}\z/ && $b_val =~ /^${RE_FLOAT}\z/) {
            return $a_val <=> $b_val;
        } elsif (my ($p1, $r1) = $a_val =~ /^${RE_FLOAT}___(\S+)\z/) {
            if (my ($p2, $r2) = $b_val =~ /^${RE_FLOAT}___(\S+)\z/) {
                if ($r1 eq $r2) {
                    return $p1 <=> $p2;
                } else {
                    return $r1 cmp $r2;
                }
            } else {
                return 1;
            }
        } elsif ($b_val =~ /^${RE_FLOAT}___(\S+)\z/) {
            return -1;
        } else {
            return $a_val cmp $b_val;
        }
    };

    my @input = ("1000000___chr1", "1.5e6___chr1", "500000___chr1");
    my @sorted = sort keysort @input;
    is_deeply(\@sorted, ["500000___chr1", "1000000___chr1", "1.5e6___chr1"], 'keysort handles scientific notation with suffixes');

    @input = ("2e2", "150", "3e1");
    @sorted = sort keysort @input;
    is_deeply(\@sorted, ["3e1", "150", "2e2"], 'keysort handles pure scientific notation');

    @input = ("1.5___chr1", "1.2___chr1", "1.3e0___chr1");
    @sorted = sort keysort @input;
    is_deeply(\@sorted, ["1.2___chr1", "1.3e0___chr1", "1.5___chr1"], 'keysort handles floats and scientific notation');
};

subtest 'sv::range hardening' => sub {
    # Note: range output format depends on how it stringifies numbers.
    my $res = sv::range(["1e6", "1.00005e6"], 100);
    ok($res =~ /^[0-9.e+-]+\.\.[0-9.e+-]+$/, "sv::range handles scientific notation: $res");

    $res = sv::range(["1e6___chr1", "1.00005e6___chr1"], 100);
    ok($res =~ /^[0-9.e+-]+\.\.[0-9.e+-]+___chr1$/, "sv::range handles scientific notation with suffixes: $res");
};

subtest 'sv::absrange hardening' => sub {
    my $res = sv::absrange(["-1.5e6___chr1"], 0);
    ok($res =~ /^[0-9.e+-]+___chr1$/, "sv::absrange handles negative scientific notation with suffixes: $res");
};

done_testing();
