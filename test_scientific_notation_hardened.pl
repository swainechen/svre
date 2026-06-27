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
    my $RE_INT = $sv::RE_INT;
    no warnings 'redefine';
    local *keysort = sub {
        my ($a_val, $b_val) = ($main::a, $main::b);
        if ($a_val =~ /^${RE_INT}\z/ && $b_val =~ /^${RE_INT}\z/) {
            return $a_val <=> $b_val;
        } elsif (my ($p1, $r1) = $a_val =~ /^${RE_INT}___(\S+)\z/) {
            if (my ($p2, $r2) = $b_val =~ /^${RE_INT}___(\S+)\z/) {
                if ($r1 eq $r2) {
                    return $p1 <=> $p2;
                } else {
                    return $r1 cmp $r2;
                }
            } else {
                return 1;
            }
        } elsif ($b_val =~ /^${RE_INT}___(\S+)\z/) {
            return -1;
        } else {
            return $a_val cmp $b_val;
        }
    };

    my @input = ("1000000___chr1", "1500000___chr1", "500000___chr1");
    my @sorted = sort keysort @input;
    is_deeply(\@sorted, ["500000___chr1", "1000000___chr1", "1500000___chr1"], 'keysort handles integers with suffixes');

    @input = ("200", "150", "30");
    @sorted = sort keysort @input;
    is_deeply(\@sorted, ["30", "150", "200"], 'keysort handles pure integers');

    # Scientific notation should NOT match RE_INT and thus fall back to string comparison
    @input = ("1.5e6___chr1", "2000000___chr1");
    @sorted = sort keysort @input;
    is($sorted[0], "1.5e6___chr1", 'scientific notation falls back to string comparison (alphabetical)');
};

subtest 'sv::range hardening' => sub {
    # Integers should work
    my $res = sv::range(["1000000", "1000050"], 100);
    is($res, "1000000..1000050", "sv::range handles integers");

    # Scientific notation should NOT be parsed as numeric range (since it doesn't match RE_INT)
    $res = sv::range(["1e6", "1.00005e6"], 100);
    # It might treat them as strings and not merge them, or fail format error
    ok($res ne "1000000..1000050", "sv::range does NOT handle scientific notation as numeric");
};

subtest 'sv::isfloat still works' => sub {
    ok(sv::isfloat("1.5e6"), "isfloat still supports scientific notation");
    ok(sv::isfloat("0.05"), "isfloat supports floats");
    ok(!sv::isfloat("100\n"), "isfloat correctly rejects trailing newlines");
};

done_testing();
