use strict;
use warnings;
use lib '.';

# Mock dependencies for sv.pm
BEGIN {
    $INC{'Math/Round.pm'} = 1;
    {
        package Math::Round;
        sub round { return int($_[0] + 0.5); }
    }
    $INC{'Math/CDF.pm'} = 1;
    {
        package Math::CDF;
        sub pnorm { return 0.5; }
    }
    $INC{'Math/Trig.pm'} = 1;
    {
        package Math::Trig;
        use constant pi => 3.14159265358979;
        sub import {
            no strict 'refs';
            *{(caller)[0].'::pi'} = \&pi;
        }
    }
    $INC{'Math/BigFloat.pm'} = 1;
    {
        package Math::BigFloat;
        sub new { return bless {}, shift; }
    }
}

require sv;

my $failed = 0;

sub assert_overlap {
    my ($r, $s, $range, $expected, $msg) = @_;
    my $result = sv::overlap($r, $s, $range);
    if ($result == $expected) {
        print "OK: $msg\n";
    } else {
        print "FAIL: $msg (Expected $expected, got $result)\n";
        $failed = 1;
    }
}

print "Starting tests for sv::overlap...\n";

# Existing cases
assert_overlap("100", "100", 0, 1, "Point 100 overlaps with 100");
assert_overlap("100___ref", "100___ref", 0, 1, "Point 100___ref overlaps with itself");
assert_overlap("100___ref", "1", 0, 0, "Point 100___ref should NOT overlap with 1");
assert_overlap("100___ref", "90..110", 0, 1, "Point 100___ref overlaps with range 90..110");
assert_overlap("100___ref", "90..110___ref", 0, 1, "Point 100___ref overlaps with range 90..110___ref");

# Regression cases from code review
assert_overlap("1.5..2.5", "2.0", 0, 1, "Float range 1.5..2.5 overlaps with 2.0");
assert_overlap("1.5..2.5", "3.0", 0, 0, "Float range 1.5..2.5 should NOT overlap with 3.0");
assert_overlap("100___ref..200___ref", "150", 0, 1, "Range 100___ref..200___ref overlaps with 150");
assert_overlap("1e2..2e2", "150", 0, 1, "Scientific notation 1e2..2e2 overlaps with 150");

if ($failed) {
    print "\nSome tests FAILED.\n";
    exit 1;
} else {
    print "\nAll tests PASSED.\n";
    exit 0;
}
