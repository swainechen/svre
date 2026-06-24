#!/usr/bin/perl
use strict;
use warnings;
use Test::More;
use File::Basename;
use File::Spec;

# Mocking external dependencies
BEGIN {
    $INC{"Math/Round.pm"} = "mock";
    *Math::Round::round = sub { return int($_[0] + 0.5); };
    $INC{"Math/CDF.pm"} = "mock";
    $INC{"Math/Trig.pm"} = "mock";
    *Math::Trig::pi = sub { return 3.14159265358979; };
    *sv::pi = sub { return 3.14159265358979; };
    $INC{"Math/BigFloat.pm"} = "mock";
}

use lib '.';
require 'sv.pm';

# Test keysort logic (extracted from svre.pl)
sub keysort_mock {
    my ($a_val, $b_val) = @_;
    local our $a = $a_val;
    local our $b = $b_val;

    if (sv::isfloat($a) && sv::isfloat($b)) {
        return $a <=> $b;
    } elsif ($a =~ /^([+-]?(?=\d|\.\d)\d*(?:\.\d*)?(?:[Ee][+-]?\d+)?)___(\S+)\z/) {
        my ($p1, $r1) = ($1, $2);
        if ($b =~ /^([+-]?(?=\d|\.\d)\d*(?:\.\d*)?(?:[Ee][+-]?\d+)?)___(\S+)\z/) {
            my ($p2, $r2) = ($1, $2);
            if ($r1 eq $r2) {
                return $p1 <=> $p2;
            } else {
                return $r1 cmp $r2;
            }
        } else {
            return 1;
        }
    } elsif ($b =~ /^([+-]?(?=\d|\.\d)\d*(?:\.\d*)?(?:[Ee][+-]?\d+)?)___(\S+)\z/) {
        return -1;
    } else {
        return $a cmp $b;
    }
}

subtest 'keysort scientific notation' => sub {
    is(keysort_mock("1.5e6___chr1", "2.5e6___chr1"), -1, "Scientific notation same chromosome");
    is(keysort_mock("2.5e6___chr1", "1.5e6___chr1"), 1, "Scientific notation same chromosome reverse");
    is(keysort_mock("100___chr1", "1e2___chr1"), 0, "Integer vs Scientific same value");
    is(keysort_mock("1.5e6___chr1", "1.5e6___chr2"), -1, "Scientific notation different chromosome");
    is(keysort_mock("1.5e6___chr2", "1.5e6___chr1"), 1, "Scientific notation different chromosome reverse");
    is(keysort_mock("1.5e6", "2.5e6"), -1, "Scientific notation no suffix");
};

subtest 'sv::range scientific notation' => sub {
    my @input = ("1.5e6___chr1..2e6___chr1", "2.5e6___chr1");
    my $res = sv::range(\@input, 100);
    like($res, qr/1500000\.\.2000000___chr1/, "sv::range parses scientific notation ranges");
    like($res, qr/2500000___chr1/, "sv::range parses scientific notation points");

    @input = ("-1.5e6___chr1..-1e6___chr1");
    $res = sv::range(\@input, 100);
    like($res, qr/-1500000\.\.-1000000___chr1/, "sv::range parses negative scientific notation ranges");
};

subtest 'sv::absrange scientific notation' => sub {
    my @input = ("-1.5e6___chr1..-1e6___chr1", "2e6___chr1");
    my $res = sv::absrange(\@input, 100);
    like($res, qr/1000000\.\.1500000___chr1/, "sv::absrange handles negative scientific notation");
    like($res, qr/2000000___chr1/, "sv::absrange handles positive scientific notation");
};

subtest 'svre.pl regex population loop' => sub {
    my $dist = "1.5e6___chr2";
    my $ywin = 1000;
    my ($location, $name, $j);
    if ($dist =~ /^([+-]?(?=\d|\.\d)\d*(?:\.\d*)?(?:[Ee][+-]?\d+)?)___(.*)\z/s) {
        $location = $1;
        $name = $2;
        $j = int ($location/$ywin) * $ywin;
    }
    is($location, "1.5e6", "Regex extracts scientific location");
    is($name, "chr2", "Regex extracts reference name");
    is($j, 1500000, "Calculation with scientific notation location works");
};

subtest 'svre.pl translocation target parsing' => sub {
    my $dist = "1.5e6___chr2";
    my $t_ref = $dist;
    $t_ref =~ s/^([+-]?(?=\d|\.\d)\d*(?:\.\d*)?(?:[Ee][+-]?\d+)?)___//;
    my $t_pos = $dist;
    $t_pos =~ s/___.*\z//s;

    is($t_ref, "chr2", "Translocation target ref extraction handles scientific notation");
    is($t_pos, "1.5e6", "Translocation target pos extraction handles scientific notation");

    my $target = "1.5e6..2e6___chr2";
    my $extracted_ref = $target;
    $extracted_ref =~ s/^([+-]?(?=\d|\.\d)\d*(?:\.\d*)?(?:[Ee][+-]?\d+)?)(?:\.\.([+-]?(?=\d|\.\d)\d*(?:\.\d*)?(?:[Ee][+-]?\d+)?))?___//;
    is($extracted_ref, "chr2", "Reference extraction from range handles scientific notation");
};

done_testing();
