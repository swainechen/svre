#!/usr/bin/perl
use strict;
use warnings;
use File::Spec;
use lib '.';

BEGIN {
    # Mocking Math::Round
    package Math::Round;
    $INC{"Math/Round.pm"} = "mock";
    sub round { int($_[0] + 0.5) }
    require Exporter;
    @Math::Round::ISA = qw(Exporter);
    @Math::Round::EXPORT = qw(round);

    # Mocking Math::Trig
    package Math::Trig;
    $INC{"Math/Trig.pm"} = "mock";
    sub pi () { 3.14159265358979 }
    @Math::Trig::ISA = qw(Exporter);
    @Math::Trig::EXPORT = qw(pi);

    # Mocking Math::CDF
    package Math::CDF;
    $INC{"Math/CDF.pm"} = "mock";

    # Mocking Math::BigFloat
    package Math::BigFloat;
    $INC{"Math/BigFloat.pm"} = "mock";
}

use sv;

my $RE_INT = $sv::RE_INT;

my @tests = (
    { input => "100___chr1", expected => "chr1" },
    { input => "100..200___chr1", expected => "chr1" },
    { input => "+100___chr2", expected => "chr2" },
    { input => "-100..-50___chr3", expected => "chr3" },
    { input => "100", expected => "100" }, # should not change if it doesn't match the full pattern including ___
);

my $pass = 1;
foreach my $test (@tests) {
    my $val = $test->{input};
    $val =~ s/^${RE_INT}(?:\.\.${RE_INT})?___//;
    if ($val eq $test->{expected}) {
        print "PASS: '$test->{input}' -> '$val'\n";
    } else {
        print "FAIL: '$test->{input}' -> '$val' (expected '$test->{expected}')\n";
        $pass = 0;
    }
}

if ($pass) {
    print "All metadata stripping tests passed!\n";
    exit 0;
} else {
    print "Some metadata stripping tests failed.\n";
    exit 1;
}
