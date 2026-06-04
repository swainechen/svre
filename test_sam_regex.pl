#!/usr/bin/perl
use strict;
use warnings;
use Test::More tests => 5;

# This script verifies that the updated regexes in svre.pl correctly
# distinguish between SAM tags and spoofed tags in the read name,
# including support for CRLF line endings.

my @test_regexes = (
    qr/\tNM:i:0(\t|\r?\n|$)/,
    qr/\tXS:i:0(\t|\r?\n|$)/,
    qr/\tXT:A:U(\t|\r?\n|$)/,
);

subtest 'Positive matches (Valid SAM lines)' => sub {
    plan tests => 3;
    my $good_line = "read1\t0\tchr1\t100\t60\t75M\t*\t0\t0\tSEQ\tQUAL\tNM:i:0\tXS:i:0\tXT:A:U\n";
    ok($good_line =~ $test_regexes[0], "NM:i:0 matches in good line");
    ok($good_line =~ $test_regexes[1], "XS:i:0 matches in good line");
    ok($good_line =~ $test_regexes[2], "XT:A:U matches in good line");
};

subtest 'CRLF support' => sub {
    plan tests => 1;
    my $crlf_line = "read1\t0\tchr1\t100\t60\t75M\t*\t0\t0\tSEQ\tQUAL\tNM:i:0\r\n";
    ok($crlf_line =~ $test_regexes[0], "NM:i:0 matches with CRLF line ending");
};

subtest 'Negative matches (Spoofed in read name)' => sub {
    plan tests => 3;
    my $bad_line = "read_NM:i:0_XS:i:0_XT:A:U\t0\tchr1\t100\t60\t75M\t*\t0\t0\tSEQ\tQUAL\tNM:i:5\n";
    ok($bad_line !~ $test_regexes[0], "NM:i:0 does NOT match in bad line");
    ok($bad_line !~ $test_regexes[1], "XS:i:0 does NOT match in bad line");
    ok($bad_line !~ $test_regexes[2], "XT:A:U does NOT match in bad line");
};

subtest 'Negative matches (Partial tag match)' => sub {
    plan tests => 2;
    my $partial_line = "read1\t0\tchr1\t100\t60\t75M\t*\t0\t0\tSEQ\tQUAL\tNM:i:00\tXS:i:0.5\n";
    ok($partial_line !~ $test_regexes[0], "NM:i:00 does NOT match NM:i:0");
    ok($partial_line !~ $test_regexes[1], "XS:i:0.5 does NOT match XS:i:0");
};

subtest 'End of line matches' => sub {
    plan tests => 1;
    my $eol_line = "read1\t0\tchr1\t100\t60\t75M\t*\t0\t0\tSEQ\tQUAL\tNM:i:0";
    ok($eol_line =~ $test_regexes[0], "NM:i:0 matches at end of string without newline");
};

done_testing();
