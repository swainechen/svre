#!/usr/bin/perl
use strict;
use warnings;
use Test::More tests => 2;

# This script demonstrates that the regex-based SAM tag validation
# in svre.pl can be bypassed if the spoofed tag appears in the QUAL field.

my $mapper = "default";
my $mapper_command_line = "";

# Current logic in svre.pl (simplified for the test)
sub check_line_vulnerable {
    my ($line) = @_;
    # This is the "else" branch in svre.pl's mapper check
    return 0 if $line !~ /\tNM:i:0(\t|\r?\n|$)/;

    my @f = split /\t/, $line;
    # FLAG 4 is unmapped
    return 0 if $f[1] & 4;
    return 1;
}

# Proposed hardened logic
sub check_line_hardened {
    my ($line) = @_;
    chomp $line;
    my @f = split /\t/, $line;
    return 0 if $f[1] & 4;

    # Check tags in fields 12+ (index 11+)
    my $found_nm = 0;
    for (my $i = 11; $i <= $#f; $i++) {
        if ($f[$i] eq "NM:i:0") {
            $found_nm = 1;
            last;
        }
    }
    return $found_nm;
}

# A SAM line where QUAL field is "NM:i:0" but there are no tags.
# Field 11 is QUAL.
my $spoofed_line = "read1\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACTG\tNM:i:0\n";

ok(check_line_vulnerable($spoofed_line), "Vulnerable logic incorrectly accepts spoofed QUAL field as NM:i:0 tag");
ok(!check_line_hardened($spoofed_line), "Hardened logic correctly rejects spoofed QUAL field");

done_testing();
