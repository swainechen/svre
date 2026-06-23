#!/usr/bin/perl
use strict;
use warnings;

sub parse_tags {
    my $line = shift;
    $line =~ s/[\r\n]+\z//;
    my @f = split /\t/, $line;
    my %tags = ();
    for (my $tag_idx = 11; $tag_idx <= $#f; $tag_idx++) {
        if ($f[$tag_idx] =~ /^([^:]+:[^:]+):(.*)\z/) {
            $tags{$1} = $2;
        }
    }
    return \%tags;
}

print "Testing SAM tag parsing with CRLF...\n";
my $line_crlf = "read1\t0\tchr1\t100\t60\t75M\t*\t0\t0\tA\tB\tNM:i:0\r\n";
my $tags = parse_tags($line_crlf);
if (exists $tags->{"NM:i"} && $tags->{"NM:i"} eq "0") {
    print "PASS: Tag matched correctly with CRLF.\n";
} else {
    print "FAIL: Tag did not match with CRLF! (Value: " . ($tags->{"NM:i"} // "undef") . ")\n";
    exit 1;
}

print "Testing SAM tag parsing with LF only...\n";
my $line_lf = "read1\t0\tchr1\t100\t60\t75M\t*\t0\t0\tA\tB\tNM:i:0\n";
$tags = parse_tags($line_lf);
if (exists $tags->{"NM:i"} && $tags->{"NM:i"} eq "0") {
    print "PASS: Tag matched correctly with LF.\n";
} else {
    print "FAIL: Tag did not match with LF!\n";
    exit 1;
}

print "Testing file extension check with \\z...\n";
my $valid_bam = "test.bam";
if ($valid_bam =~ /\.bam\z/i) {
    print "PASS: Valid filename matched extension.\n";
} else {
    print "FAIL: Valid filename did not match extension!\n";
    exit 1;
}

my $bam_with_newline = "test.bam\n";
if ($bam_with_newline =~ /\.bam\z/i) {
    print "FAIL: Filename with newline matched extension! (Potential bypass)\n";
    exit 1;
} else {
    print "PASS: Filename with newline rejected.\n";
}

print "ALL NEWLINE VALIDATION TESTS PASSED.\n";
exit 0;
