use strict;
use warnings;
my $script = 'svre.pl';
my @tests = (
    { name => "bootstrap high", args => ["-bootstrap", "100000001"], expected => qr/Error: bootstrap must be a valid number between 1 and 100,000,000/ },
    { name => "cov high", args => ["-cov", "1000001"], expected => qr/Error: cov must be a valid number between 1 and 1,000,000/ },
    { name => "bootstrap ok", args => ["-bootstrap", "100000000"], expected => qr/Error: samtools not found in PATH|Usage:/ },
    { name => "cov ok", args => ["-cov", "1000000"], expected => qr/Error: samtools not found in PATH|Usage:/ },
);

foreach my $t (@tests) {
    print "Testing $t->{name}... ";
    # Using a one-liner to mock dependencies at BEGIN time for svre.pl
    my $mock = 'BEGIN { package sv; sub isfloat { return $_[0] =~ /^([+-])?(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/ } $INC{"sv.pm"} = "mock"; $INC{"Math/Round.pm"} = "mock"; $INC{"Math/CDF.pm"} = "mock"; $INC{"Math/Trig.pm"} = "mock"; $INC{"Math/BigFloat.pm"} = "mock"; }';
    my $cmd = "perl -e '$mock' -S $script " . join(" ", @{$t->{args}});
    # Actually, svre.pl is a script, so I should probably do it differently
    my $full_cmd = "perl -I. -e '$mock do \"$script\" or die \$\@' -- " . join(" ", @{$t->{args}});
    my $out = `$full_cmd 2>&1`;
    if ($out =~ $t->{expected}) {
        print "SUCCESS\n";
    } else {
        print "FAILURE\n";
        print "Output: $out\n";
    }
}
