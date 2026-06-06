use strict;
use warnings;
my $script = 'svre.pl';
my @tests = (
    { name => "start non-numeric", args => ["-start", "abc"], expected => qr/Value "abc" invalid for option start|Error parsing command line options/ },
    { name => "end non-numeric", args => ["-end", "def"], expected => qr/Value "def" invalid for option end|Error parsing command line options/ },
    { name => "bootstrap iterations over limit", args => ["-bootstrap", "1000001", "-cov", "1000"], expected => qr/Error: total bootstrap iterations \(bootstrap \* cov\) exceeds limit of 1,000,000,000/ },
    { name => "bootstrap iterations at limit", args => ["-bootstrap", "1000000", "-cov", "1000"], expected => qr/Error: samtools not found in PATH|Usage:/ },
);

foreach my $t (@tests) {
    print "Testing $t->{name}... ";
    # Using a one-liner to mock dependencies at BEGIN time for svre.pl
    my $mock = 'BEGIN { package sv; sub isfloat { return $_[0] =~ /^([+-])?(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/ } $INC{"sv.pm"} = "mock"; $INC{"Math/Round.pm"} = "mock"; $INC{"Math/CDF.pm"} = "mock"; $INC{"Math/Trig.pm"} = "mock"; $INC{"Math/BigFloat.pm"} = "mock"; }';
    my $full_cmd = "perl -I. -e '$mock do \"$script\" or die \$\@' -- " . join(" ", @{$t->{args}});
    my $out = `$full_cmd 2>&1`;
    if ($out =~ $t->{expected}) {
        print "SUCCESS\n";
    } else {
        print "FAILURE\n";
        print "Output: $out\n";
    }
}
