use strict;
use warnings;

# Mocking dependencies
BEGIN {
    $INC{"Math/Round.pm"} = 1;
    package Math::Round;
    use Exporter 'import';
    our @EXPORT = qw(round);
    sub round { return int($_[0] + 0.5); }

    $INC{"Math/CDF.pm"} = 1;
    package Math::CDF;
    sub pnorm { return 0.5; }

    $INC{"Math/Trig.pm"} = 1;
    package Math::Trig;
    use Exporter 'import';
    our @EXPORT = qw(pi);
    sub pi { return 3.14159; }

    $INC{"Math/BigFloat.pm"} = 1;
    package Math::BigFloat;
    sub new { return bless {}, $_[0]; }
    sub bpow { return $_[0]; }
    sub bfac { return $_[0]; }
    sub bmul { return $_[0]; }
    sub bdiv { return $_[0]; }
}

use lib '.';
require 'sv.pm';
use Data::Dumper;

my $refh = { "chr1" => 1000 };
# We need to make sure we only have ONE bin.
# ric bins by coverage. If we have few reads, we get one bin.
# Note: ric expects precount->{$ref}->{$i}->{pair} to be an array ref if defined.
my $precount = { "chr1" => { 10 => { 100 => 10, "pair" => [] } } };
my $cov_bin = 100; # rcount 10 < cov_bin/2 (50)
my $global_dist = { 100 => 1 };

my ($ri, $count) = sv::ric($refh, $precount, $cov_bin, $global_dist);
print "RI for chr1: ", Dumper($ri->{chr1});
print "Count for chr1: ", Dumper($count->{chr1});
