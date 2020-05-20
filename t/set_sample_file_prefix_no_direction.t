#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw { :all };
use Modern::Perl qw{ 2018 };

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.01;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::File_info}      => [qw{ set_sample_file_prefix_no_direction }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File_info qw{ set_sample_file_prefix_no_direction };

diag(   q{Test set_sample_file_prefix_no_direction from File_info.pm v}
      . $MIP::File_info::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a sample_id and an file_prefix_no_direction
my %file_info;
my $file_prefix_no_direction = q{7_161011_HHJJCCCXY_ADM1059A1_NAATGCGC};
my $sample_id                = q{sample_id};

## When sequence run type is single-end
my $sequence_run_type = q{single-end};

## Then do not set file_prefix_no_direction in file_info hash
set_sample_file_prefix_no_direction(
    {
        file_info_href    => \%file_info,
        mip_file_format   => $file_prefix_no_direction,
        sample_id         => $sample_id,
        sequence_run_type => $sequence_run_type,
    }
);

## Then set file_prefix_no_direction in file_info hash
is( $file_info{$sample_id}{file_prefix_no_direction}{$file_prefix_no_direction},
    $sequence_run_type, q{Set sequence run type for file_prefix_no_direction} );

done_testing();