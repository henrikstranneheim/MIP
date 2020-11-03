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
use Readonly;
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $COMMA => q{,};
Readonly my $SPACE => q{ };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Get::File}      => [qw{ get_exome_target_bed_file }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Get::File qw{ get_exome_target_bed_file };

diag(   q{Test get_exome_target_bed_file from File.pm v}
      . $MIP::Get::File::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Create log
my $log = test_log( {} );

## Base arguments
my $sample_1         = q{sample_1};
my $sample_2         = q{sample_2};
my $sample_3         = q{sample_3};
my %exome_target_bed = (
    q{exome_target_file.bed}   => $sample_1 . $COMMA . $sample_2,
    q{exome_target_file_2.bed} => $sample_3,
);

## Given a sample_id
my $exome_target_bed_file = get_exome_target_bed_file(
    {
        exome_target_bed_href => \%exome_target_bed,
        sample_id             => $sample_1,
    }
);

## Then return exome target bed file for sample
is( q{exome_target_file.bed}, $exome_target_bed_file,
    q{Get exom target bed file for sample_1} );

## Given a sample id and a file ending
$exome_target_bed_file = get_exome_target_bed_file(
    {
        exome_target_bed_href => \%exome_target_bed,
        sample_id             => $sample_1,
        file_ending           => q{.interval_list},
    }
);

## Then return exome target file with interval list file ending
is( q{exome_target_file.bed.interval_list},
    $exome_target_bed_file, q{Get exom target bed.interval_list file for sample_1} );

## Given another sample id
$exome_target_bed_file = get_exome_target_bed_file(
    {
        exome_target_bed_href => \%exome_target_bed,
        sample_id             => $sample_3,
    }
);

## Then return exome target bed file for sample id
is( q{exome_target_file_2.bed}, $exome_target_bed_file,
    q{Get exom target bed file for sample_3} );

## Given a not exisitng sample id
trap {
    get_exome_target_bed_file(
        {
            exome_target_bed_href => \%exome_target_bed,
            sample_id             => q{not_a_sample_id},
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if sample id cannot be found} );
like(
    $trap->stderr,
    qr/Could \s+ not \s+ detect/xms,
    q{Throw fatal log message if sample id cannot be found}
);

done_testing();
