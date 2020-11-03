#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Modern::Perl qw{ 2018 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Commands qw{ test_function };
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.04;

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
        q{MIP::Program::Bcftools} => [qw{ bcftools_base }],
        q{MIP::Test::Fixtures}    => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Bcftools qw{ bcftools_base };

diag(   q{Test bcftools_base from Bcftools.pm v}
      . $MIP::Program::Bcftools::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $NR_THREADS_TO_USE => 12;

## Base arguments
my @function_base_commands = qw{ bcftools };

my %base_argument = (
    filehandle => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
);

## Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument
my %required_argument = (
    commands_ref => {
        inputs_ref      => [qw{ bcftools mpileup }],
        expected_output => q{bcftools mpileup},
    },
    filehandle => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
);

my %specific_argument = (
    commands_ref => {
        inputs_ref      => [qw{ bcftools mpileup }],
        expected_output => q{mpileup},
    },
    outfile_path => {
        input           => catfile(qw{ a test file }),
        expected_output => q{--output} . $SPACE . catfile(qw{ a test file }),
    },
    output_type => {
        input           => q{b},
        expected_output => q{--output-type} . $SPACE . q{b},
    },
    regions_file_path => {
        input           => catfile(qw{ a test.bed }),
        expected_output => q{--regions-file} . $SPACE . catfile(qw{ a test.bed}),
    },
    regions_ref => {
        inputs_ref      => [qw{ 1 2 }],
        expected_output => q{--regions 1,2},
    },
    samples_file_path => {
        input           => catfile(qw{ a test sample_file }),
        expected_output => q{--samples-file} . $SPACE . catfile(qw{ a test sample_file }),
    },
    samples_ref => {
        inputs_ref      => [qw{ ^sample_1 sample_2 }],
        expected_output => q{--samples ^sample_1,sample_2},
    },
    threads => {
        input           => $NR_THREADS_TO_USE,
        expected_output => q{--threads} . $SPACE . $NR_THREADS_TO_USE,
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&bcftools_base;

## Test both base and function specific arguments
my @arguments = ( \%base_argument, \%specific_argument );

ARGUMENT_HASH_REF:
foreach my $argument_href (@arguments) {
    my @commands = test_function(
        {
            argument_href              => $argument_href,
            do_test_base_command       => 1,
            function_base_commands_ref => \@function_base_commands,
            module_function_cref       => $module_function_cref,
            required_argument_href     => \%required_argument,
        }
    );
}

done_testing();
