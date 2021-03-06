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
use autodie qw{ :all };
use Modern::Perl qw{ 2018 };

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $SPACE $COMMA };
use MIP::Test::Commands qw{ test_function };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Program::Bcftools} => [qw{ bcftools_query }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Bcftools qw{ bcftools_query };

diag(   q{Test bcftools_query from Bcftools.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ bcftools query };

my %base_argument = (
    filehandle => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
    stderrfile_path => {
        input           => q{stderrfile.test},
        expected_output => q{2> stderrfile.test},
    },
    stderrfile_path_append => {
        input           => q{stderrfile.test},
        expected_output => q{2>> stderrfile.test},
    },
    stdoutfile_path => {
        input           => q{stdoutfile.test},
        expected_output => q{1> stdoutfile.test},
    },
);

## Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument
my %specific_argument = (
    exclude => {
        input           => q{%QUAL<10 || (RPB<0.1 && %QUAL<15)},
        expected_output => q{--exclude %QUAL<10 || (RPB<0.1 && %QUAL<15)},
    },
    format => {
        input           => q{'%CHROM\t%POS\t%REF,%ALT\t%INFO/TAG\n'},
        expected_output => q{--format '%CHROM\t%POS\t%REF,%ALT\t%INFO/TAG\n'},
    },
    include => {
        input           => q{INFO/CSQ[*]~":p[.]"},
        expected_output => q{--include INFO/CSQ[*]~":p[.]"},
    },
    infile_paths_ref => {
        inputs_ref      => [qw{ infile.test }],
        expected_output => q{infile.test},
    },
    outfile_path => {
        input           => q{outfile.txt},
        expected_output => q{--output outfile.txt},
    },
    print_header => {
        input           => 1,
        expected_output => q{--print-header},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&bcftools_query;

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
        }
    );
}

done_testing();
