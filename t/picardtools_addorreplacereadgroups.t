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


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Program::Picardtools} => [qw{ picardtools_addorreplacereadgroups }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Picardtools qw{ picardtools_addorreplacereadgroups };

diag(   q{Test picardtools_addorreplacereadgroups from Picardtools.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ picard AddOrReplaceReadGroups };

my %base_argument = (
    filehandle => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
);

## Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument
my %required_argument = (
    infile_path => {
        input           => catfile(qw{ dir test.bam }),
        expected_output => q{-INPUT} . $SPACE . catfile(qw{ dir test.bam }),
    },
    outfile_path => {
        input           => catfile(qw{ dir test.rg.bam }),
        expected_output => q{-OUTPUT} . $SPACE . catfile(qw{ dir test.rg.bam }),
    },
    readgroup_id => {
        input           => 1,
        expected_output => q{-RGID 1},
    },
    readgroup_library => {
        input           => q{lib},
        expected_output => q{-RGLB lib},
    },
    readgroup_platform => {
        input           => q{IlluminaHiseqXten},
        expected_output => q{-RGPL IlluminaHiseqXten},
    },
    readgroup_platform_unit => {
        input           => q{UnitX},
        expected_output => q{-RGPU UnitX},
    },
    readgroup_sample => {
        input           => q{Sven},
        expected_output => q{-RGSM Sven},
    },

);

my %specific_argument = (
    infile_path => {
        input           => catfile(qw{ dir test.bam }),
        expected_output => q{-INPUT} . $SPACE . catfile(qw{ dir test.bam }),
    },
    outfile_path => {
        input           => catfile(qw{ dir test.rg.bam }),
        expected_output => q{-OUTPUT} . $SPACE . catfile(qw{ dir test.rg.bam }),
    },
    readgroup_id => {
        input           => 1,
        expected_output => q{-RGID 1},
    },
    readgroup_library => {
        input           => q{lib},
        expected_output => q{-RGLB lib},
    },
    readgroup_platform => {
        input           => q{IlluminaHiseqXten},
        expected_output => q{-RGPL IlluminaHiseqXten},
    },
    readgroup_platform_unit => {
        input           => q{UnitX},
        expected_output => q{-RGPU UnitX},
    },
    readgroup_sample => {
        input           => q{Sven},
        expected_output => q{-RGSM Sven},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&picardtools_addorreplacereadgroups;

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

## Base arguments
@function_base_commands = qw{ picard java };

my %specific_java_argument = (
    java_jar => {
        input           => q{gatk.jar},
        expected_output => q{-jar gatk.jar},
    },
);

## Test both base and function specific arguments
@arguments = ( \%specific_java_argument );

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
