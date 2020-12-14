#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname  };
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
        q{MIP::Program::Gatk}  => [qw{ gatk_leftalignandtrimvariants }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Gatk qw{ gatk_leftalignandtrimvariants };

diag(   q{Test gatk_leftalignandtrimvariants from Gatk.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ gatk LeftAlignAndTrimVariants };

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
        input           => catfile(qw{ path_to_analysis_dir infile.vcf }),
        expected_output => q{--variant}
          . $SPACE
          . catfile(qw{ path_to_analysis_dir infile.vcf }),
    },
    outfile_path => {
        input           => catfile(qw{ path_to_analysis_dir outfile.vcf }),
        expected_output => q{--output}
          . $SPACE
          . catfile(qw{ path_to_analysis_dir outfile.vcf }),
    },
    referencefile_path => {
        input           => catfile(qw{ my genome }),
        expected_output => q{--reference } . catdir(qw{ my genome }),
    },
);

my %specific_argument = (
    infile_path => {
        input           => catfile(qw{ path_to_analysis_dir infile.vcf }),
        expected_output => q{--variant}
          . $SPACE
          . catfile(qw{ path_to_analysis_dir infile.vcf }),
    },
    outfile_path => {
        input           => catfile(qw{ path_to_analysis_dir outfile.vcf }),
        expected_output => q{--output}
          . $SPACE
          . catfile(qw{ path_to_analysis_dir outfile.vcf }),
    },
    referencefile_path => {
        input           => catfile(qw{ my genome }),
        expected_output => q{--reference } . catdir(qw{ my genome }),
    },
    split_multiallelics => {
        input           => 1,
        expected_output => q{--split-multi-allelics},
    },

);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&gatk_leftalignandtrimvariants;

## Test both base and function specific arguments
my @arguments = ( \%base_argument, \%specific_argument );

ARGUMENT_HASH_REF:
foreach my $argument_href (@arguments) {
    my @commands = test_function(
        {
            argument_href              => $argument_href,
            base_commands_index        => 1,
            do_test_base_command       => 1,
            function_base_commands_ref => \@function_base_commands,
            module_function_cref       => $module_function_cref,
            required_argument_href     => \%required_argument,
        }
    );
}

done_testing();
