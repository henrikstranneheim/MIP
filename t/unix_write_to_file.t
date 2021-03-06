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
use Modern::Perl qw{ 2014 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Unix::Write_to_file} => [qw{ unix_write_to_file }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Unix::Write_to_file qw{ unix_write_to_file };
use MIP::Test::Writefile qw{ test_write_to_file };

diag(   q{Test unix_write_to_file from Write_to_file.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Coderef - enables generalized use of generate call
my $module_function_cref = \&unix_write_to_file;

### Tests
## In line test
my @commands = ( q{commands_ref}, [qw{ test in-line }] );
test_write_to_file(
    {
        args_ref             => \@commands,
        base_commands_ref    => [qw{ test in-line }],
        module_function_cref => $module_function_cref,
    }
);

## Per line test and testing with singularity
my %active_parameter = ( pedigree_file => q{a_pedigre_file.yaml}, );

@commands = ( q{commands_ref}, [qw{ test per-line }] );
test_write_to_file(
    {
        args_ref             => \@commands,
        base_commands_ref    => [qw{ test }],
        module_function_cref => $module_function_cref,
        separator            => q{\n},
    }
);

done_testing();
