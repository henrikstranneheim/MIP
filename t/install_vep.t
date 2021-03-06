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
use autodie qw { :all };
use Modern::Perl qw{ 2018 };
use Readonly;
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COLON $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_constants test_log test_mip_hashes };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Recipes::Install::Vep} => [qw{ install_vep }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_mip_hashes }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Install::Vep qw{ install_vep };

diag(   q{Test install_vep from Vep.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

test_log( { no_screen => 1, } );

test_constants( {} );

## Given install parameters
my %active_parameter =
  test_mip_hashes( { mip_hash_name => q{install_active_parameter}, } );
$active_parameter{reference_dir} = catdir(qw{ reference dir });
my $is_ok = install_vep(
    {
        active_parameter_href => \%active_parameter,
        container_href        => $active_parameter{container}{vep},
    }
);

## Then return TRUE
ok( $is_ok, q{Executed install vep recipe } );

## Given no auto flag 'p'
$active_parameter{vep_auto_flag} = q{acf};
$is_ok = install_vep(
    {
        active_parameter_href => \%active_parameter,
        container_href        => $active_parameter{container}{vep},
    }
);

## Then return TRUE
ok( $is_ok, q{Catch case with no plugins} );

## Given no reference or cache directory
$active_parameter{reference_dir} = undef;
trap {
    install_vep(
        {
            active_parameter_href => \%active_parameter,
            container_href        => $active_parameter{container}{vep},
        }
    )
};

## Then exit and print fatal message
ok( $trap->exit, q{Exit without a reference or cache directory} );

## Given auto flag 'a'
$active_parameter{vep_auto_flag} = q{a};
$is_ok = install_vep(
    {
        active_parameter_href => \%active_parameter,
        container_href        => $active_parameter{container}{vep},
    }
);

## Then return undef
is( $is_ok, undef, q{Return on auto flag a} );

## Given a failure during installation
my %process_return = (
    buffers_ref   => [],
    error_message => q{Error message},
    stderrs_ref   => [],
    stdouts_ref   => [],
    success       => 0,
);
test_constants( { test_process_return_href => \%process_return } );
$active_parameter{reference_dir} = catdir(qw{ reference dir });
$active_parameter{vep_auto_flag} = q{ac};

trap {
    install_vep(
        {
            active_parameter_href => \%active_parameter,
            container_href        => $active_parameter{container}{vep},
        }
    );
};

## Then die
is( $trap->leaveby, q{die}, q{Die on failure} );

done_testing();
