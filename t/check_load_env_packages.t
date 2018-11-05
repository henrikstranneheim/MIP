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
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Test::Fixtures qw{ test_mip_hashes test_standard_cli };

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
        q{MIP::Check::Parameter} => [qw{ check_load_env_packages }],
        q{MIP::Test::Fixtures}   => [qw{ test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Check::Parameter qw{ check_load_env_packages };

diag(   q{Test check_load_env_packages from Parameter.pm v}
      . $MIP::Check::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given load_env packages with existing program executables
my %parameter = test_mip_hashes( { mip_hash_name => q{define_parameter}, } );
push @{ $parameter{dynamic_parameter}{program_name_path} }, qw{ bwa samtools sambamba };

my %active_parameter = test_mip_hashes( { mip_hash_name => q{active_parameter}, } );

my $is_ok = check_load_env_packages(
    {
        active_parameter_href => \%active_parameter,
        parameter_href        => \%parameter,
    }
);

## Then all packages should be found
ok( $is_ok, q{Found all packages} );

## Given a not existing recipe
$active_parameter{load_env}{test}{not_a_package} = undef;

trap {
    check_load_env_packages(
        {
            active_parameter_href => \%active_parameter,
            parameter_href        => \%parameter,
        }
    );
};

## Then croak and exist
is( $trap->leaveby, q{die}, q{Exit if the package cannot be found} );
like(
    $trap->die,
    qr/Could\s+not\s+find\s+load_env\spackage/xms,
    q{Throw error if package cannot be found}
);

done_testing();
