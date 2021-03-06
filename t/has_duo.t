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
use MIP::Test::Fixtures qw{ test_mip_hashes };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Pedigree}       => [qw{ has_duo }],
        q{MIP::Test::Fixtures} => [qw{ test_mip_hashes }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Pedigree qw{ has_duo };

diag(   q{Test has_duo from Pedigre.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my %active_parameter = test_mip_hashes( { mip_hash_name => q{active_parameter}, } );
my %sample_info      = test_mip_hashes( { mip_hash_name => q{qc_sample_info}, } );

## Given a trio where a sample id has phenotype affected

## When checking if case constellation is duo and has samples with phenotype "affected"
my $has_duo = has_duo(
    {
        active_parameter_href => \%active_parameter,
        sample_info_href      => \%sample_info,
    }
);

## Then return 0
is( $has_duo, 0, q{Detected non duo constellation for trio} );

## Given a duo where a sample id has phenotype affected
@{ $active_parameter{sample_ids} } = (qw{ADM1059A1 ADM1059A2});

## When checking if case constellation is duo and has samples with phenotype "affected"
$has_duo = has_duo(
    {
        active_parameter_href => \%active_parameter,
        sample_info_href      => \%sample_info,
    }
);

## Then return 1
ok( $has_duo, q{Detected duo constellation} );

$sample_info{sample}{ADM1059A1}{phenotype} = q{unaffected};

## When checking if case constellation is duo and has samples with phenotype "affected"
$has_duo = has_duo(
    {
        active_parameter_href => \%active_parameter,
        sample_info_href      => \%sample_info,
    }
);

## Then return 0
is( $has_duo, 0, q{Detected no samples with phenotype affected} );

## Given a single sample case
@{ $active_parameter{sample_ids} } = (q{ADM1059A1});

## When checking if case constellation is duo and has samples with phenotype "affected"
$has_duo = has_duo(
    {
        active_parameter_href => \%active_parameter,
        sample_info_href      => \%sample_info,
    }
);
## Then return 0
is( $has_duo, 0, q{Detected non duo constellation for single sample} );

done_testing();
