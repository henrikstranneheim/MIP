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
        q{MIP::Sample_info}    => [qw{ get_case_members_attributes_in_duos }],
        q{MIP::Test::Fixtures} => [qw{ test_mip_hashes }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Sample_info qw{ get_case_members_attributes_in_duos };

diag(   q{Test get_case_members_attributes_in_duos from Sample_info.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given sample_info hash
my %sample_info = test_mip_hashes( { mip_hash_name => q{qc_sample_info}, } );
$sample_info{sample}{ADM1059A2}{phenotype} = q{unknown};
$sample_info{sample}{ADM1059A1}{mother}    = q{0};

delete $sample_info{sample}{ADM1059A3};
my %family_member_id =
  get_case_members_attributes_in_duos( { sample_info_href => \%sample_info, } );

## Then return family members id
my %expected = (
    affected => [qw{ ADM1059A1 }],
    children => [qw{ ADM1059A1 }],
    father   => q{ADM1059A2},
    mother   => q{0},
    unknown  => [qw{ ADM1059A2 }],
);

is_deeply( \%family_member_id, \%expected, q{Got family hash} );

done_testing();
