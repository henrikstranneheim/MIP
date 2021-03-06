#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Path qw{ rmtree };
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
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Environment::Child_process} => [qw{child_process}],
        q{MIP::Test::Fixtures}             => [qw{test_log}],
        q{MIP::Toml}                       => [qw{ load_toml write_toml }],
        q{MIP::Vcfanno}                    => [qw{check_vcfanno_toml}],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Environment::Child_process qw{ child_process };
use MIP::Toml qw{ load_toml write_toml };
use MIP::Vcfanno qw{ check_vcfanno_toml };

diag(   q{Test check_vcfanno_toml from Vcfanno.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Creates log object
my $log = test_log( {} );

## Replace file path depending on location - required for TRAVIS
my $test_reference_dir = catfile( $Bin, qw{ data references } );

### Prepare temporary file for testing
my $vcfanno_config =
  catfile( $test_reference_dir,
    qw{ grch37_frequency_vcfanno_filter_config_-v1.0-.toml  } );

# For the actual test
my $test_vcfanno_config = catfile( $test_reference_dir,
    qw{ grch37_frequency_vcfanno_filter_config_test_check_vcfanno_-v1.0-.toml  } );

my $toml_href = load_toml(
    {
        path => $vcfanno_config,
    }
);

## Set test file paths
$toml_href->{functions}{file} =
  catfile( $Bin, qw{ data references vcfanno_functions_-v1.0-.lua } );
$toml_href->{annotation}[0]{file} =
  catfile( $Bin, qw{ data references grch37_gnomad.genomes_-r2.0.1-.vcf.gz } );
$toml_href->{annotation}[1]{file} =
  catfile( $Bin, qw{ data references grch37_gnomad_reformated_-r2.1.1_sv-.vcf.gz } );
$toml_href->{annotation}[2]{file} =
  catfile( $Bin, qw{ data references grch37_cadd_whole_genome_snvs_-v1.4-.tsv.gz } );

write_toml(
    {
        data_href => $toml_href,
        path      => $test_vcfanno_config,
    }
);

my %active_parameter = ( vcfanno_config => $test_vcfanno_config, );

## Given a toml config file with a file path
my $is_ok = check_vcfanno_toml(
    {
        active_parameter_href => \%active_parameter,
        vcfanno_config_name   => q{vcfanno_config},
        vcfanno_functions     => q{vcfanno_functions},
    }
);

## Then return true
ok( $is_ok, q{Passed check for toml file} );

## Clean-up
rmtree($test_vcfanno_config);

## Given a toml config file, when mandatory features are absent
my $faulty_vcfanno_config_file = catfile( $Bin,
    qw{ data references grch37_frequency_vcfanno_filter_config_bad_data_-v1.0-.toml } );

$active_parameter{vcfanno_config} = $faulty_vcfanno_config_file;
trap {
    check_vcfanno_toml(
        {
            active_parameter_href => \%active_parameter,
            vcfanno_config_name   => q{vcfanno_config},
            vcfanno_functions     => q{vcfanno_functions},
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if the record does not match} );
like( $trap->stderr, qr/FATAL/xms,
    q{Throw fatal log message for non matching reference} );

done_testing();
