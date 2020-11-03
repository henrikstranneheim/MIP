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
use MIP::Test::Fixtures qw{ test_add_io_for_recipe test_log test_mip_hashes test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.05;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $GENOME_BUILD_VERSION_HG   => 20;
Readonly my $GENOME_BUILD_VERSION_GRCH => 37;

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Recipes::Analysis::Plink} => [qw{ analysis_plink }],
        q{MIP::Test::Fixtures} => [qw{ test_add_io_for_recipe test_log test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Analysis::Plink qw{ analysis_plink };

diag(   q{Test analysis_plink from Plink.pm v}
      . $MIP::Recipes::Analysis::Plink::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( { no_screen => 1 } );

## Given analysis parameters
my $recipe_name    = q{plink};
my $slurm_mock_cmd = catfile( $Bin, qw{ data modules slurm-mock.pl } );

my %active_parameter = test_mip_hashes(
    {
        mip_hash_name => q{active_parameter},
        recipe_name   => $recipe_name,
    }
);
$active_parameter{$recipe_name}                     = 1;
$active_parameter{recipe_core_number}{$recipe_name} = 1;
$active_parameter{recipe_time}{$recipe_name}        = 1;
my $case_id = $active_parameter{case_id};
$active_parameter{gender}{others} = [qw{ ADM1059A3 }];
$active_parameter{exome_target_bed}{q{target_capture.bed}} = q{ADM1059A3};

my %file_info = test_mip_hashes(
    {
        mip_hash_name => q{file_info},
        recipe_name   => $recipe_name,
    }
);

$file_info{human_genome_reference_version} = $GENOME_BUILD_VERSION_HG;
$file_info{human_genome_reference_source}  = q{hg};

my %job_id;
my %parameter = test_mip_hashes(
    {
        mip_hash_name => q{recipe_parameter},
        recipe_name   => $recipe_name,
    }
);

test_add_io_for_recipe(
    {
        file_info_href    => \%file_info,
        id                => $case_id,
        parameter_href    => \%parameter,
        recipe_name       => $recipe_name,
        step              => q{vcf},
    }
);

$parameter{cache}{consensus_analysis_type} = q{wes};

my %sample_info = test_mip_hashes(
    {
        mip_hash_name => q{qc_sample_info},
        recipe_name   => $recipe_name,
    }
);

my $is_ok = analysis_plink(
    {
        active_parameter_href => \%active_parameter,
        case_id               => $case_id,
        file_info_href        => \%file_info,
        job_id_href           => \%job_id,
        parameter_href        => \%parameter,
        profile_base_command  => $slurm_mock_cmd,
        recipe_name           => $recipe_name,
        sample_info_href      => \%sample_info,
    }
);

## Then return TRUE
ok( $is_ok, q{Executed analysis recipe } . $recipe_name );

## Given a grch reference
$file_info{human_genome_reference_version} = $GENOME_BUILD_VERSION_GRCH;
$file_info{human_genome_reference_source}  = q{grch};
$is_ok                                     = analysis_plink(
    {
        active_parameter_href => \%active_parameter,
        case_id               => $case_id,
        file_info_href        => \%file_info,
        job_id_href           => \%job_id,
        parameter_href        => \%parameter,
        profile_base_command  => $slurm_mock_cmd,
        recipe_name           => $recipe_name,
        sample_info_href      => \%sample_info,
    }
);

## Then return TRUE
ok( $is_ok, q{Executed analysis recipe } . $recipe_name );

## Given no eliglbe test to run
# Reduce to single sample
@{ $active_parameter{sample_ids} } = $active_parameter{sample_ids}[0];

my $return = analysis_plink(
    {
        active_parameter_href => \%active_parameter,
        case_id               => $case_id,
        file_info_href        => \%file_info,
        job_id_href           => \%job_id,
        parameter_href        => \%parameter,
        profile_base_command  => $slurm_mock_cmd,
        recipe_name           => $recipe_name,
        sample_info_href      => \%sample_info,
    }
);

## Then return undef, since there are no tests to run
is( $return, undef, q{ Skipped analysis recipe } . $recipe_name );

done_testing();
