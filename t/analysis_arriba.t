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
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COLON $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_add_io_for_recipe test_log test_mip_hashes test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.02;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Recipes::Analysis::Arriba} => [qw{ analysis_arriba }],
        q{MIP::Test::Fixtures} => [qw{ test_add_io_for_recipe test_log test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Analysis::Arriba qw{ analysis_arriba };

diag(   q{Test analysis_arriba from Arriba.pm v}
      . $MIP::Recipes::Analysis::Arriba::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( { log_name => q{MIP}, no_screen => 1, } );

## Given analysis parameters
my $recipe_name    = q{arriba_ar};
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
$active_parameter{star_aln_reference_genome}        = q{a_dir};
$active_parameter{transcript_annotation}            = q{transcripts.gtf};
$active_parameter{arriba_blacklist_path}            = q{blacklist.tsv};
$active_parameter{arriba_cytoband_path}             = q{cytobands.tsv};
$active_parameter{arriba_proteindomain_path}        = q{proteindomains.gff};
$active_parameter{platform}                         = q{ILLUMINA};
my $sample_id = $active_parameter{sample_ids}[0];

my %file_info = test_mip_hashes(
    {
        mip_hash_name => q{file_info},
        recipe_name   => $recipe_name,
    }
);
@{ $file_info{$sample_id}{lanes} } = ( 1, 2 );
$file_info{star_aln_reference_genome} = [q{reference_genome}];
$file_info{$sample_id}{$recipe_name}{file_tag} = q{trim};

my %job_id;
my %parameter = test_mip_hashes(
    {
        mip_hash_name => q{recipe_parameter},
        recipe_name   => $recipe_name,
    }
);
test_add_io_for_recipe(
    {
        file_info_href => \%file_info,
        id             => $sample_id,
        outfile_suffix => q{.tsv},
        parameter_href => \%parameter,
        recipe_name    => $recipe_name,
        step           => q{bam},
    }
);

my %sample_info = (
    sample => {
        $sample_id => {
            file => {
                ADM1059A1_161011_HHJJCCCXY_NAATGCGC_lane7 => {
                    sequence_run_type   => q{single-end},
                    read_direction_file => {
                        ADM1059A1_161011_HHJJCCCXY_NAATGCGC_lane7_1 => {
                            flowcell       => q{HHJJCCCXY},
                            lane           => q{7},
                            sample_barcode => q{NAATGCGC},
                        },
                    },
                },
            },
        },
    },
);

my $is_ok = analysis_arriba(
    {
        active_parameter_href => \%active_parameter,
        file_info_href        => \%file_info,
        job_id_href           => \%job_id,
        parameter_href        => \%parameter,
        profile_base_command  => $slurm_mock_cmd,
        recipe_name           => $recipe_name,
        sample_id             => $sample_id,
        sample_info_href      => \%sample_info,
    }
);

## Then return TRUE
ok( $is_ok, q{ Executed analysis recipe } . $recipe_name );

done_testing();
