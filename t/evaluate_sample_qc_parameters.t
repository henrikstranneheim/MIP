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
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COLON $COMMA $SPACE $UNDERSCORE};
use MIP::Test::Fixtures qw{ test_mip_hashes };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Qccollect}      => [qw{ evaluate_sample_qc_parameters }],
        q{MIP::Test::Fixtures} => [qw{ test_mip_hashes }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Qccollect qw{ evaluate_sample_qc_parameters };

diag(   q{Test evaluate_sample_qc_parameters from Qccollect.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $PERCENTAGE_MAPPED_READS_EVAL => 95;
Readonly my $PERCENTAGE_MAPPED_READS_PASS => 99;

## Given sample qc data when metric lacks a header
my $infile    = q{an_infile};
my $metric    = q{percentage_mapped_reads};
my $recipe    = q{bamstats};
my $sample_id = q{ADM1059A1};
my %qc_data   = (
    sample => {
        $sample_id =>
          { $infile => { $recipe => { $metric => $PERCENTAGE_MAPPED_READS_PASS }, }, },
    }
);

## Add test for skipping evaluation key
$qc_data{sample}{$sample_id}{an_EvalUation_infile} = q{PASS};

## Add test for skipping infile names with undetermined key
$qc_data{sample}{$sample_id}{an_Undetermined_infile} = q{who knows};

## Add dummy entry for plink relation checking
$qc_data{sample}{$sample_id}{plink_relation_check} = q{PASS};

## Alias
my $qc_data_recipe_href =
  \%{ $qc_data{sample}{$sample_id}{$infile}{$recipe} };

my %sample_info = test_mip_hashes(
    {
        mip_hash_name => q{qc_sample_info},
        recipe_name   => $recipe,
    }
);

## Defines recipes, metrics and thresholds to evaluate
my %evaluate_metric = (
    $sample_id => {
        $infile => {
            $recipe => {
                $metric => $PERCENTAGE_MAPPED_READS_EVAL,
            },
        },
    },
);

my $is_ok = evaluate_sample_qc_parameters(
    {
        evaluate_metric_href => \%evaluate_metric,
        qc_data_href         => \%qc_data,
    }
);

ok( $is_ok, q{Evaluated sample qc data metrics} );

done_testing();
