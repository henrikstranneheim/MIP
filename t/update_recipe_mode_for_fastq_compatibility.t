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
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log test_mip_hashes };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Analysis} => [qw{ update_recipe_mode_for_fastq_compatibility }],
        q{MIP::Test::Fixtures}   => [qw{ test_log test_mip_hashes }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Analysis qw{ update_recipe_mode_for_fastq_compatibility };
use MIP::Test::Fixtures qw{ test_log test_mip_hashes };

diag(   q{Test update_recipe_mode_for_fastq_compatibility from Analysis.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my %dependency_tree = test_mip_hashes(
    {
        mip_hash_name => q{dependency_tree_rna},
    }
);
my %parameter = test_mip_hashes(
    {
        mip_hash_name => q{define_parameter},
    }
);
my %active_parameter = test_mip_hashes(
    {
        mip_hash_name => q{active_parameter},
    }
);

$parameter{dependency_tree_href} = \%dependency_tree;

my $log = test_log( {} );

my %file_info         = test_mip_hashes( { mip_hash_name => q{file_info}, } );
my $mip_file_format_2 = q{ADM1059A1_161011_HHJJCCCXY_NAATGCGC_lane1};
my $sample_id         = q{ADM1059A1};

## Given that both lanes have been sequenced the same way
my $recipe_fastq_compatibility = update_recipe_mode_for_fastq_compatibility(
    {
        active_parameter_href => \%active_parameter,
        file_info_href        => \%file_info,
        parameter_href        => \%parameter,
        recipe_name           => q{salmon_quant},
    }
);
## Then return true
ok( $recipe_fastq_compatibility, q{Compatible} );

## Given a lane difference
push @{ $file_info{$sample_id}{no_direction_infile_prefixes} }, $mip_file_format_2;
$file_info{$sample_id}{$mip_file_format_2}{sequence_run_type} = q{paired-end};

trap {
    $recipe_fastq_compatibility = update_recipe_mode_for_fastq_compatibility(
        {
            active_parameter_href => \%active_parameter,
            file_info_href        => \%file_info,
            parameter_href        => \%parameter,
            recipe_name           => q{salmon_quant},
        }
    )
};
## Then not compatible
is( $recipe_fastq_compatibility, 0, q{Identify non compatible sequence types} );
like( $trap->stderr, qr/Multiple\ssequence\srun\stypes\sdetected/xms, q{Log warning for non compatible sequence types with recipe} );

done_testing();
