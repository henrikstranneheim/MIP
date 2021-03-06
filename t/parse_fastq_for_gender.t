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
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log test_mip_hashes };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Recipes::Analysis::Estimate_gender} => [qw{ parse_fastq_for_gender }],
        q{MIP::Test::Fixtures}                     => [qw{ test_log test_mip_hashes }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Analysis::Estimate_gender qw{ parse_fastq_for_gender };

diag(   q{Test parse_fastq_for_gender from Estimate_gender.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

test_log( { no_screen => 0, } );

## Given no other gender
my %active_parameter = test_mip_hashes(
    {
        mip_hash_name => q{active_parameter},
        recipe_name   => q{bwa_mem},
    }
);
my $consensus_analysis_type = q{wgs};
my $sample_id               = $active_parameter{sample_ids}[2];
my %file_info               = test_mip_hashes(
    {
        mip_hash_name => q{file_info},
    }
);
push @{ $file_info{$sample_id}{mip_infiles} },                  q{ADM1059A3.fastq};
push @{ $file_info{$sample_id}{no_direction_infile_prefixes} }, q{ADM1059A3};
$file_info{$sample_id}{ADM1059A3}{sequence_run_type} = q{paired-end};
my %sample_info = test_mip_hashes(
    {
        mip_hash_name => q{qc_sample_info},
    }
);

my $estimated_unknown_gender = parse_fastq_for_gender(
    {
        active_parameter_href   => \%active_parameter,
        consensus_analysis_type => $consensus_analysis_type,
        file_info_href          => \%file_info,
        sample_info_href        => \%sample_info,
    }
);

## Then skip estimation using reads
is( $estimated_unknown_gender, undef, q{No unknown gender - skip} );

## Given a sample when gender unknown
$active_parameter{gender}{others} = [$sample_id];

$estimated_unknown_gender = parse_fastq_for_gender(
    {
        active_parameter_href   => \%active_parameter,
        consensus_analysis_type => $consensus_analysis_type,
        file_info_href          => \%file_info,
        sample_info_href        => \%sample_info,
    }
);

## Then skip estimation using reads
is( $estimated_unknown_gender, 1, q{Estimated unknown gender} );

done_testing();
