#! /usr/bin/env perl

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
use MIP::Constants qw{ $COMMA $SPACE $UNDERSCORE };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module =
      ( q{MIP::Processmanagement::Processes} => [qw{ clear_case_id_job_id_dependency_tree }], );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Processmanagement::Processes qw{ clear_case_id_job_id_dependency_tree };

diag(   q{Test clear_case_id_job_id_dependency_tree from Processes.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my $case_id             = q{case1};
my $sample_id           = q{sample1};
my $path                = q{MAIN};
my $case_id_chain_key   = $case_id . $UNDERSCORE . $path;
my $sample_id_chain_key = $sample_id . $UNDERSCORE . $path;
my $pan_chain_key       = $case_id_chain_key . $UNDERSCORE . $sample_id_chain_key;

my %job_id = (
    $case_id_chain_key => {
        q{sample1} . $UNDERSCORE . $path => [qw{job_id_1 job_id_2}],
        q{sample2} . $UNDERSCORE . $path => [qw{job_id_3}],
        q{sample3} . $UNDERSCORE . $path => [qw{job_id_4 job_id_5 job_id_8}],
        q{sample4} . $UNDERSCORE . $path => [undef],
        $pan_chain_key                   => [qw{job_id_1 job_id_2}],
        $case_id_chain_key               => [qw{job_id_6}],
    },
);

### Clear pan job id dependency for case id

## Add job_ids from MAIN chain to job_id_string
clear_case_id_job_id_dependency_tree(
    {
        case_id_chain_key => $case_id_chain_key,
        job_id_href       => \%job_id,
    }
);

my $result_ref      = $job_id{$case_id_chain_key}{$case_id_chain_key};
my $expected_result = 0;
is( @{$result_ref}, $expected_result, q{Cleared job_ids from case_id job_id chain} );

done_testing();
