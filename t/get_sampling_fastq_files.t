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
        q{MIP::File_info}      => [qw{ get_sampling_fastq_files }],
        q{MIP::Test::Fixtures} => [qw{ test_mip_hashes }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File_info qw{ get_sampling_fastq_files };

diag(   q{Test get_sampling_fastq_files from File_info.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a sample id
my $sample_id = q{ADM1059A3};

## Given no files to sample from
my %file_info;

## When getting sample of fastq file
my ( $return, ) = get_sampling_fastq_files(
    {
        file_info_sample_href => \%{ $file_info{$sample_id} },
    }
);

## Then return undef
is( $return, undef, q{No file to sample from} );

## Given an interleaved sequence_run_type
%file_info = test_mip_hashes(
    {
        mip_hash_name => q{file_info},
    }
);

push @{ $file_info{$sample_id}{no_direction_infile_prefixes} }, q{ADM1059A3};
$file_info{$sample_id}{ADM1059A3}{sequence_run_type} = q{interleaved};

## When getting sample of fastq file
my ( $is_interleaved_fastq, @fastq_files ) = get_sampling_fastq_files(
    {
        file_info_sample_href => \%{ $file_info{$sample_id} },
    }
);
my @expected_fastq_files = qw{ ADM1059A3.fastq };

## Then return true for interleaved
is( $is_interleaved_fastq, 1, q{Is interleaved file} );

## Then return fastq file
is_deeply( \@fastq_files, \@expected_fastq_files, q{Got interleaved fastq file} );

## Given a paired-end sequence_run_type
push @{ $file_info{$sample_id}{mip_infiles} }, q{ADM1059A3.fastq};
$file_info{$sample_id}{ADM1059A3}{sequence_run_type} = q{paired-end};

## When getting sample of fastq file
( $is_interleaved_fastq, @fastq_files ) = get_sampling_fastq_files(
    {
        file_info_sample_href => \%{ $file_info{$sample_id} },
    }
);
push @expected_fastq_files, q{ADM1059A3.fastq};

## Then return undef for interleaved
is( $is_interleaved_fastq, 0, q{No interleaved files} );

## Then return fastq files
is_deeply( \@fastq_files, \@expected_fastq_files, q{Got fastq files} );

done_testing();
