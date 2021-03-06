#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename dirname };
use File::Spec::Functions qw{ catdir catfile };
use File::Temp;
use FindBin qw{ $Bin };
use Getopt::Long;
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use Test::More;
use Time::Piece;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw { :all };
use Modern::Perl qw{ 2018 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log };

## Constants
Readonly my $READ_LENGTH => 151;

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Sample_info}    => [qw{ set_infile_info }],
        q{MIP::Test::Fixtures} => [qw{ test_log }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Fastq qw{ define_mip_fastq_file_features };
use MIP::File_info qw{ get_sample_file_attribute };
use MIP::Sample_info qw{ set_infile_info };

diag(   q{Test set_infile_info from Sample_info.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Create temp logger
my $log = test_log( {} );

## Set file info parameters
my $fastq_file_read_one = q{file-1};
my $fastq_file_read_two = q{file-2};
my $sample_id           = q{sample-1};

my %file_info = (
    $sample_id => {
        $fastq_file_read_one => {
            date             => q{150703},
            direction        => 1,
            flowcell         => q{Undetermined-flow-rider},
            index            => q{ATCG},
            infile_sample_id => $sample_id,
            is_interleaved   => 0,
            lane             => 1,
            read_length      => $READ_LENGTH,
        },
        $fastq_file_read_two => {
            date             => q{150703},
            direction        => 2,
            flowcell         => q{Undetermined-flow-rider},
            index            => q{ATCG},
            infile_sample_id => $sample_id,
            lane             => 1,
            read_length      => $READ_LENGTH,
        },
        lanes       => [1],
        mip_infiles => [ $fastq_file_read_one, $fastq_file_read_two ],
    },
);
my %sample_info;

my %attribute = get_sample_file_attribute(
    {
        file_info_href => \%file_info,
        file_name      => $fastq_file_read_one,
        sample_id      => $sample_id,
    }
);

## Define file formats
my ( $mip_file_format, $mip_file_format_with_direction,
    $original_file_name_prefix, $run_barcode )
  = define_mip_fastq_file_features(
    {
        date               => $attribute{date},
        direction          => $attribute{direction},
        flowcell           => $attribute{flowcell},
        index              => $attribute{index},
        lane               => $attribute{lane},
        original_file_name => $fastq_file_read_one,
        sample_id          => $sample_id,
    }
  );

my %attribute_2 = get_sample_file_attribute(
    {
        file_info_href => \%file_info,
        file_name      => $fastq_file_read_two,
        sample_id      => $sample_id,
    }
);
my ( $mip_file_format_2, $mip_file_format_with_direction_2 ) =
  define_mip_fastq_file_features(
    {
        date               => $attribute_2{date},
        direction          => $attribute_2{direction},
        flowcell           => $attribute_2{flowcell},
        index              => $attribute_2{index},
        lane               => $attribute_2{lane},
        original_file_name => $fastq_file_read_two,
        sample_id          => $sample_id,
    }
  );

my $parsed_date = Time::Piece->strptime( $attribute{date}, q{%y%m%d} );
$parsed_date = $parsed_date->ymd;
my $lane_tracker = 0;

## Given single file, when Undetermined in flowcell name
SAMPLE_ID:
for my $sample_id ( keys %file_info ) {

    $lane_tracker = 0;

  INFILE:
    foreach my $file_name ( @{ $file_info{$sample_id}{mip_infiles} } ) {

        get_sample_file_attribute(
            {
                file_info_href => \%file_info,
                file_name      => $file_name,
                sample_id      => $sample_id,
            }
        );
        $lane_tracker = set_infile_info(
            {
                file_info_href   => \%file_info,
                file_name        => $file_name,
                lane_tracker     => $lane_tracker,
                sample_id        => $sample_id,
                sample_info_href => \%sample_info,
            }
        );
    }
}

## Collect what to expect in one hash
my %expected_result = (
    lane => {
        $sample_id => {
            lanes => [1],
        },
    },
    sample_info => {
        sample => {
            $sample_id => {
                file => {
                    $mip_file_format => {
                        sequence_run_type   => q{paired-end},
                        sequence_length     => $READ_LENGTH,
                        interleaved         => 0,
                        read_direction_file => {
                            $mip_file_format_with_direction => {
                                date                      => $parsed_date,
                                flowcell                  => $attribute{flowcell},
                                lane                      => $attribute{lane},
                                original_file_name        => $fastq_file_read_one,
                                original_file_name_prefix => $original_file_name_prefix,
                                read_direction            => $attribute{direction},
                                run_barcode               => $run_barcode,
                                sample_barcode            => $attribute{index},
                            },
                            $mip_file_format_with_direction_2 => {
                                date                      => $parsed_date,
                                flowcell                  => $attribute_2{flowcell},
                                lane                      => $attribute_2{lane},
                                original_file_name        => $fastq_file_read_two,
                                original_file_name_prefix => $original_file_name_prefix,
                                sample_barcode            => $attribute_2{index},
                                read_direction            => $attribute_2{direction},
                                run_barcode               => $run_barcode,
                            },
                        },
                    },
                },
            },
        },
    },
);

## Then lane tracker should be one
is( $lane_tracker, 1, q{Tracked lane} );

## Then add the lane info
is_deeply(
    \@{ $file_info{$sample_id}{lanes} },
    \@{ $expected_result{lane}{$sample_id}{lanes} },
    q{Added lane info for single-end read}
);

## Then add no_direction_infile_prefixes to file_info
is( $file_info{$sample_id}{no_direction_infile_prefixes}[0],
    $mip_file_format, q{Added no_direction_infile_prefixes to file_info } );

## Then add no_direction_infile_prefixes with sequence type
is( $file_info{$sample_id}{$mip_file_format}{sequence_run_type},
    q{paired-end},
    q{Added sequence run type to no_direction_infile_prefixes in file_info } );

## Then add single-end read info from file name
is_deeply(
    \%sample_info,
    \%{ $expected_result{sample_info} },
    q{Added sample info for paired-end read}
);

done_testing();
