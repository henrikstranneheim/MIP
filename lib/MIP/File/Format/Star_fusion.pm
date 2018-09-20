package MIP::File::Format::Star_fusion;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ create_star_fusion_sample_file };
}

## Constants
Readonly my $NEWLINE => q{\n};
Readonly my $SPACE   => q{ };
Readonly my $TAB     => q{\t};

sub create_star_fusion_sample_file {

## Function : Create the samples file for STAR-fusion.
## Returns  :
## Arguments: $FILEHANDLE              => Filehandle to write to
##          : $infiles_ref             => Infiles for sample {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $insample_directory      => Directory for fastq files
##          : $sample_id               => Sample id
##          : $samples_file_path       => The family file path
##          : $sample_info_href        => Info on samples and family hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infiles_ref;
    my $infile_lane_prefix_href;
    my $insample_directory;
    my $sample_id;
    my $samples_file_path;
    my $sample_info_href;

    my $tmpl = {
        FILEHANDLE => {
            defined  => 1,
            required => 1,
            store    => \$FILEHANDLE,
        },
        infiles_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$infiles_ref,
            strict_type => 1,
        },
        infile_lane_prefix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_lane_prefix_href,
            strict_type => 1,
        },
        insample_directory => {
            defined     => 1,
            required    => 1,
            store       => \$insample_directory,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
        samples_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$samples_file_path,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Gnu::Coreutils qw{ gnu_echo };

    say {$FILEHANDLE} q{# Generating STAR-fusion 'samples_file'};

    my %sample_line;
    my @strings;

    # Too adjust infile_index for paired-ends
    my $paired_end_tracker = 0;

    ## Perform per single-end or read pair
  INFILE_PREFIX:
    while ( my ( $infile_index, $infile_prefix ) =
        each @{ $infile_lane_prefix_href->{$sample_id} } )
    {

        ## Add sample id file index array
        push @{ $sample_line{$infile_index} }, $sample_id;

        # Collect paired-end or single-end sequence run mode
        my $sequence_run_mode =
          $sample_info_href->{sample}{$sample_id}{file}{$infile_prefix}
          {sequence_run_type};

        my $insample_dir_fastqc_path_read_one =
          catfile( $insample_directory, $infiles_ref->[$paired_end_tracker] );

        ## Add read one to file index array
        push @{ $sample_line{$infile_index} },
          $insample_dir_fastqc_path_read_one;

        # If second read direction is present
        if ( $sequence_run_mode eq q{paired-end} ) {

            # Increment to collect correct read 2 from %infile
            $paired_end_tracker = $paired_end_tracker + 1;
            my $insample_dir_fastqc_path_read_two =
              catfile( $insample_directory,
                $infiles_ref->[$paired_end_tracker] );

            ## Add read two file index array
            push @{ $sample_line{$infile_index} },
              $insample_dir_fastqc_path_read_two;

        }

        ## Increment paired end tracker
        $paired_end_tracker++;

        ## Add tab to each element and add as string to array
        push @strings, join $TAB, @{ $sample_line{$infile_index} };

    }

    ## Add newline to each string (line)
    @strings = map { $_ . $NEWLINE } @strings;

    gnu_echo(
        {
            enable_interpretation => 1,
            FILEHANDLE            => $FILEHANDLE,
            no_trailing_newline   => 1,
            outfile_path          => $samples_file_path,
            strings_ref           => \@strings,
        }
    );
    return 1;
}

1;
