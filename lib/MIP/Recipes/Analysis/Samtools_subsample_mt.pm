package MIP::Recipes::Analysis::Samtools_subsample_mt;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use List::MoreUtils qw{ first_value };
use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.04;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_samtools_subsample_mt };

}

## Constants
Readonly my $BACKTICK                    => q{`};
Readonly my $DOT                         => q{.};
Readonly my $NEWLINE                     => qq{\n};
Readonly my $MAX_DEPTH_TRESHOLD          => 500_000;
Readonly my $MAX_LIMIT_SEED              => 100;
Readonly my $PIPE                        => q{|};
Readonly my $SAMTOOLS_UNMAPPED_READ_FLAG => 4;
Readonly my $SPACE                       => q{ };
Readonly my $UNDERSCORE                  => q{_};

sub analysis_samtools_subsample_mt {

## Function : Creates a BAM file containing a subset of the MT alignments
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $family_id               => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $program_name            => Program name
##          : $sample_id               => Sample id
##          : $sample_info_href        => Info on samples and family hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $program_name;
    my $sample_id;
    my $sample_info_href;

    ## Default(s)
    my $family_id;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            store       => \$family_id,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        infile_lane_prefix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_lane_prefix_href,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        program_name => {
            defined     => 1,
            required    => 1,
            store       => \$program_name,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
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

    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_module_parameters get_program_attributes };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Program::Alignment::Samtools
      qw{ samtools_depth samtools_index samtools_view };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::QC::Record qw{ add_program_outfile_to_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Unpack parameters
    ## Get the io infiles per chain and id
    my %io = get_io_files(
        {
            id             => $sample_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            program_name   => $program_name,
            stream         => q{in},
        }
    );
    my $infile_name_prefix   = $io{in}{file_name_prefix};
    my @infile_name_prefixes = @{ $io{in}{file_name_prefixes} };
    my @infile_paths         = @{ $io{in}{file_paths} };

    ## Find Mitochondrial contig infile_path
    my $infile_path = first_value { / $infile_name_prefix [.]M /sxm } @infile_paths;

    my $job_id_chain = get_program_attributes(
        {
            parameter_href => $parameter_href,
            program_name   => $program_name,
            attribute      => q{chain},
        }
    );
    my $mt_subsample_depth = $active_parameter_href->{samtools_subsample_mt_depth};
    my $program_mode       = $active_parameter_href->{$program_name};
    my ( $core_number, $time, @source_environment_cmds ) = get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name          => $program_name,
        }
    );

    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id               => $job_id_chain,
                id                     => $sample_id,
                file_info_href         => $file_info_href,
                file_name_prefixes_ref => \@infile_name_prefixes,
                outdata_dir            => $active_parameter_href->{outdata_dir},
                parameter_href         => $parameter_href,
                program_name           => $program_name,
            }
        )
    );

    my $outfile_name_prefix = $io{out}{file_name_prefix};
    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my $outfile_suffix      = $io{out}{file_suffix};
    my $outfile_path        = catfile( $outfile_path_prefix . $outfile_suffix );

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $recipe_file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            directory_id                    => $sample_id,
            FILEHANDLE                      => $FILEHANDLE,
            job_id_href                     => $job_id_href,
            log                             => $log,
            process_time                    => $time,
            program_directory               => $program_name,
            program_name                    => $program_name,
            source_environment_commands_ref => \@source_environment_cmds,
        }
    );

    ### SHELL:

    ## Set up seed and fraction combination
    say {$FILEHANDLE} q{## Creating subsample filter for samtools view};

    ## Get average coverage over MT bases
    print {$FILEHANDLE} q{MT_COVERAGE=} . $BACKTICK;

    # Get depth per base
    samtools_depth(
        {
            FILEHANDLE         => $FILEHANDLE,
            infile_path        => $infile_path,
            max_depth_treshold => $MAX_DEPTH_TRESHOLD,
        }
    );

    # Pipe to AWK
    print {$FILEHANDLE} $PIPE . $SPACE;

    # Add AWK statment for calculation of avgerage coverage
    print {$FILEHANDLE} _awk_calculate_average_coverage();

    # Close statment
    say {$FILEHANDLE} $BACKTICK;

    ## Get random seed
    my $seed = int rand $MAX_LIMIT_SEED;

    ## Add seed to fraction for ~100x
    # Create bash variable
    say {$FILEHANDLE} q{SEED_FRACTION=}

      # Open statment
      . $BACKTICK

      # Lauch perl and print
      . q?perl -e "print ?

      # Add the random seed number to..
      . $seed . q{ + }

      # ...the subsample fraction, consisting of the desired subsample coverag...
      . $mt_subsample_depth

      # ...divided by the starting coverage
      . q? / $MT_COVERAGE"?

      # Close statment
      . $BACKTICK . $NEWLINE;

    ## Filter the bam file to only include a subset of reads that maps to the MT
    say {$FILEHANDLE} q{## Filter the BAM file};
    samtools_view(
        {
            exclude_reads_with_these_flags => $SAMTOOLS_UNMAPPED_READ_FLAG,
            FILEHANDLE                     => $FILEHANDLE,
            fraction                       => q{"$SEED_FRACTION"},
            infile_path                    => $infile_path,
            outfile_path                   => $outfile_path,
            with_header                    => 1,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Index new bam file
    say {$FILEHANDLE} q{## Index the subsampled BAM file};
    samtools_index(
        {
            bai_format  => 1,
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $outfile_path,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Close FILEHANDLES
    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});

    if ( $program_mode == 1 ) {

        ## Collect QC metadata info for later use
        add_program_outfile_to_sample_info(
            {
                infile           => $outfile_name_prefix,
                path             => $outfile_path,
                program_name     => q{samtools_subsample_mt},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        submit_recipe(
            {
                dependency_method       => q{sample_to_island},
                family_id               => $family_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                job_id_chain            => $job_id_chain,
                sample_id               => $sample_id,
                recipe_file_path        => $recipe_file_path,
                submission_profile      => $active_parameter_href->{submission_profile},
            }
        );
    }
    return;
}

sub _awk_calculate_average_coverage {

## Function : Writes an awk expression to an open filehandle. The awk expression calculates the average coverage based on input from samtools depth and prints it.
## Returns  : $awk_statment

    my $awk_statment =

      # Start awk
      q?awk '?

      # Sum the coverage data for each base ()
      . q?{cov += $3}?

      # Add end rule
      . q?END?

      # Divide the total coverage sum with the number of covered
      # bases (rows of output from samtools depth),
      # stored in the awk built in "NR"
      . q?{ if (NR > 0) print cov / NR }'?;

    return $awk_statment;
}

1;
