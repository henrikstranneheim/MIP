package MIP::Recipes::Analysis::Salmon_quant;

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
use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_salmon_quant };

}

## Constants
Readonly my $ASTERIX    => q{*};
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $UNDERSCORE => q{_};

sub analysis_salmon_quant {

## Function : Transcript quantification using salmon quant
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
##          : $temp_directory          => Temporary directory

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
    my $temp_directory;

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
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::File qw{ get_file_suffix };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::IO::Files qw{ migrate_file };
    use MIP::Program::Variantcalling::Salmon qw{ salmon_quant };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_step_in_parallel };
    use MIP::QC::Record
      qw{ add_program_outfile_to_sample_info add_processing_metafile_to_sample_info };
    use MIP::Set::File qw{ set_file_suffix };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program mode
    my $program_mode = $active_parameter_href->{$program_name};

    ## Get job_id_chain
    my $job_id_chain = $parameter_href->{$program_name}{chain};

    my ( $core_number, $time, @source_environment_cmds ) =
      get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name          => $program_name,
        }
      );

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Get infiles
    my @infiles = @{ $file_info_href->{$sample_id}{mip_infiles} };

    ## Directories
    my $insample_directory  = $file_info_href->{$sample_id}{mip_infiles_dir};
    my $outsample_directory = catdir( $active_parameter_href->{outdata_dir},
        $sample_id, $active_parameter_href->{outaligner_dir} );

    ## Assign file_tags
    my $outfile_tag =
      $file_info_href->{$sample_id}{$program_name}{file_tag};

    ### Assign suffix
    ## Set file suffix for next module within jobid chain
    my $outfile_suffix = set_file_suffix(
        {
            file_suffix    => $parameter_href->{$program_name}{outfile_suffix},
            job_id_chain   => $job_id_chain,
            parameter_href => $parameter_href,
            suffix_key     => q{sailfish_quantification_file_suffix},
        }
    );

    # Too avoid adjusting infile_index in submitting to jobs
    my $paired_end_tracker = 0;

    ## Perform per single-end or read pair
  INFILE_PREFIX:
    while ( my ( $infile_index, $infile_prefix ) =
        each @{ $infile_lane_prefix_href->{$sample_id} } )
    {

        ## Assign file tags
        my $file_path_prefix = catfile( $temp_directory, $infile_prefix );
        my $outfile_path_prefix = $file_path_prefix . $outfile_tag;

        # Collect paired-end or single-end sequence run mode
        my $sequence_run_mode =
          $sample_info_href->{sample}{$sample_id}{file}{$infile_prefix}
          {sequence_run_type};

        # Collect interleaved info
        my $interleaved_fastq_file =
          $sample_info_href->{sample}{$sample_id}{file}{$infile_prefix}
          {interleaved};

        ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
        my ( $file_name, $program_info_path ) = setup_script(
            {
                active_parameter_href => $active_parameter_href,
                core_number           => $core_number,
                directory_id          => $sample_id,
                FILEHANDLE            => $FILEHANDLE,
                job_id_href           => $job_id_href,
                log                   => $log,
                program_directory =>
                  lc $active_parameter_href->{outaligner_dir},
                program_name                    => $program_name,
                process_time                    => $time,
                source_environment_commands_ref => \@source_environment_cmds,
                temp_directory                  => $temp_directory,
            }
        );

        ## Copies file to temporary directory.
        say {$FILEHANDLE} q{## Copy file(s) to temporary directory};

        # Read 1
        my $insample_dir_fastqc_path_read_one =
          catfile( $insample_directory, $infiles[$paired_end_tracker] );
        migrate_file(
            {
                FILEHANDLE   => $FILEHANDLE,
                infile_path  => $insample_dir_fastqc_path_read_one,
                outfile_path => $temp_directory,
            }
        );

        # If second read direction is present
        if ( $sequence_run_mode eq q{paired-end} ) {

            my $insample_dir_fastqc_path_read_two =
              catfile( $insample_directory,
                $infiles[ $paired_end_tracker + 1 ] );

            # Read 2
            migrate_file(
                {
                    FILEHANDLE   => $FILEHANDLE,
                    infile_path  => $insample_dir_fastqc_path_read_two,
                    outfile_path => $temp_directory,
                }
            );
        }
        say {$FILEHANDLE} q{wait}, $NEWLINE;

        ## Salmon quant
        say {$FILEHANDLE} q{## Quantifying transcripts using } . $program_name;

        ### Get parameters

        ## Infile(s)
        my @fastq_files =
          ( catfile( $temp_directory, $infiles[$paired_end_tracker] ) );

        # If second read direction is present
        if ( $sequence_run_mode eq q{paired-end} ) {

            # Increment to collect correct read 2 from %infile
            $paired_end_tracker = $paired_end_tracker + 1;
            push @fastq_files,
              catfile( $temp_directory, $infiles[$paired_end_tracker] );

            #catfile( $temp_directory, $infiles[$paired_end_tracker] );
        }
        my $referencefile_dir_path = $active_parameter_href->{reference_dir};

        # If second read direction is present
        if ( $sequence_run_mode ne q{paired-end} ) {
            salmon_quant(
                {
                    FILEHANDLE        => $FILEHANDLE,
                    index_path        => $referencefile_dir_path,
                    outfile_path      => $outfile_path_prefix,
                    read_1_fastq_path => $fastq_files[0],
                },
            );
            say {$FILEHANDLE} $NEWLINE;
        }
        else {
            salmon_quant(
                {
                    FILEHANDLE        => $FILEHANDLE,
                    index_path        => $referencefile_dir_path,
                    outfile_path      => $outfile_path_prefix,
                    read_1_fastq_path => $fastq_files[0],
                    read_2_fastq_path => $fastq_files[1],
                },
            );
            say {$FILEHANDLE} $NEWLINE;
        }

        ## Increment paired end tracker
        $paired_end_tracker++;

        ## Copies file from temporary directory.
        say {$FILEHANDLE} q{## Copy file from temporary directory};
        migrate_file(
            {
                FILEHANDLE   => $FILEHANDLE,
                infile_path  => $outfile_path_prefix . $ASTERIX,
                outfile_path => $outsample_directory,
            }
        );
        say {$FILEHANDLE} q{wait}, $NEWLINE;

        ## Close FILEHANDLES
        close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});

        if ( $program_mode == 1 ) {

            my $program_outfile_path =
              $outfile_path_prefix . $UNDERSCORE . q{sf};

            ## Collect QC metadata info for later use
            add_program_outfile_to_sample_info(
                {
                    path             => $program_outfile_path,
                    program_name     => $program_name,
                    sample_id        => $sample_id,
                    sample_info_href => $sample_info_href,
                }
            );

            slurm_submit_job_sample_id_dependency_step_in_parallel(
                {
                    family_id               => $family_id,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    log                     => $log,
                    path                    => $job_id_chain,
                    sample_id               => $sample_id,
                    sbatch_file_name        => $file_name,
                    sbatch_script_tracker   => $infile_index
                }
            );
        }
    }
    return;
}

1;
