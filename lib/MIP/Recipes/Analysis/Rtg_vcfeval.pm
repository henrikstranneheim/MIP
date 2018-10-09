package MIP::Recipes::Analysis::Rtg_vcfeval;

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
    our @EXPORT_OK = qw{ analysis_rtg_vcfeval };

}

## Constants
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $UNDERSCORE => q{_};

sub analysis_rtg_vcfeval {

## Function : Evaluation of vcf variants using rtg
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
            strict_type => 1,
            store       => \$temp_directory,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_module_parameters get_program_attributes };
    use MIP::Gnu::Coreutils qw{ gnu_rm  };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Program::Qc::Rtg qw{ rtg_vcfeval };
    use MIP::Program::Variantcalling::Bcftools
      qw{ bcftools_rename_vcf_samples bcftools_view_and_index_vcf };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_family_dead_end };
    use MIP::QC::Record qw{ add_program_outfile_to_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Return if not a nist_id sample
    return if ( not $sample_id =~ /$active_parameter_href->{nist_id}/sxm );

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Unpack parameters
## Get the io infiles per chain and id
    my %io = get_io_files(
        {
            id             => $family_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            program_name   => $program_name,
            stream         => q{in},
            temp_directory => $temp_directory,
        }
    );
    my $infile_name_prefix = $io{in}{file_name_prefix};
    my $infile_path        = $io{in}{file_path};

    my $job_id_chain = get_program_attributes(
        {
            parameter_href => $parameter_href,
            program_name   => $program_name,
            attribute      => q{chain},
        }
    );
    my $nist_id      = $active_parameter_href->{nist_id};
    my $program_mode = $active_parameter_href->{$program_name};
    my ( $core_number, $time, @source_environment_cmds ) =
      get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name          => $program_name,
        }
      );

    ## Set and get the io files per chain, id and stream
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id         => $job_id_chain,
                id               => $family_id,
                file_info_href   => $file_info_href,
                file_name_prefix => $infile_name_prefix,
                iterators_ref    => [$nist_id],
                outdata_dir      => $active_parameter_href->{outdata_dir},
                parameter_href   => $parameter_href,
                program_name     => $program_name,
                temp_directory   => $temp_directory,
            }
        )
    );

    my $outdir_path_prefix  = $io{out}{dir_path_prefix};
    my $outfile_path_prefix = $io{out}{file_path_prefix};

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    my $nist_file_path = catfile( $temp_directory, q{nist} );

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            directory_id                    => $family_id,
            FILEHANDLE                      => $FILEHANDLE,
            job_id_href                     => $job_id_href,
            log                             => $log,
            process_time                    => $time,
            program_directory               => $program_name,
            program_name                    => $program_name,
            source_environment_commands_ref => \@source_environment_cmds,
        }
    );

    say {$FILEHANDLE} q{## Adding sample name to baseline calls};
    bcftools_rename_vcf_samples(
        {
            FILEHANDLE => $FILEHANDLE,
            index      => 1,
            index_type => q{tbi},
            infile => $active_parameter_href->{nist_high_confidence_call_set},
            outfile_path_prefix => $nist_file_path . $UNDERSCORE . q{refrm},
            output_type         => q{z},
            temp_directory      => $temp_directory,
            sample_ids_ref      => [$nist_id],
        }
    );

    say {$FILEHANDLE} q{## Compressing and indexing sample calls};
    bcftools_view_and_index_vcf(
        {
            FILEHANDLE          => $FILEHANDLE,
            index               => 1,
            index_type          => q{tbi},
            infile_path         => $infile_path,
            outfile_path_prefix => $outfile_path_prefix,
            output_type         => q{z},
        }
    );

    say {$FILEHANDLE} q{## Remove potential old Rtg vcfeval outdir};
    my $rtg_outdirectory_path = catfile( $outdir_path_prefix, $sample_id );
    gnu_rm(
        {
            FILEHANDLE  => $FILEHANDLE,
            force       => 1,
            infile_path => $rtg_outdirectory_path,
            recursive   => 1,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    say {$FILEHANDLE} q{## Rtg vcfeval};
    rtg_vcfeval(
        {
            baselinefile_path => $nist_file_path
              . $UNDERSCORE
              . q{refrm.vcf.gz},
            callfile_path => $outfile_path_prefix . $DOT . q{vcf.gz},
            eval_region_file_path =>
              $active_parameter_href->{nist_high_confidence_call_set_bed},
            FILEHANDLE           => $FILEHANDLE,
            outputdirectory_path => $rtg_outdirectory_path,
            sample_id            => $nist_id,
            sdf_template_file_path =>
              $active_parameter_href->{rtg_vcfeval_reference_genome}
              . $file_info_href->{rtg_vcfeval_reference_genome}[0]
            ,    # Only one directory for sdf
        }
    );

    ## Close FILEHANDLES
    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});

    if ( $program_mode == 1 ) {

## Collect QC metadata info for later use
        add_program_outfile_to_sample_info(
            {
                outdirectory     => $rtg_outdirectory_path,
                program_name     => $program_name,
                sample_info_href => $sample_info_href,
            }
        );

        slurm_submit_job_sample_id_dependency_family_dead_end(
            {
                family_id               => $family_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                path                    => $job_id_chain,
                sample_ids_ref   => \@{ $active_parameter_href->{sample_ids} },
                sbatch_file_name => $file_path,
            }
        );
    }
    return;
}

1;
