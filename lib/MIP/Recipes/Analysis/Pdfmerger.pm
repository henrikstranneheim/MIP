package MIP::Recipes::Analysis::Pdfmerger;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile devnull };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $DOT $LOG_NAME $NEWLINE $PIPE $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_merge_fusion_reports };

}

sub analysis_merge_fusion_reports {

## Function : Generate a single fusion report from arriba fusion calls
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $case_id               => Family id
##          : $file_info_href        => File_info hash {REF}
##          : $job_id_href           => Job id hash {REF}
##          : $parameter_href        => Parameter hash {REF}
##          : $profile_base_command  => Submission profile base command
##          : $recipe_name           => Recipe name
##          : $sample_info_href      => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_info_href;

    ## Default(s)
    my $case_id;
    my $profile_base_command;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        case_id => {
            default     => $arg_href->{active_parameter_href}{case_id},
            store       => \$case_id,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
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
        profile_base_command => {
            default     => q{sbatch},
            store       => \$profile_base_command,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
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

    use MIP::Analysis qw{ get_vcf_parser_analysis_suffix };
    use MIP::File_info qw{ get_io_files parse_io_outfiles };
    use MIP::Program::Pdfmerger qw{ pdfmerger };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Sample_info qw{
      get_pedigree_sample_id_attributes
      set_file_path_to_store
      set_recipe_metafile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my @report_types = get_vcf_parser_analysis_suffix(
        {
            vcfparser_outfile_count => $active_parameter_href->{fusion_outfile_count},
        }
    );

    my %infile;
  SAMPLE_ID:
    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        ## Get the io infiles per chain and id
        my %sample_io = get_io_files(
            {
                id             => $sample_id,
                file_info_href => $file_info_href,
                parameter_href => $parameter_href,
                recipe_name    => q{fusion_report},
                stream         => q{out},
            }
        );
        push @{ $infile{research} }, $sample_io{out}{file_paths}->[0];

        if ( @report_types > 1 ) {

            push @{ $infile{clinical} }, $sample_io{out}{file_paths}->[1];
        }
    }

    my %recipe = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    ## Set and get the io files per chain, id and stream
    my %io = (
        parse_io_outfiles(
            {
                chain_id         => $recipe{job_id_chain},
                id               => $case_id,
                file_info_href   => $file_info_href,
                file_name_prefix => $case_id,
                iterators_ref    => \@report_types,
                outdata_dir      => $active_parameter_href->{outdata_dir},
                parameter_href   => $parameter_href,
                recipe_name      => $recipe_name,
            }
        )
    );
    my $outdir_path = $io{out}{dir_path};
    my %outfile     = (
        research => {
            file_path => $io{out}{file_paths}->[0],
            file_name => $io{out}{file_names}->[0],
        },
    );
    if ( @report_types > 1 ) {

        $outfile{clinical}{file_path} = $io{out}{file_paths}->[1];
        $outfile{clinical}{file_name} = $io{out}{file_names}->[1];
    }

    ### Filehandles
    ## Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $recipe{core_number},
            directory_id          => $case_id,
            filehandle            => $filehandle,
            job_id_href           => $job_id_href,
            memory_allocation     => $recipe{memory},
            process_time          => $recipe{time},
            recipe_directory      => $recipe_name,
            recipe_name           => $recipe_name,
        }
    );

    #### SHELL:

    say {$filehandle} q{## } . $recipe_name;

  TAG:
    foreach my $tag ( keys %infile ) {

        pdfmerger(
            {
                filehandle       => $filehandle,
                infile_paths_ref => $infile{$tag},
                orientation      => q{landscape},
                outdir_path      => $outdir_path,
                outfile_name     => $outfile{$tag}{file_name},
                write_filenames  => 1,
            }
        );
        say {$filehandle} $NEWLINE;
    }

    ## Close filehandle
    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe{mode} == 1 ) {

      TAG:
        foreach my $tag ( keys %outfile ) {

            ## Collect QC metadata info for later use
            set_recipe_metafile_in_sample_info(
                {
                    sample_info_href => $sample_info_href,
                    recipe_name      => $recipe_name,
                    metafile_tag     => $tag,
                    path             => $outfile{$tag}{file_path},
                }
            );

            set_file_path_to_store(
                {
                    format           => q{meta},
                    id               => $case_id,
                    path             => $outfile{$tag}{file_path},
                    recipe_name      => $recipe_name,
                    sample_info_href => $sample_info_href,
                    tag              => $tag,
                }
            );

        }

        submit_recipe(
            {
                base_command                      => $profile_base_command,
                case_id                           => $case_id,
                dependency_method                 => q{case_to_island},
                job_id_chain                      => $recipe{job_id_chain},
                job_id_href                       => $job_id_href,
                job_reservation_name              => $active_parameter_href->{job_reservation_name},
                log                               => $log,
                max_parallel_processes_count_href =>
                  $file_info_href->{max_parallel_processes_count},
                recipe_file_path   => $recipe_file_path,
                sample_ids_ref     => \@{ $active_parameter_href->{sample_ids} },
                submission_profile => $active_parameter_href->{submission_profile},
            }
        );
    }
    return 1;
}

1;
