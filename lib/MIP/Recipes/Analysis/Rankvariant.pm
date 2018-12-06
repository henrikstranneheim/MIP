package MIP::Recipes::Analysis::Rankvariant;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname basename };
use File::Spec::Functions qw{ catdir catfile devnull };
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
    our $VERSION = 1.06;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ analysis_rankvariant analysis_rankvariant_unaffected analysis_rankvariant_sv analysis_rankvariant_sv_unaffected };

}

## Constants
Readonly my $AMPERSAND  => q{&};
Readonly my $ASTERISK   => q{*};
Readonly my $DASH       => q{-};
Readonly my $DOT        => q{.};
Readonly my $EMPTY_STR  => q{};
Readonly my $NEWLINE    => qq{\n};
Readonly my $PIPE       => q{|};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

sub analysis_rankvariant {

## Function : Annotate and score variants depending on mendelian inheritance, frequency and phenotype etc.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $recipe_name             => Program name
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_info_href;

    ## Default(s)
    my $case_id;
    my $temp_directory;
    my $xargs_file_counter;

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
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            store       => \$temp_directory,
            strict_type => 1,
        },
        xargs_file_counter => {
            allow       => qr{ \A\d+\z }xsm,
            default     => 0,
            store       => \$xargs_file_counter,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::Cluster qw{ check_max_core_number };
    use MIP::Cluster qw{ get_core_number };
    use MIP::File::Format::Pedigree qw{ create_fam_file };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_parameters get_recipe_attributes };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Variantcalling::Genmod
      qw{ genmod_annotate genmod_compound genmod_models genmod_score };
    use MIP::QC::Record
      qw{ add_recipe_metafile_to_sample_info add_recipe_outfile_to_sample_info };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Constant
    Readonly my $CORE_NUMBER_REQUESTED => $active_parameter_href->{max_cores_per_node};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Unpack parameters
## Get the io infiles per chain and id
    my %io = get_io_files(
        {
            id             => $case_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            stream         => q{in},
            temp_directory => $temp_directory,
        }
    );
    my $infile_name_prefix = $io{in}{file_name_prefix};
    my %infile_path        = %{ $io{in}{file_path_href} };

    my $consensus_analysis_type = $parameter_href->{cache}{consensus_analysis_type};
    my $job_id_chain            = get_recipe_attributes(
        {
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            attribute      => q{chain},
        }
    );
    my $recipe_mode = $active_parameter_href->{$recipe_name};
    my $xargs_file_path_prefix;
    my ( $core_number, $time, @source_environment_cmds ) = get_recipe_parameters(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    ## Set and get the io files per chain, id and stream
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id         => $job_id_chain,
                id               => $case_id,
                file_info_href   => $file_info_href,
                file_name_prefix => $infile_name_prefix,
                iterators_ref    => [ ( keys %infile_path ) ],
                outdata_dir      => $active_parameter_href->{outdata_dir},
                parameter_href   => $parameter_href,
                recipe_name      => $recipe_name,
                temp_directory   => $temp_directory,
            }
        )
    );

    my $outdir_path_prefix = $io{out}{dir_path_prefix};
    my %outfile_path       = %{ $io{out}{file_path_href} };
    my @outfile_paths      = @{ $io{out}{file_paths} };

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE      = IO::Handle->new();
    my $XARGSFILEHANDLE = IO::Handle->new();

    ## Get core number depending on user supplied input exists or not and max number of cores
    $core_number = get_core_number(
        {
            max_cores_per_node   => $active_parameter_href->{max_cores_per_node},
            modifier_core_number => scalar keys %infile_path,
            recipe_core_number   => $core_number,
        }
    );

    ### Detect the number of cores to use per genmod process.
    ## Limit number of cores requested to the maximum number of cores available per node
    my $genmod_core_number = check_max_core_number(
        {
            core_number_requested => $CORE_NUMBER_REQUESTED,
            max_cores_per_node    => $active_parameter_href->{max_cores_per_node},
        }
    );

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            directory_id                    => $case_id,
            log                             => $log,
            FILEHANDLE                      => $FILEHANDLE,
            job_id_href                     => $job_id_href,
            process_time                    => $time,
            recipe_directory                => $recipe_name,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => \@source_environment_cmds,
            temp_directory                  => $temp_directory,
        }
    );

    ### SHELL:

    my $case_file_path = catfile( $outdir_path_prefix, $case_id . $DOT . q{fam} );

    ## Create .fam file to be used in variant calling analyses
    create_fam_file(
        {
            active_parameter_href => $active_parameter_href,
            fam_file_path         => $case_file_path,
            FILEHANDLE            => $FILEHANDLE,
            log                   => $log,
            parameter_href        => $parameter_href,
            sample_info_href      => $sample_info_href,
        }
    );

    ## Genmod
    say {$FILEHANDLE} q{## Genmod};

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            core_number        => $genmod_core_number,
            FILEHANDLE         => $FILEHANDLE,
            file_path          => $recipe_file_path,
            recipe_info_path   => $recipe_info_path,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            xargs_file_counter => $xargs_file_counter,
        }
    );

    ## Get parameters
    ## Track which genmod modules has been processed
    my $genmod_module = $EMPTY_STR;

    my $use_vep;
    ## Use VEP annotations in compound models
    if ( $active_parameter_href->{varianteffectpredictor}
        and not $active_parameter_href->{genmod_annotate_regions} )
    {

        $use_vep = 1;
    }

    ## Process per contig
    while ( my ( $contig_index, $infile_path ) = each %infile_path ) {

        ## Get parameters
        # Restart for next contig
        $genmod_module = $EMPTY_STR;

        ## Infile
        my $genmod_indata = $infile_path;

        ## Output stream
        my $genmod_outfile_path =
          catfile( dirname( devnull() ), q{stdout} );

        my $stderrfile_path_prefix = $xargs_file_path_prefix . $DOT . $contig_index;

        ## Genmod Annotate
        $genmod_module = $UNDERSCORE . q{annotate};

        my $annotate_stderrfile_path =
          $stderrfile_path_prefix . $genmod_module . $DOT . q{stderr.txt};

        genmod_annotate(
            {
                annotate_region => $active_parameter_href->{genmod_annotate_regions},
                cadd_file_paths_ref =>
                  \@{ $active_parameter_href->{genmod_annotate_cadd_files} },
                FILEHANDLE       => $XARGSFILEHANDLE,
                infile_path      => $genmod_indata,
                outfile_path     => $genmod_outfile_path,
                spidex_file_path => $active_parameter_href->{genmod_annotate_spidex_file},
                stderrfile_path  => $annotate_stderrfile_path,
                temp_directory_path => $temp_directory,
                verbosity           => q{v},
            }
        );

        ## Pipe
        print {$XARGSFILEHANDLE} $PIPE . $SPACE;

        # Preparation for next module
        $genmod_indata = $DASH;

        ## Genmod Models
        $genmod_module .= $UNDERSCORE . q{models};

        my $models_stderrfile_path =
          $stderrfile_path_prefix . $genmod_module . $DOT . q{stderr.txt};
        genmod_models(
            {
                FILEHANDLE   => $XARGSFILEHANDLE,
                case_file    => $case_file_path,
                case_type    => $active_parameter_href->{genmod_models_case_type},
                infile_path  => $genmod_indata,
                outfile_path => catfile( dirname( devnull() ), q{stdout} ),
                reduced_penetrance_file_path =>
                  $active_parameter_href->{genmod_models_reduced_penetrance_file},
                stderrfile_path     => $models_stderrfile_path,
                temp_directory_path => $temp_directory,
                thread_number       => 4,
                vep                 => $use_vep,
                verbosity           => q{v},
                whole_gene          => $active_parameter_href->{genmod_models_whole_gene},
            }
        );
        ## Pipe
        print {$XARGSFILEHANDLE} $PIPE . $SPACE;

        ## Genmod Score
        $genmod_module .= $UNDERSCORE . q{score};

        my $score_stderrfile_path =
          $stderrfile_path_prefix . $genmod_module . $DOT . q{stderr.txt};
        genmod_score(
            {
                case_file            => $case_file_path,
                case_type            => $active_parameter_href->{genmod_models_case_type},
                FILEHANDLE           => $XARGSFILEHANDLE,
                infile_path          => $genmod_indata,
                outfile_path         => catfile( dirname( devnull() ), q{stdout} ),
                rank_result          => 1,
                rank_model_file_path => $active_parameter_href->{rank_model_file},
                stderrfile_path      => $score_stderrfile_path,
                verbosity            => q{v},
            }
        );
        ## Pipe
        print {$XARGSFILEHANDLE} $PIPE . $SPACE;

        ## Genmod Compound
        $genmod_module .= q{_compound};

        my $compound_stderrfile_path =
          $stderrfile_path_prefix . $genmod_module . $DOT . q{stderr.txt};
        genmod_compound(
            {
                FILEHANDLE          => $XARGSFILEHANDLE,
                infile_path         => $genmod_indata,
                outfile_path        => $outfile_path{$contig_index},
                stderrfile_path     => $compound_stderrfile_path,
                temp_directory_path => $temp_directory,
                verbosity           => q{v},
                vep                 => $use_vep,
            }
        );
        say {$XARGSFILEHANDLE} $NEWLINE;
    }

    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});
    close $XARGSFILEHANDLE
      or $log->logcroak(q{Could not close XARGSFILEHANDLE});

    if ( $recipe_mode == 1 ) {

        add_recipe_outfile_to_sample_info(
            {
                path             => $outfile_paths[0],
                recipe_name      => q{genmod},
                sample_info_href => $sample_info_href,
            }
        );

        # Add to Sample_info
        if ( defined $active_parameter_href->{rank_model_file} ) {

            my ($rank_model_version) =
              $active_parameter_href->{rank_model_file} =~
              m/ v(\d+[.]\d+[.]\d+ | \d+[.]\d+) /sxm;

            add_recipe_metafile_to_sample_info(
                {
                    file         => basename( $active_parameter_href->{rank_model_file} ),
                    metafile_tag => q{rank_model},
                    path         => $active_parameter_href->{rank_model_file},
                    recipe_name  => q{genmod},
                    sample_info_href => $sample_info_href,
                    version          => $rank_model_version,
                }
            );
        }
        submit_recipe(
            {
                dependency_method       => q{sample_to_case},
                case_id                 => $case_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                job_id_chain            => $job_id_chain,
                recipe_file_path        => $recipe_file_path,
                sample_ids_ref          => \@{ $active_parameter_href->{sample_ids} },
                submission_profile      => $active_parameter_href->{submission_profile},
            }
        );
    }
    return;
}

sub analysis_rankvariant_unaffected {

## Function : Annotate variants but do not score.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $recipe_name             => Program name
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_info_href;

    ## Default(s)
    my $case_id;
    my $temp_directory;
    my $xargs_file_counter;

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
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            store       => \$temp_directory,
            strict_type => 1,
        },
        xargs_file_counter => {
            allow       => qr{ \A\d+\z }xsm,
            default     => 0,
            store       => \$xargs_file_counter,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::Cluster qw{ check_max_core_number };
    use MIP::Cluster qw{ get_core_number };
    use MIP::File::Format::Pedigree qw{ create_fam_file };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_parameters get_recipe_attributes };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Variantcalling::Genmod
      qw{ genmod_annotate genmod_compound genmod_models genmod_score };
    use MIP::QC::Record
      qw{ add_recipe_metafile_to_sample_info add_recipe_outfile_to_sample_info };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Constant
    Readonly my $CORE_NUMBER_REQUESTED => $active_parameter_href->{max_cores_per_node};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Unpack parameters
## Get the io infiles per chain and id
    my %io = get_io_files(
        {
            id             => $case_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            stream         => q{in},
            temp_directory => $temp_directory,
        }
    );
    my $infile_name_prefix = $io{in}{file_name_prefix};
    my %infile_path        = %{ $io{in}{file_path_href} };

    my $consensus_analysis_type = $parameter_href->{cache}{consensus_analysis_type};
    my $job_id_chain            = get_recipe_attributes(
        {
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            attribute      => q{chain},
        }
    );
    my $recipe_mode = $active_parameter_href->{$recipe_name};
    my $xargs_file_path_prefix;
    my ( $core_number, $time, @source_environment_cmds ) = get_recipe_parameters(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    ## Set and get the io files per chain, id and stream
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id         => $job_id_chain,
                id               => $case_id,
                file_info_href   => $file_info_href,
                file_name_prefix => $infile_name_prefix,
                iterators_ref    => [ ( keys %infile_path ) ],
                outdata_dir      => $active_parameter_href->{outdata_dir},
                parameter_href   => $parameter_href,
                recipe_name      => $recipe_name,
                temp_directory   => $temp_directory,
            }
        )
    );

    my $outdir_path_prefix = $io{out}{dir_path_prefix};
    my %outfile_path       = %{ $io{out}{file_path_href} };
    my @outfile_paths      = @{ $io{out}{file_paths} };

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE      = IO::Handle->new();
    my $XARGSFILEHANDLE = IO::Handle->new();

    ## Get core number depending on user supplied input exists or not and max number of cores
    $core_number = get_core_number(
        {
            max_cores_per_node   => $active_parameter_href->{max_cores_per_node},
            modifier_core_number => scalar keys %infile_path,
            recipe_core_number   => $core_number,
        }
    );

    ### Detect the number of cores to use per genmod process.
    ## Limit number of cores requested to the maximum number of cores available per node
    my $genmod_core_number = check_max_core_number(
        {
            core_number_requested => $CORE_NUMBER_REQUESTED,
            max_cores_per_node    => $active_parameter_href->{max_cores_per_node},
        }
    );

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            directory_id                    => $case_id,
            FILEHANDLE                      => $FILEHANDLE,
            job_id_href                     => $job_id_href,
            log                             => $log,
            process_time                    => $time,
            recipe_directory                => $recipe_name,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => \@source_environment_cmds,
            temp_directory                  => $temp_directory,
        }
    );

    ### SHELL:

    my $case_file_path = catfile( $outdir_path_prefix, $case_id . $DOT . q{fam} );

    ## Create .fam file to be used in variant calling analyses
    create_fam_file(
        {
            active_parameter_href => $active_parameter_href,
            fam_file_path         => $case_file_path,
            FILEHANDLE            => $FILEHANDLE,
            log                   => $log,
            parameter_href        => $parameter_href,
            sample_info_href      => $sample_info_href,
        }
    );

    ## Genmod
    say {$FILEHANDLE} q{## Genmod};

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            core_number        => $genmod_core_number,
            FILEHANDLE         => $FILEHANDLE,
            file_path          => $recipe_file_path,
            recipe_info_path   => $recipe_info_path,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            xargs_file_counter => $xargs_file_counter,
        }
    );

    ## Track which genmod modules has been processed
    my $genmod_module = $EMPTY_STR;

    ## Process per contig
    while ( my ( $contig_index, $infile_path ) = each %infile_path ) {

        ## Get parameters
        # Restart for next contig
        $genmod_module = $EMPTY_STR;

        ## Infile
        my $genmod_indata = $infile_path;

        ## Output file
        my $genmod_outfile_path = $outfile_path{$contig_index};

        ## Genmod Annotate
        $genmod_module = $UNDERSCORE . q{annotate};

        my $stderrfile_path_prefix = $xargs_file_path_prefix . $DOT . $contig_index;
        my $annotate_stderrfile_path =
          $stderrfile_path_prefix . $genmod_module . $DOT . q{stderr.txt};

        genmod_annotate(
            {
                annotate_region => $active_parameter_href->{genmod_annotate_regions},
                cadd_file_paths_ref =>
                  \@{ $active_parameter_href->{genmod_annotate_cadd_files} },
                FILEHANDLE       => $XARGSFILEHANDLE,
                infile_path      => $genmod_indata,
                outfile_path     => $genmod_outfile_path,
                spidex_file_path => $active_parameter_href->{genmod_annotate_spidex_file},
                stderrfile_path  => $annotate_stderrfile_path,
                temp_directory_path => $temp_directory,
                verbosity           => q{v},
            }
        );
        say {$XARGSFILEHANDLE} $NEWLINE;
    }

    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});
    close $XARGSFILEHANDLE
      or $log->logcroak(q{Could not close XARGSFILEHANDLE});

    if ( $recipe_mode == 1 ) {

        add_recipe_outfile_to_sample_info(
            {
                path             => $outfile_paths[0],
                recipe_name      => q{genmod},
                sample_info_href => $sample_info_href,
            }
        );

        # Add to Sample_info
        if ( defined $active_parameter_href->{rank_model_file} ) {

            my ($rank_model_version) =
              $active_parameter_href->{rank_model_file} =~
              m/ v(\d+[.]\d+[.]\d+ | \d+[.]\d+) /sxm;

            add_recipe_metafile_to_sample_info(
                {
                    file         => basename( $active_parameter_href->{rank_model_file} ),
                    metafile_tag => q{rank_model},
                    path         => $active_parameter_href->{rank_model_file},
                    recipe_name  => q{genmod},
                    sample_info_href => $sample_info_href,
                    version          => $rank_model_version,
                }
            );
        }
        submit_recipe(
            {
                dependency_method       => q{sample_to_case},
                case_id                 => $case_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                job_id_chain            => $job_id_chain,
                recipe_file_path        => $recipe_file_path,
                sample_ids_ref          => \@{ $active_parameter_href->{sample_ids} },
                submission_profile      => $active_parameter_href->{submission_profile},
            }
        );
    }
    return;
}

sub analysis_rankvariant_sv {

## Function : Annotate and score SV variants depending on mendelian inheritance, frequency and phenotype etc.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $FILEHANDLE              => Sbatch filehandle to write to
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $recipe_name             => Program name
##          : $reference_dir_ref       => MIP reference directory {REF}
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_info_href;

    ## Default(s)
    my $case_id;
    my $reference_dir_ref;
    my $temp_directory;

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
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
        reference_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{reference_dir},
            store       => \$reference_dir_ref,
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

    use MIP::File::Format::Pedigree qw{ create_fam_file };
    use MIP::Get::Analysis qw{ get_vcf_parser_analysis_suffix };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_parameters get_recipe_attributes };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Variantcalling::Genmod
      qw{ genmod_annotate genmod_compound genmod_models genmod_score };
    use MIP::QC::Record
      qw{ add_recipe_outfile_to_sample_info add_recipe_metafile_to_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Unpack parameters
    my %io = get_io_files(
        {
            id             => $case_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            stream         => q{in},
            temp_directory => $temp_directory,
        }
    );

    my $infile_name_prefix = $io{in}{file_name_prefix};
    my @infile_paths       = @{ $io{in}{file_paths} };

    my $consensus_analysis_type = $parameter_href->{cache}{consensus_analysis_type};
    my $job_id_chain            = get_recipe_attributes(
        {
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            attribute      => q{chain},
        }
    );
    my $recipe_mode = $active_parameter_href->{$recipe_name};
    my ( $core_number, $time, @source_environment_cmds ) = get_recipe_parameters(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    my @vcfparser_analysis_types = get_vcf_parser_analysis_suffix(
        {
            vcfparser_outfile_count =>
              $active_parameter_href->{sv_vcfparser_outfile_count},
        }
    );

    ## Set and get the io files per chain, id and stream
    my @set_outfile_name_prefixes =
      map { $infile_name_prefix . $_ } @vcfparser_analysis_types;
    ## Set and get the io files per chain, id and stream
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id         => $job_id_chain,
                id               => $case_id,
                file_info_href   => $file_info_href,
                file_name_prefix => $infile_name_prefix,
                iterators_ref    => \@vcfparser_analysis_types,
                outdata_dir      => $active_parameter_href->{outdata_dir},
                parameter_href   => $parameter_href,
                recipe_name      => $recipe_name,
                temp_directory   => $temp_directory,
            }
        )
    );

    my $outdir_path_prefix = $io{out}{dir_path_prefix};
    my @outfile_suffixes   = @{ $io{out}{file_suffixes} };
    my @outfile_paths      = @{ $io{out}{file_paths} };

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Limit number of cores requested to the maximum number of cores available per node
    my $genmod_core_number = check_max_core_number(
        {
            core_number_requested => $core_number,
            max_cores_per_node    => $active_parameter_href->{max_cores_per_node},
        }
    );

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            directory_id                    => $case_id,
            FILEHANDLE                      => $FILEHANDLE,
            job_id_href                     => $job_id_href,
            log                             => $log,
            process_time                    => $time,
            recipe_directory                => $recipe_name,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => \@source_environment_cmds,
            temp_directory                  => $temp_directory,
        }
    );

    ### SHELL:

    my $fam_file_path = catfile( $outdir_path_prefix, $case_id . $DOT . q{fam} );

    ## Create .fam file
    create_fam_file(
        {
            active_parameter_href => $active_parameter_href,
            fam_file_path         => $fam_file_path,
            FILEHANDLE            => $FILEHANDLE,
            log                   => $log,
            parameter_href        => $parameter_href,
            sample_info_href      => $sample_info_href,
        }
    );

    say {$FILEHANDLE} q{## Genmod};

    ## Get parameters
    my $use_vep;
    ## Use VEP annotations in compound models
    if ( $active_parameter_href->{sv_varianteffectpredictor}
        and not $active_parameter_href->{sv_genmod_annotate_regions} )
    {

        $use_vep = 1;
    }

  INFILE:
    while ( my ( $infile_index, $infile_path ) = each @infile_paths ) {

        ## Alias
        # For pipe
        my $genmod_indata = $infile_path;

        # Outfile
        my $genmod_outfile_path =
          catfile( dirname( devnull() ), q{stdout} );

        ## Genmod annotate
        my $genmod_module = $UNDERSCORE . q{annotate};
        genmod_annotate(
            {
                annotate_region => $active_parameter_href->{sv_genmod_annotate_regions},
                FILEHANDLE      => $FILEHANDLE,
                infile_path     => $genmod_indata,
                outfile_path    => $genmod_outfile_path,
                stderrfile_path => $recipe_info_path
                  . $genmod_module
                  . $DOT
                  . q{stderr.txt},
                temp_directory_path => $temp_directory,
                verbosity           => q{v},
            }
        );

        # Pipe
        print {$FILEHANDLE} $PIPE . $SPACE;

        ## Get parameters
        # Preparation for next module
        $genmod_indata = $DASH;

        ## Genmod models
        $genmod_module .= $UNDERSCORE . q{models};

        genmod_models(
            {
                FILEHANDLE   => $FILEHANDLE,
                case_file    => $fam_file_path,
                case_type    => $active_parameter_href->{sv_genmod_models_case_type},
                infile_path  => $genmod_indata,
                outfile_path => catfile( dirname( devnull() ), q{stdout} ),
                reduced_penetrance_file_path =>
                  $active_parameter_href->{sv_genmod_models_reduced_penetrance_file},
                stderrfile_path => $recipe_info_path
                  . $genmod_module
                  . $outfile_suffixes[$infile_index]
                  . $DOT
                  . q{stderr.txt},
                temp_directory_path => $temp_directory,
                thread_number       => 4,
                vep                 => $use_vep,
                verbosity           => q{v},
                whole_gene => $active_parameter_href->{sv_genmod_models_whole_gene},
            }
        );

        # Pipe
        print {$FILEHANDLE} $PIPE . $SPACE;

        ## Genmod score
        $genmod_module .= $UNDERSCORE . q{score};
        genmod_score(
            {
                FILEHANDLE   => $FILEHANDLE,
                case_file    => $fam_file_path,
                case_type    => $active_parameter_href->{sv_genmod_models_case_type},
                infile_path  => $genmod_indata,
                outfile_path => catfile( dirname( devnull() ), q{stdout} ),
                rank_model_file_path => $active_parameter_href->{sv_rank_model_file},
                rank_result          => 1,
                stderrfile_path      => $recipe_info_path
                  . $genmod_module
                  . $outfile_suffixes[$infile_index]
                  . $DOT
                  . q{stderr.txt},
                verbosity => q{v},
            }
        );

        # Pipe
        print {$FILEHANDLE} $PIPE . $SPACE;

        ## Genmod compound
        $genmod_module .= $UNDERSCORE . q{compound};
        genmod_compound(
            {
                FILEHANDLE      => $FILEHANDLE,
                infile_path     => $genmod_indata,
                outfile_path    => $outfile_paths[$infile_index],
                stderrfile_path => $recipe_info_path
                  . $genmod_module
                  . $outfile_suffixes[$infile_index]
                  . $DOT
                  . q{stderr.txt},
                temp_directory_path => $temp_directory,
                vep                 => $use_vep,
                verbosity           => q{v},
            }
        );

        say {$FILEHANDLE} $AMPERSAND . $NEWLINE;

        if ( $recipe_mode == 1 ) {

            add_recipe_outfile_to_sample_info(
                {
                    path             => $outfile_paths[$infile_index],
                    recipe_name      => q{sv_genmod},
                    sample_info_href => $sample_info_href,
                }
            );
        }
    }
    say {$FILEHANDLE} q{wait} . $NEWLINE;

    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});

    if ( $recipe_mode == 1 ) {

        ## Add to sample_info
        if ( defined $active_parameter_href->{sv_rank_model_file} ) {

            my ($sv_rank_model_version) =
              $active_parameter_href->{sv_rank_model_file} =~
              m/ v(\d+[.]\d+[.]\d+ | \d+[.]\d+) /sxm;

            add_recipe_metafile_to_sample_info(
                {
                    file => basename( $active_parameter_href->{sv_rank_model_file} ),
                    metafile_tag     => q{sv_rank_model},
                    path             => $active_parameter_href->{sv_rank_model_file},
                    recipe_name      => q{sv_genmod},
                    sample_info_href => $sample_info_href,
                    version          => $sv_rank_model_version,
                }
            );
        }

        submit_recipe(
            {
                dependency_method       => q{sample_to_case},
                case_id                 => $case_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                job_id_chain            => $job_id_chain,
                recipe_file_path        => $recipe_file_path,
                sample_ids_ref          => \@{ $active_parameter_href->{sample_ids} },
                submission_profile      => $active_parameter_href->{submission_profile},
            }
        );
    }
    return;
}

sub analysis_rankvariant_sv_unaffected {

## Function : Annotate variants using genmod anotate only.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id               => Family id
##          : $FILEHANDLE              => Sbatch filehandle to write to
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $recipe_name            => Program name
##          : $reference_dir_ref       => MIP reference directory {REF}
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_info_href;

    ## Default(s)
    my $case_id;
    my $reference_dir_ref;
    my $temp_directory;

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
            strict_type => 1,
            store       => \$parameter_href,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
        reference_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{reference_dir},
            store       => \$reference_dir_ref,
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

    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_parameters get_recipe_attributes };
    use MIP::IO::Files qw{ migrate_files };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Variantcalling::Genmod qw{ genmod_annotate };
    use MIP::QC::Record
      qw{ add_recipe_outfile_to_sample_info add_recipe_metafile_to_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Unpack parameters
    my %io = get_io_files(
        {
            id             => $case_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            stream         => q{in},
            temp_directory => $temp_directory,
        }
    );

    my $indir_path_prefix  = $io{in}{dir_path_prefix};
    my $infile_name_prefix = $io{in}{file_name_prefix};
    my @infile_paths       = @{ $io{in}{file_paths} };

    my $consensus_analysis_type = $parameter_href->{cache}{consensus_analysis_type};
    my $job_id_chain            = get_recipe_attributes(
        {
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            attribute      => q{chain},
        }
    );
    my $recipe_mode = $active_parameter_href->{$recipe_name};
    my ( $core_number, $time, @source_environment_cmds ) = get_recipe_parameters(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    my @vcfparser_analysis_types = get_vcf_parser_analysis_suffix(
        {
            vcfparser_outfile_count =>
              $active_parameter_href->{sv_vcfparser_outfile_count},
        }
    );

    ## Set and get the io files per chain, id and stream
    my @set_outfile_name_prefixes =
      map { $infile_name_prefix . $_ } @vcfparser_analysis_types;
    ## Set and get the io files per chain, id and stream
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id         => $job_id_chain,
                id               => $case_id,
                file_info_href   => $file_info_href,
                file_name_prefix => $infile_name_prefix,
                iterators_ref    => \@vcfparser_analysis_types,
                outdata_dir      => $active_parameter_href->{outdata_dir},
                parameter_href   => $parameter_href,
                recipe_name      => $recipe_name,
                temp_directory   => $temp_directory,
            }
        )
    );

    my @outfile_suffixes = @{ $io{out}{file_suffixes} };
    my @outfile_paths    = @{ $io{out}{file_paths} };

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Limit number of cores requested to the maximum number of cores available per node
    my $genmod_core_number = check_max_core_number(
        {
            core_number_requested => $core_number,
            max_cores_per_node    => $active_parameter_href->{max_cores_per_node},
        }
    );

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            directory_id                    => $case_id,
            FILEHANDLE                      => $FILEHANDLE,
            job_id_href                     => $job_id_href,
            log                             => $log,
            process_time                    => $time,
            recipe_directory                => $recipe_name,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => \@source_environment_cmds,
            temp_directory                  => $temp_directory,
        }
    );

    ### SHELL:

    say {$FILEHANDLE} q{## Genmod};
    ## Get parameters

  INFILE:
    while ( my ( $infile_index, $infile_path ) = each @infile_paths ) {

        ## Annotate
        my $genmod_module = q{_annotate};
        genmod_annotate(
            {
                annotate_region => $active_parameter_href->{sv_genmod_annotate_regions},
                FILEHANDLE      => $FILEHANDLE,
                infile_path     => $infile_path,
                outfile_path    => $outfile_paths[$infile_index],
                stderrfile_path => $recipe_info_path
                  . $genmod_module
                  . $outfile_suffixes[$infile_index]
                  . $DOT
                  . q{stderr.txt},
                temp_directory_path => $temp_directory,
                verbosity           => q{v},
            }
        );

        say {$FILEHANDLE} $AMPERSAND . $NEWLINE;

        if ( $recipe_mode == 1 ) {

            add_recipe_outfile_to_sample_info(
                {
                    path             => $outfile_paths[$infile_index],
                    recipe_name      => q{sv_genmod},
                    sample_info_href => $sample_info_href,
                }
            );
        }
    }
    say {$FILEHANDLE} q{wait} . $NEWLINE;

    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});

    if ( $recipe_mode == 1 ) {

        ## Add to Sample_info
        if ( defined $active_parameter_href->{sv_rank_model_file} ) {

            my ($sv_rank_model_version) =
              $active_parameter_href->{sv_rank_model_file} =~
              m/ v(\d+[.]\d+[.]\d+ | \d+[.]\d+) /sxm;

            add_recipe_metafile_to_sample_info(
                {
                    file => basename( $active_parameter_href->{sv_rank_model_file} ),
                    metafile_tag     => q{sv_rank_model},
                    path             => $active_parameter_href->{sv_rank_model_file},
                    recipe_name      => q{sv_genmod},
                    sample_info_href => $sample_info_href,
                    version          => $sv_rank_model_version,
                }
            );

        }

        submit_recipe(
            {
                dependency_method       => q{sample_to_case},
                case_id                 => $case_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                job_id_chain            => $job_id_chain,
                recipe_file_path        => $recipe_file_path,
                sample_ids_ref          => \@{ $active_parameter_href->{sample_ids} },
                submission_profile      => $active_parameter_href->{submission_profile},
            }
        );
    }
    return;
}

1;
