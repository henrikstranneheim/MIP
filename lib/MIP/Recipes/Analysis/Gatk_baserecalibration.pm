package MIP::Recipes::Analysis::Gatk_baserecalibration;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use File::Spec::Functions qw{ catdir catfile };
use Params::Check qw{ allow check last_error };
use POSIX;
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ %ANALYSIS $ASTERISK $DOT $LOG_NAME $NEWLINE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.25;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ analysis_gatk_baserecalibration analysis_gatk_baserecalibration_panel analysis_gatk_baserecalibration_rna };

}

## Constants
Readonly my $JAVA_MEMORY_ALLOCATION => 6;
Readonly my $JAVA_GUEST_OS_MEMORY   => $ANALYSIS{JAVA_GUEST_OS_MEMORY};
Readonly my $MINUS_ONE              => -1;

sub analysis_gatk_baserecalibration {

## Function : GATK baserecalibrator/GatherBQSRReports/ApplyBQSR to recalibrate bases before variant calling. BaseRecalibrator/GatherBQSRReports/ApplyBQSR will be executed within the same sbatch script.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File info hash {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Program name
##          : $sample_id               => Sample id
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_info_href;
    my $sample_id;

    ## Default(s)
    my $case_id;
    my $profile_base_command;
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
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            store       => \$temp_directory,
            strict_type => 1,
        },
        xargs_file_counter => {
            allow       => qr/ \A \d+ \z /xsm,
            default     => 0,
            store       => \$xargs_file_counter,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Cluster qw{ get_parallel_processes };
    use MIP::Get::File qw{ get_merged_infile_prefix get_io_files };
    use MIP::Get::Parameter
      qw{ get_gatk_intervals get_recipe_attributes get_recipe_resources };
    use MIP::Program::Gnu::Coreutils qw{ gnu_cp };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Gatk
      qw{ gatk_applybqsr gatk_baserecalibrator gatk_gatherbqsrreports };
    use MIP::Program::Picardtools qw{ picardtools_gatherbamfiles };
    use MIP::Program::Samtools qw{ samtools_index samtools_view };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Sample_info qw{
      set_file_path_to_store
      set_recipe_outfile_in_sample_info
      set_recipe_metafile_in_sample_info
    };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Unpack parameters
    ## Get the io infiles per chain and id
    my %io = get_io_files(
        {
            id             => $sample_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            stream         => q{in},
        }
    );
    my $indir_path_prefix  = $io{in}{dir_path_prefix};
    my $infile_suffix      = $io{in}{file_suffix};
    my $infile_name_prefix = $io{in}{file_name_prefix};
    my %infile_path        = %{ $io{in}{file_path_href} };

    my %rec_atr = get_recipe_attributes(
        {
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
        }
    );
    my $analysis_type      = $active_parameter_href->{analysis_type}{$sample_id};
    my $job_id_chain       = $rec_atr{chain};
    my $recipe_mode        = $active_parameter_href->{$recipe_name};
    my $referencefile_path = $active_parameter_href->{human_genome_reference};
    my $xargs_file_path_prefix;
    my %recipe_resource = get_recipe_resources(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );
    my $core_number = $recipe_resource{core_number};

    ## Add merged infile name prefix after merging all BAM files per sample_id
    my $merged_infile_prefix = get_merged_infile_prefix(
        {
            file_info_href => $file_info_href,
            sample_id      => $sample_id,
        }
    );

    ## Outpaths
    ## Assign suffix
    my $outfile_suffix = $rec_atr{outfile_suffix};
    my $outsample_directory =
      catdir( $active_parameter_href->{outdata_dir}, $sample_id, $recipe_name );
    my $outfile_tag =
      $file_info_href->{$sample_id}{$recipe_name}{file_tag};
    my @outfile_paths =
      map {
        catdir( $outsample_directory,
            $merged_infile_prefix . $outfile_tag . $DOT . $_ . $outfile_suffix )
      } @{ $file_info_href->{bam_contigs_size_ordered} };

    ## Set and get the io files per chain, id and stream
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id       => $job_id_chain,
                id             => $sample_id,
                file_info_href => $file_info_href,
                file_paths_ref => \@outfile_paths,
                parameter_href => $parameter_href,
                recipe_name    => $recipe_name,
                temp_directory => $temp_directory,
            }
        )
    );

    my $outdir_path_prefix  = $io{out}{dir_path_prefix};
    my $outfile_name_prefix = $io{out}{file_name_prefix};
    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my %outfile_path        = %{ $io{out}{file_path_href} };

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle      = IO::Handle->new();
    my $xargsfilehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            directory_id                    => $sample_id,
            filehandle                      => $filehandle,
            job_id_href                     => $job_id_href,
            memory_allocation               => $recipe_resource{memory},
            process_time                    => $recipe_resource{time},
            recipe_directory                => $recipe_name,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => $recipe_resource{load_env_ref},
            temp_directory                  => $temp_directory,
        }
    );

    ### SHELL:

    ## Generate gatk intervals. Chromosomes for WGS/WTS and paths to contig_bed_files for WES
    my %gatk_intervals = get_gatk_intervals(
        {
            analysis_type         => $analysis_type,
            contigs_ref           => \@{ $file_info_href->{bam_contigs_size_ordered} },
            exome_target_bed_href => $active_parameter_href->{exome_target_bed},
            filehandle            => $filehandle,
            file_ending           => $file_info_href->{exome_target_bed}[0],
            max_cores_per_node    => $core_number,
            log                   => $log,
            outdirectory          => $outdir_path_prefix,
            reference_dir         => $active_parameter_href->{reference_dir},
            sample_id             => $sample_id,
        }
    );

    my $process_memory_allocation = $JAVA_MEMORY_ALLOCATION + $JAVA_GUEST_OS_MEMORY;

    # Constrain parallelization to match available memory
    my $parallel_processes = get_parallel_processes(
        {
            core_number               => $core_number,
            process_memory_allocation => $process_memory_allocation,
            recipe_memory_allocation  => $recipe_resource{memory},
        }
    );

    ## GATK BaseRecalibrator
    say {$filehandle} q{## GATK BaseRecalibrator};

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            core_number        => $parallel_processes,
            filehandle         => $filehandle,
            file_path          => $recipe_file_path,
            recipe_info_path   => $recipe_info_path,
            xargsfilehandle    => $xargsfilehandle,
            xargs_file_counter => $xargs_file_counter,
        }
    );

    my @base_quality_score_recalibration_files;
  CONTIG:
    foreach my $contig ( @{ $file_info_href->{bam_contigs_size_ordered} } ) {

        my $base_quality_score_recalibration_file =
          $outfile_path_prefix . $DOT . $contig . $DOT . q{grp};

        ## Add for gathering base recal files later
        push @base_quality_score_recalibration_files,
          $base_quality_score_recalibration_file;
        my $stderrfile_path =
          $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};
        gatk_baserecalibrator(
            {
                filehandle           => $xargsfilehandle,
                infile_path          => $infile_path{$contig},
                intervals_ref        => $gatk_intervals{$contig},
                java_use_large_pages => $active_parameter_href->{java_use_large_pages},
                known_sites_ref =>
                  \@{ $active_parameter_href->{gatk_baserecalibration_known_sites} },
                memory_allocation  => q{Xmx} . $JAVA_MEMORY_ALLOCATION . q{g},
                outfile_path       => $base_quality_score_recalibration_file,
                referencefile_path => $referencefile_path,
                stderrfile_path    => $stderrfile_path,
                temp_directory     => $temp_directory,
                verbosity          => $active_parameter_href->{gatk_logging_level},
                xargs_mode         => 1,
            }
        );
        say {$xargsfilehandle} $NEWLINE;
    }

    ## GATK GatherBQSRReports
    say {$filehandle} q{## GATK GatherBQSRReports};
    my $gatk_gatherbqsr_outfile_path =
      $outfile_path_prefix . $DOT . $sample_id . $DOT . q{grp};
    gatk_gatherbqsrreports(
        {
            base_quality_score_recalibration_files_ref =>
              \@base_quality_score_recalibration_files,
            filehandle   => $filehandle,
            outfile_path => $gatk_gatherbqsr_outfile_path,
        }
    );
    say {$filehandle} $NEWLINE;

    ## GATK ApplyBQSR
    say {$filehandle} q{## GATK ApplyBQSR};

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            core_number        => $parallel_processes,
            filehandle         => $filehandle,
            file_path          => $recipe_file_path,
            recipe_info_path   => $recipe_info_path,
            xargsfilehandle    => $xargsfilehandle,
            xargs_file_counter => $xargs_file_counter,
        }
    );

  CONTIG:
    foreach my $contig ( @{ $file_info_href->{bam_contigs_size_ordered} } ) {

        my $stderrfile_path =
          $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};
        gatk_applybqsr(
            {
                base_quality_score_recalibration_file => $gatk_gatherbqsr_outfile_path,
                filehandle                            => $xargsfilehandle,
                infile_path                           => $infile_path{$contig},
                intervals_ref                         => $gatk_intervals{$contig},
                java_use_large_pages => $active_parameter_href->{java_use_large_pages},
                memory_allocation    => q{Xmx} . $JAVA_MEMORY_ALLOCATION . q{g},
                verbosity            => $active_parameter_href->{gatk_logging_level},
                read_filters_ref =>
                  \@{ $active_parameter_href->{gatk_baserecalibration_read_filters} },
                referencefile_path         => $referencefile_path,
                static_quantized_quals_ref => \@{
                    $active_parameter_href
                      ->{gatk_baserecalibration_static_quantized_quals}
                },
                outfile_path       => $outfile_path{$contig},
                referencefile_path => $referencefile_path,
                stderrfile_path    => $stderrfile_path,
                temp_directory     => $temp_directory,
                xargs_mode         => 1,
            }
        );
        say {$xargsfilehandle} $NEWLINE;
    }

    ## Gather BAM files
    say {$filehandle} q{## Gather BAM files};

    ## Assemble infile paths in contig order and not per size
    my @gather_infile_paths =
      map { $outfile_path{$_} } @{ $file_info_href->{bam_contigs} };
    my $store_outfile_path = $outfile_path_prefix . $outfile_suffix;

    picardtools_gatherbamfiles(
        {
            create_index     => q{true},
            filehandle       => $filehandle,
            infile_paths_ref => \@gather_infile_paths,
            java_jar =>
              catfile( $active_parameter_href->{picardtools_path}, q{picard.jar} ),
            java_use_large_pages => $active_parameter_href->{java_use_large_pages},
            memory_allocation    => q{Xmx4g},
            outfile_path         => $outfile_path_prefix . $outfile_suffix,
            referencefile_path   => $referencefile_path,
            temp_directory       => $temp_directory,
        }
    );
    say {$filehandle} $NEWLINE;

    ## Rename the bam file index file so that Expansion Hunter can find it
    say {$filehandle}
      q{## Copy index file to ".bam.bai" so that Expansionhunter can find it downstream};

    gnu_cp(
        {
            filehandle   => $filehandle,
            force        => 1,
            infile_path  => $outfile_path_prefix . q{.bai},
            outfile_path => $outfile_path_prefix . $outfile_suffix . q{.bai},
        }
    );
    say {$filehandle} $NEWLINE;

    ## Create BAM to CRAM for long term storage

    $store_outfile_path = $outfile_path_prefix . $DOT . q{cram};

    say {$filehandle} q{## Convert BAM to CRAM for long term storage};
    samtools_view(
        {
            filehandle         => $filehandle,
            infile_path        => $outfile_path_prefix . $outfile_suffix,
            outfile_path       => $store_outfile_path,
            output_format      => q{cram},
            referencefile_path => $referencefile_path,
            thread_number      => $parallel_processes,
        }
    );
    say {$filehandle} $NEWLINE;

    ## Index CRAM
    samtools_index(
        {
            filehandle  => $filehandle,
            infile_path => $store_outfile_path,
        }
    );

    close $xargsfilehandle;
    close $filehandle;

    if ( $recipe_mode == 1 ) {

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                infile           => $outfile_name_prefix,
                path             => $store_outfile_path,
                recipe_name      => $recipe_name,
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        set_file_path_to_store(
            {
                format           => q{cram},
                id               => $sample_id,
                path             => $store_outfile_path,
                path_index       => $store_outfile_path . $DOT . q{crai},
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
            }
        );

        submit_recipe(
            {
                base_command         => $profile_base_command,
                case_id              => $case_id,
                dependency_method    => q{sample_to_sample},
                job_id_chain         => $job_id_chain,
                job_id_href          => $job_id_href,
                job_reservation_name => $active_parameter_href->{job_reservation_name},
                log                  => $log,
                max_parallel_processes_count_href =>
                  $file_info_href->{max_parallel_processes_count},
                recipe_file_path   => $recipe_file_path,
                sample_id          => $sample_id,
                submission_profile => $active_parameter_href->{submission_profile},
            }
        );
    }
    return 1;
}

sub analysis_gatk_baserecalibration_panel {

## Function : GATK baserecalibrator/GatherBQSRReports/ApplyBQSR to recalibrate bases before variant calling. BaseRecalibrator/GatherBQSRReports/ApplyBQSR will be executed within the same sbatch script.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File info hash {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Program name
##          : $sample_id               => Sample id
##          : $sample_info_href        => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_info_href;
    my $sample_id;

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
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::File qw{ get_exome_target_bed_file get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_attributes get_recipe_resources };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Gatk qw{ gatk_applybqsr gatk_baserecalibrator };
    use MIP::Program::Samtools qw{ samtools_index samtools_view };
    use MIP::Sample_info qw{ set_file_path_to_store set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Unpack parameters
    ## Get the io infiles per chain and id
    my %io = get_io_files(
        {
            id             => $sample_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            stream         => q{in},
        }
    );
    my $infile_name_prefix = $io{in}{file_name_prefix};
    my $infile_path        = $io{in}{file_path};

    my %rec_atr = get_recipe_attributes(
        {
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
        }
    );
    my $job_id_chain       = $rec_atr{chain};
    my $recipe_mode        = $active_parameter_href->{$recipe_name};
    my $referencefile_path = $active_parameter_href->{human_genome_reference};
    my %recipe_resource    = get_recipe_resources(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );
    my $core_number = $recipe_resource{core_number};

    ## Outpaths
    ## Set and get the io files per chain, id and stream
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id               => $job_id_chain,
                id                     => $sample_id,
                file_info_href         => $file_info_href,
                file_name_prefixes_ref => [$infile_name_prefix],
                outdata_dir            => $active_parameter_href->{outdata_dir},
                parameter_href         => $parameter_href,
                recipe_name            => $recipe_name,
            }
        )
    );

    my $outfile_name_prefix = $io{out}{file_name_prefix};
    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my $outfile_path        = $io{out}{file_path};

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            directory_id                    => $sample_id,
            filehandle                      => $filehandle,
            job_id_href                     => $job_id_href,
            memory_allocation               => $recipe_resource{memory},
            process_time                    => $recipe_resource{time},
            recipe_directory                => $recipe_name,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => $recipe_resource{load_env_ref},
            temp_directory                  => $active_parameter_href->{temp_directory},
        }
    );

    ### SHELL:

    my $process_memory_allocation = $JAVA_MEMORY_ALLOCATION + $JAVA_GUEST_OS_MEMORY;

    # Constrain parallelization to match available memory
    my $parallel_processes = get_parallel_processes(
        {
            core_number               => $core_number,
            process_memory_allocation => $process_memory_allocation,
            recipe_memory_allocation  => $recipe_resource{memory},
        }
    );

    ## GATK BaseRecalibrator
    say {$filehandle} q{## GATK BaseRecalibrator};

    ## Get exome_target_bed file for specfic sample_id and add file_ending from file_info hash
    my $exome_target_bed_file = get_exome_target_bed_file(
        {
            exome_target_bed_href => $active_parameter_href->{exome_target_bed},
            sample_id             => $sample_id,
        }
    );
    my $padded_interval_list_ending = $file_info_href->{exome_target_bed}[1];
    my $padded_exome_target_bed_file =
      $exome_target_bed_file . $padded_interval_list_ending;

    my $base_quality_score_recalibration_file = $outfile_path_prefix . $DOT . q{grp};
    gatk_baserecalibrator(
        {
            filehandle           => $filehandle,
            infile_path          => $infile_path,
            intervals_ref        => [$padded_exome_target_bed_file],
            java_use_large_pages => $active_parameter_href->{java_use_large_pages},
            known_sites_ref =>
              \@{ $active_parameter_href->{gatk_baserecalibration_known_sites} },
            memory_allocation  => q{Xmx} . $JAVA_MEMORY_ALLOCATION . q{g},
            outfile_path       => $base_quality_score_recalibration_file,
            referencefile_path => $referencefile_path,
            temp_directory     => $active_parameter_href->{temp_directory},
            verbosity          => $active_parameter_href->{gatk_logging_level},
        }
    );
    say {$filehandle} $NEWLINE;

    ## GATK ApplyBQSR
    say {$filehandle} q{## GATK ApplyBQSR};

    gatk_applybqsr(
        {
            base_quality_score_recalibration_file =>
              $base_quality_score_recalibration_file,
            filehandle           => $filehandle,
            infile_path          => $infile_path,
            intervals_ref        => [$padded_exome_target_bed_file],
            java_use_large_pages => $active_parameter_href->{java_use_large_pages},
            memory_allocation    => q{Xmx} . $JAVA_MEMORY_ALLOCATION . q{g},
            verbosity            => $active_parameter_href->{gatk_logging_level},
            read_filters_ref =>
              \@{ $active_parameter_href->{gatk_baserecalibration_read_filters} },
            referencefile_path => $referencefile_path,
            static_quantized_quals_ref =>
              \@{ $active_parameter_href->{gatk_baserecalibration_static_quantized_quals}
              },
            outfile_path       => $outfile_path,
            referencefile_path => $referencefile_path,
            temp_directory     => $active_parameter_href->{temp_directory},
        }
    );
    say {$filehandle} $NEWLINE;

    ## Create BAM to CRAM for long term storage
    my $store_outfile_path = $outfile_path_prefix . $DOT . q{cram};

    say {$filehandle} q{## Convert BAM to CRAM for long term storage};
    samtools_view(
        {
            filehandle         => $filehandle,
            infile_path        => $outfile_path,
            outfile_path       => $store_outfile_path,
            output_format      => q{cram},
            referencefile_path => $referencefile_path,
            thread_number      => $parallel_processes,
        }
    );
    say {$filehandle} $NEWLINE;

    ## Index CRAM
    samtools_index(
        {
            filehandle  => $filehandle,
            infile_path => $store_outfile_path,
        }
    );

    close $filehandle;

    if ( $recipe_mode == 1 ) {

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                infile           => $outfile_name_prefix,
                path             => $store_outfile_path,
                recipe_name      => $recipe_name,
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        set_file_path_to_store(
            {
                format           => q{cram},
                id               => $sample_id,
                path             => $store_outfile_path,
                path_index       => $store_outfile_path . $DOT . q{crai},
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
            }
        );

        submit_recipe(
            {
                base_command         => $profile_base_command,
                case_id              => $case_id,
                dependency_method    => q{sample_to_sample},
                job_id_chain         => $job_id_chain,
                job_id_href          => $job_id_href,
                job_reservation_name => $active_parameter_href->{job_reservation_name},
                log                  => $log,
                max_parallel_processes_count_href =>
                  $file_info_href->{max_parallel_processes_count},
                recipe_file_path   => $recipe_file_path,
                sample_id          => $sample_id,
                submission_profile => $active_parameter_href->{submission_profile},
            }
        );
    }
    return 1;
}

sub analysis_gatk_baserecalibration_rna {

## Function : GATK baserecalibrator/GatherBQSRReports/ApplyBQSR to recalibrate bases before variant calling.
##          : BaseRecalibrator/GatherBQSRReports/ApplyBQSR will be executed within the same sbatch script.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File info hash {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Program name
##          : $sample_id               => Sample id
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_info_href;
    my $sample_id;

    ## Default(s)
    my $case_id;
    my $profile_base_command;
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
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            store       => \$temp_directory,
            strict_type => 1,
        },
        xargs_file_counter => {
            allow       => qr/ \A \d+ \z /xsm,
            default     => 0,
            store       => \$xargs_file_counter,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Cluster qw{ get_parallel_processes };
    use MIP::Get::File qw{ get_merged_infile_prefix get_io_files };
    use MIP::Get::Parameter
      qw{ get_gatk_intervals get_recipe_attributes get_recipe_resources };
    use MIP::Program::Gnu::Coreutils qw{ gnu_cp };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Gatk
      qw{ gatk_applybqsr gatk_baserecalibrator gatk_gatherbqsrreports };
    use MIP::Program::Picardtools qw{ picardtools_gatherbamfiles };
    use MIP::Program::Samtools qw{ samtools_index samtools_view };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Sample_info qw{ set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Unpack parameters
    ## Get the io infiles per chain and id
    my %io = get_io_files(
        {
            id             => $sample_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            stream         => q{in},
        }
    );
    my $indir_path_prefix  = $io{in}{dir_path_prefix};
    my $infile_suffix      = $io{in}{file_suffix};
    my $infile_name_prefix = $io{in}{file_name_prefix};
    my %infile_path        = %{ $io{in}{file_path_href} };

    my %rec_atr = get_recipe_attributes(
        {
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
        }
    );
    my $analysis_type      = $active_parameter_href->{analysis_type}{$sample_id};
    my $job_id_chain       = $rec_atr{chain};
    my $recipe_mode        = $active_parameter_href->{$recipe_name};
    my $referencefile_path = $active_parameter_href->{human_genome_reference};
    my $xargs_file_path_prefix;
    my %recipe_resource = get_recipe_resources(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );
    my $core_number = $recipe_resource{core_number};

    ## Add merged infile name prefix after merging all BAM files per sample_id
    my $merged_infile_prefix = get_merged_infile_prefix(
        {
            file_info_href => $file_info_href,
            sample_id      => $sample_id,
        }
    );

    ## Outpaths
    ## Assign suffix
    my $outfile_suffix = $rec_atr{outfile_suffix};
    my $outsample_directory =
      catdir( $active_parameter_href->{outdata_dir}, $sample_id, $recipe_name );
    my $outfile_tag =
      $file_info_href->{$sample_id}{$recipe_name}{file_tag};
    my @outfile_paths =
      map {
        catdir( $outsample_directory,
            $merged_infile_prefix . $outfile_tag . $DOT . $_ . $outfile_suffix )
      } @{ $file_info_href->{bam_contigs_size_ordered} };

    ## Set and get the io files per chain, id and stream
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id       => $job_id_chain,
                id             => $sample_id,
                file_info_href => $file_info_href,
                file_paths_ref => \@outfile_paths,
                parameter_href => $parameter_href,
                recipe_name    => $recipe_name,
                temp_directory => $temp_directory,
            }
        )
    );

    my $outdir_path_prefix  = $io{out}{dir_path_prefix};
    my $outfile_name_prefix = $io{out}{file_name_prefix};
    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my %outfile_path        = %{ $io{out}{file_path_href} };

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle      = IO::Handle->new();
    my $xargsfilehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            directory_id                    => $sample_id,
            filehandle                      => $filehandle,
            job_id_href                     => $job_id_href,
            memory_allocation               => $recipe_resource{memory},
            process_time                    => $recipe_resource{time},
            recipe_directory                => $recipe_name,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => $recipe_resource{load_env_ref},
            temp_directory                  => $temp_directory,
        }
    );

    ### SHELL:

    ## Generate gatk intervals. Chromosomes for WGS/WTS and paths to contig_bed_files for WES
    my %gatk_intervals = get_gatk_intervals(
        {
            analysis_type         => $analysis_type,
            contigs_ref           => \@{ $file_info_href->{bam_contigs_size_ordered} },
            exome_target_bed_href => $active_parameter_href->{exome_target_bed},
            filehandle            => $filehandle,
            file_ending           => $file_info_href->{exome_target_bed}[0],
            max_cores_per_node    => $core_number,
            log                   => $log,
            outdirectory          => $outdir_path_prefix,
            reference_dir         => $active_parameter_href->{reference_dir},
            sample_id             => $sample_id,
        }
    );

    my $process_memory_allocation = $JAVA_MEMORY_ALLOCATION + $JAVA_GUEST_OS_MEMORY;

    # Constrain parallelization to match available memory
    my $parallel_processes = get_parallel_processes(
        {
            core_number               => $core_number,
            process_memory_allocation => $process_memory_allocation,
            recipe_memory_allocation  => $recipe_resource{memory},
        }
    );

    ## GATK BaseRecalibrator
    say {$filehandle} q{## GATK BaseRecalibrator};

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            core_number        => $parallel_processes,
            filehandle         => $filehandle,
            file_path          => $recipe_file_path,
            recipe_info_path   => $recipe_info_path,
            xargsfilehandle    => $xargsfilehandle,
            xargs_file_counter => $xargs_file_counter,
        }
    );

    my @base_quality_score_recalibration_files;
  CONTIG:
    foreach my $contig ( @{ $file_info_href->{bam_contigs_size_ordered} } ) {

        my $base_quality_score_recalibration_file =
          $outfile_path_prefix . $DOT . $contig . $DOT . q{grp};

        ## Add for gathering base recal files later
        push @base_quality_score_recalibration_files,
          $base_quality_score_recalibration_file;
        my $stderrfile_path =
          $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};
        gatk_baserecalibrator(
            {
                filehandle           => $xargsfilehandle,
                infile_path          => $infile_path{$contig},
                intervals_ref        => $gatk_intervals{$contig},
                java_use_large_pages => $active_parameter_href->{java_use_large_pages},
                known_sites_ref =>
                  \@{ $active_parameter_href->{gatk_baserecalibration_known_sites} },
                memory_allocation  => q{Xmx} . $JAVA_MEMORY_ALLOCATION . q{g},
                outfile_path       => $base_quality_score_recalibration_file,
                referencefile_path => $referencefile_path,
                stderrfile_path    => $stderrfile_path,
                temp_directory     => $temp_directory,
                verbosity          => $active_parameter_href->{gatk_logging_level},
                xargs_mode         => 1,
            }
        );
        say {$xargsfilehandle} $NEWLINE;
    }

    ## GATK GatherBQSRReports
    say {$filehandle} q{## GATK GatherBQSRReports};
    my $gatk_gatherbqsr_outfile_path =
      $outfile_path_prefix . $DOT . $sample_id . $DOT . q{grp};
    gatk_gatherbqsrreports(
        {
            base_quality_score_recalibration_files_ref =>
              \@base_quality_score_recalibration_files,
            filehandle   => $filehandle,
            outfile_path => $gatk_gatherbqsr_outfile_path,
        }
    );
    say {$filehandle} $NEWLINE;

    ## GATK ApplyBQSR
    say {$filehandle} q{## GATK ApplyBQSR};

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            core_number        => $parallel_processes,
            filehandle         => $filehandle,
            file_path          => $recipe_file_path,
            recipe_info_path   => $recipe_info_path,
            xargsfilehandle    => $xargsfilehandle,
            xargs_file_counter => $xargs_file_counter,
        }
    );

  CONTIG:
    foreach my $contig ( @{ $file_info_href->{bam_contigs_size_ordered} } ) {

        my $stderrfile_path =
          $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};
        gatk_applybqsr(
            {
                base_quality_score_recalibration_file => $gatk_gatherbqsr_outfile_path,
                filehandle                            => $xargsfilehandle,
                infile_path                           => $infile_path{$contig},
                intervals_ref                         => $gatk_intervals{$contig},
                java_use_large_pages => $active_parameter_href->{java_use_large_pages},
                memory_allocation    => q{Xmx} . $JAVA_MEMORY_ALLOCATION . q{g},
                verbosity            => $active_parameter_href->{gatk_logging_level},
                read_filters_ref =>
                  \@{ $active_parameter_href->{gatk_baserecalibration_read_filters} },
                referencefile_path         => $referencefile_path,
                static_quantized_quals_ref => \@{
                    $active_parameter_href
                      ->{gatk_baserecalibration_static_quantized_quals}
                },
                outfile_path       => $outfile_path{$contig},
                referencefile_path => $referencefile_path,
                stderrfile_path    => $stderrfile_path,
                temp_directory     => $temp_directory,
                xargs_mode         => 1,
            }
        );
        say {$xargsfilehandle} $NEWLINE;
    }

    ## Gather BAM files
    say {$filehandle} q{## Gather BAM files};

    ## Assemble infile paths in contig order and not per size
    my @gather_infile_paths =
      map { $outfile_path{$_} } @{ $file_info_href->{bam_contigs} };
    my $store_outfile_path = $outfile_path_prefix . $outfile_suffix;

    picardtools_gatherbamfiles(
        {
            create_index     => q{true},
            filehandle       => $filehandle,
            infile_paths_ref => \@gather_infile_paths,
            java_jar =>
              catfile( $active_parameter_href->{picardtools_path}, q{picard.jar} ),
            java_use_large_pages => $active_parameter_href->{java_use_large_pages},
            memory_allocation    => q{Xmx4g},
            outfile_path         => $outfile_path_prefix . $outfile_suffix,
            referencefile_path   => $referencefile_path,
            temp_directory       => $temp_directory,
        }
    );
    say {$filehandle} $NEWLINE;

    close $xargsfilehandle;
    close $filehandle;

    if ( $recipe_mode == 1 ) {

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                infile           => $outfile_name_prefix,
                path             => $store_outfile_path,
                recipe_name      => $recipe_name,
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        submit_recipe(
            {
                base_command         => $profile_base_command,
                case_id              => $case_id,
                dependency_method    => q{sample_to_sample},
                job_id_chain         => $job_id_chain,
                job_id_href          => $job_id_href,
                job_reservation_name => $active_parameter_href->{job_reservation_name},
                log                  => $log,
                max_parallel_processes_count_href =>
                  $file_info_href->{max_parallel_processes_count},
                recipe_file_path   => $recipe_file_path,
                sample_id          => $sample_id,
                submission_profile => $active_parameter_href->{submission_profile},
            }
        );
    }
    return 1;
}

1;
