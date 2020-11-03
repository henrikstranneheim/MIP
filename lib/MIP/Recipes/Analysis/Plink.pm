package MIP::Recipes::Analysis::Plink;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile splitpath};
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $ASTERISK $DASH $DOT $LOG_NAME $NEWLINE $PIPE $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.22;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_plink };

}

## Constants
Readonly my $CHR_X_NUMBER        => 23;
Readonly my $CHR_Y_NUMBER        => 24;
Readonly my $FEMALE_MAX_F        => 0.2;
Readonly my $INDEP_WINDOW_SIZE   => 50;
Readonly my $INDEP_STEP_SIZE     => 5;
Readonly my $INDEP_VIF_THRESHOLD => 2;
Readonly my $MALE_MIN_F          => 0.75;

sub analysis_plink {

## Function : Tests sample for correct relatives (only performed for samples with relatives defined in pedigree file) performed on sequence data.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Program name
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory

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
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Pedigree qw{ create_fam_file };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_attributes get_recipe_resources };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Bcftools qw(bcftools_view bcftools_annotate);
    use MIP::Program::Plink
      qw{ plink_calculate_inbreeding plink_check_sex_chroms plink_create_mibs plink_fix_fam_ped_map_freq plink_sex_check plink_variant_pruning };
    use MIP::Program::Vt qw(vt_uniq);
    use MIP::Sample_info
      qw{ set_recipe_outfile_in_sample_info set_recipe_metafile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

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
    my $infile_path_prefix = $io{in}{file_path_prefix};
    my $infile_suffix      = $io{in}{file_suffix};
    my $infile_path        = $io{in}{file_path};

    my $consensus_analysis_type = $parameter_href->{cache}{consensus_analysis_type};
    my $human_genome_reference_version =
      $file_info_href->{human_genome_reference_version};
    my $human_genome_reference_source = $file_info_href->{human_genome_reference_source};
    my $job_id_chain                  = get_recipe_attributes(
        {
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            attribute      => q{chain},
        }
    );
    my $recipe_mode     = $active_parameter_href->{$recipe_name};
    my @sample_ids      = @{ $active_parameter_href->{sample_ids} };
    my %recipe_resource = get_recipe_resources(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    ## Set outfiles depending on sample data
    my %plink_outanalysis_prefix;
    my %plink_outfile_analysis_map = (
        check_for_sex    => { plink_sexcheck => q{sexcheck}, },
        multiple_samples => {
            inbreeding_factor => q{het},
            relation_check    => q{mibs},
        },
    );

    my @plink_outfiles;
  MODE:
    while ( my ( $mode, $program_href ) = each %plink_outfile_analysis_map ) {

      PLINK_PROGRAM:
        while ( my ( $file_name_prefix, $file_suffix ) = each %{$program_href} ) {

            if ( @sample_ids > 1 and $mode eq q{multiple_samples} ) {

                $plink_outanalysis_prefix{$file_name_prefix} = $file_name_prefix;
                push @plink_outfiles, $file_name_prefix . $DOT . $file_suffix;
                next;
            }
            if (    defined $active_parameter_href->{gender}{others}
                and @{ $active_parameter_href->{gender}{others} } != @sample_ids
                and $mode eq q{check_for_sex} )
            {
                $plink_outanalysis_prefix{$file_name_prefix} = $file_name_prefix;
                push @plink_outfiles, $file_name_prefix . $DOT . $file_suffix;
            }
        }
    }

    ## No eligible test to run
    if ( not @plink_outfiles ) {
        $log->warn(
            q{No eligible Plink test to run for pedigree and sample(s) - skipping 'plink'}
        );
        return;
    }

    ## Set and get the io files per chain, id and stream
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id         => $job_id_chain,
                id               => $case_id,
                file_info_href   => $file_info_href,
                file_name_prefix => $infile_name_prefix,
                iterators_ref    => \@plink_outfiles,
                outdata_dir      => $active_parameter_href->{outdata_dir},
                parameter_href   => $parameter_href,
                recipe_name      => $recipe_name,
                temp_directory   => $temp_directory,
            }
        )
    );

    my $outdir_path_prefix  = $io{out}{dir_path_prefix};
    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my %outfile_path        = %{ $io{out}{file_path_href} };

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $recipe_resource{core_number},
            directory_id                    => $case_id,
            filehandle                      => $filehandle,
            job_id_href                     => $job_id_href,
            memory_allocation               => $recipe_resource{memory},
            process_time                    => $recipe_resource{time},
            recipe_directory                => $recipe_name,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => $recipe_resource{load_env_ref},
        }
    );

    # Split to enable submission to &sample_info_qc later
    my ( $volume, $directory, $recipe_info_file ) = splitpath($recipe_info_path);

    # To enable submission to %sample_info_qc later
    my $stdout_file_path = catfile( $directory, $recipe_info_file . q{.stdout.txt} );

    ### SHELL:

    my $plink_outfile_prefix = catfile($outfile_path_prefix);

    my $case_file_path = catfile( $outdir_path_prefix, $case_id . $DOT . q{fam} );

    ## Create .fam file to be used in plink analyses
    create_fam_file(
        {
            case_id          => $case_id,
            fam_file_path    => $case_file_path,
            filehandle       => $filehandle,
            parameter_href   => $parameter_href,
            sample_ids_ref   => $active_parameter_href->{sample_ids},
            sample_info_href => $sample_info_href,
        }
    );

    ## Prepare input

    ## Define regions to analyse
    my $regions_file_path = _setup_plink_for_analysis_type(
        {
            consensus_analysis_type => $consensus_analysis_type,
            exome_target_bed_href   => $active_parameter_href->{exome_target_bed},
            filehandle              => $filehandle,
            outdir_path_prefix      => $outdir_path_prefix,
            sample_ids_ref          => \@sample_ids,
        }
    );

    say {$filehandle} q{## Remove indels using bcftools};
    my $view_outfile_path =
      $infile_path_prefix . $UNDERSCORE . q{no_indels} . $infile_suffix;
    bcftools_view(
        {
            exclude_types_ref => [qw{indels}],
            filehandle        => $filehandle,
            infile_path       => $infile_path,
            outfile_path      => $view_outfile_path,
            output_type       => q{v},
            regions_file_path => $regions_file_path,
        }
    );
    say {$filehandle} $NEWLINE;

    say {$filehandle} q{## Create uniq IDs and remove duplicate variants};
    bcftools_annotate(
        {
            filehandle     => $filehandle,
            infile_path    => $view_outfile_path,
            output_type    => q{v},
            remove_ids_ref => [q{ID}],
            set_id         => q?+'%CHROM:%POS:%REF:%ALT'?,
        }
    );

    print {$filehandle} $PIPE . $SPACE;

    ## Drops duplicate variants that appear later in the the VCF file
    my $uniq_outfile_path =
      $infile_path_prefix . $UNDERSCORE . q{no_indels_ann_uniq} . $infile_suffix;
    vt_uniq(
        {
            filehandle   => $filehandle,
            infile_path  => $DASH,
            outfile_path => $uniq_outfile_path,
        }
    );
    say {$filehandle} $NEWLINE;

    ### Plink variant pruning and creation of unique Ids
    say {$filehandle} q{## Create pruning set and uniq IDs};
    plink_variant_pruning(
        {
            const_fid           => $case_id,
            filehandle          => $filehandle,
            make_bed            => 1,
            outfile_prefix      => $plink_outfile_prefix,
            indep               => 1,
            indep_step_size     => $INDEP_STEP_SIZE,
            indep_vif_threshold => $INDEP_VIF_THRESHOLD,
            indep_window_size   => $INDEP_WINDOW_SIZE,
            set_missing_var_ids => q?@:#[?
              . $human_genome_reference_version
              . q?]\$1,\$2?,
            vcffile_path   => $uniq_outfile_path,
            vcf_half_call  => q{haploid},
            vcf_require_gt => 1,
        }
    );
    say {$filehandle} $NEWLINE;

    say {$filehandle}
      q{## Update Plink fam. Create ped and map file and frequency report};
    ## Get parameters
    my $allow_no_sex;

    ## If not all samples have a known sex
    if ( $active_parameter_href->{gender}{others}
        and @{ $active_parameter_href->{gender}{others} } )
    {

        $allow_no_sex = 1;
    }
    my $binary_fileset_prefix = $plink_outfile_prefix;

    plink_fix_fam_ped_map_freq(
        {
            allow_no_sex          => $allow_no_sex,
            binary_fileset_prefix => $binary_fileset_prefix,
            fam_file_path         => $case_file_path,
            filehandle            => $filehandle,
            freqx                 => 1,
            make_just_fam         => 1,
            outfile_prefix        => $plink_outfile_prefix,
            recode                => 1,
        }
    );
    say {$filehandle} $NEWLINE;

    # Only perform if more than 1 sample
    if ( @sample_ids > 1 ) {

        my $inbreeding_outfile_prefix_hets =
          $binary_fileset_prefix . $DOT . $plink_outanalysis_prefix{inbreeding_factor};

        say {$filehandle} q{## Calculate inbreeding coefficients per case};
        plink_calculate_inbreeding(
            {
                binary_fileset_prefix   => $binary_fileset_prefix,
                extract_file            => $binary_fileset_prefix . $DOT . q{prune.in},
                filehandle              => $filehandle,
                inbreeding_coefficients => 1,
                het                     => 1,
                outfile_prefix          => $inbreeding_outfile_prefix_hets,
                small_sample            => 1,
            }
        );
        say {$filehandle} $NEWLINE;

        my $inbreeding_outfile_prefix_mibs =
          $binary_fileset_prefix . $DOT . $plink_outanalysis_prefix{relation_check};
        say {$filehandle} q{## Create Plink .mibs per case};
        plink_create_mibs(
            {
                cluster        => 1,
                filehandle     => $filehandle,
                map_file_path  => $binary_fileset_prefix . $DOT . q{map},
                matrix         => 1,
                outfile_prefix => $inbreeding_outfile_prefix_mibs,
                ped_file_path  => $binary_fileset_prefix . $DOT . q{ped},
            }
        );
        say {$filehandle} $NEWLINE;
    }

    if ( defined $active_parameter_href->{gender}{others}
        and @{ $active_parameter_href->{gender}{others} } != @sample_ids )
    {

        ## Only if not all samples have unknown sex
        ### Plink sex-check
        ## Get parameters
        my $genome_build;

        ## Set correct build prefix
        if ( $human_genome_reference_source eq q{grch} ) {

            $genome_build = q{b} . $human_genome_reference_version;
        }
        else {

            $genome_build = q{hg} . $human_genome_reference_version;
        }

        plink_check_sex_chroms(
            {
                binary_fileset_prefix => $binary_fileset_prefix,
                filehandle            => $filehandle,
                no_fail               => 1,
                make_bed              => 1,
                outfile_prefix        => $plink_outfile_prefix . $UNDERSCORE . q{unsplit},
                regions_ref           => [ $CHR_X_NUMBER, $CHR_Y_NUMBER ],
                split_x               => $genome_build,
            }
        );
        say {$filehandle} $NEWLINE;

        ## Get parameters
        my $sex_check_min_f;
        if ( $consensus_analysis_type eq q{wes} ) {

            $sex_check_min_f = $FEMALE_MAX_F . $SPACE . $MALE_MIN_F;
        }
        my $extract_file;
        my $read_freqfile_path;

        if ( @sample_ids > 1 ) {

            $extract_file       = $binary_fileset_prefix . $DOT . q{prune.in};
            $read_freqfile_path = $binary_fileset_prefix . $DOT . q{frqx};
        }

        my $sex_check_outfile_prefix =
          $binary_fileset_prefix . $DOT . $plink_outanalysis_prefix{plink_sexcheck};
        plink_sex_check(
            {
                binary_fileset_prefix => $binary_fileset_prefix
                  . $UNDERSCORE
                  . q{unsplit},
                extract_file       => $extract_file,
                filehandle         => $filehandle,
                outfile_prefix     => $sex_check_outfile_prefix,
                read_freqfile_path => $read_freqfile_path,
                sex_check_min_f    => $sex_check_min_f,
            }
        );
        say {$filehandle} $NEWLINE;
    }

    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe_mode == 1 ) {

        while ( my ( $outfile_tag, $outfile_path ) = each %outfile_path ) {

            ## Collect QC metadata info for later use
            set_recipe_outfile_in_sample_info(
                {
                    path             => $outfile_path,
                    recipe_name      => $outfile_tag,
                    sample_info_href => $sample_info_href,
                }
            );

        }

        submit_recipe(
            {
                base_command         => $profile_base_command,
                case_id              => $case_id,
                dependency_method    => q{case_to_island},
                job_id_chain         => $job_id_chain,
                job_id_href          => $job_id_href,
                job_reservation_name => $active_parameter_href->{job_reservation_name},
                log                  => $log,
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

sub _setup_plink_for_analysis_type {

## Function : Restrict analysis to capture kit depending on analysis type.
## Returns  :
## Arguments: $consensus_analysis_type => Consensus analysis type
##          : $exome_target_bed_href   => Bed files {REF}
##          : $filehandle              => Filehandle
##          : $outdir_path_prefix      => Outdir path prefix
##          : $sample_ids_ref          => Sample ID:s {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $consensus_analysis_type;
    my $exome_target_bed_href;
    my $filehandle;
    my $outdir_path_prefix;
    my $sample_ids_ref;

    my $tmpl = {
        consensus_analysis_type => {
            defined     => 1,
            required    => 1,
            store       => \$consensus_analysis_type,
            strict_type => 1,
        },
        exome_target_bed_href => {
            default     => {},
            required    => 1,
            store       => \$exome_target_bed_href,
            strict_type => 1,
        },
        filehandle => {
            required    => 1,
            store       => \$filehandle,
        },
        outdir_path_prefix => {
            defined     => 1,
            required    => 1,
            store       => \$outdir_path_prefix,
            strict_type => 1,
        },
        sample_ids_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$sample_ids_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use List::MoreUtils qw{ uniq };
    use MIP::Program::Bedtools qw{ bedtools_intersect };
    use MIP::Program::Gnu::Coreutils qw{ gnu_mv };
    use MIP::Get::File qw{ get_exome_target_bed_file };

    ## Single sample case or all samples are wgs
    return if ( @{$sample_ids_ref} == 1 or $consensus_analysis_type eq q{wgs} );

    my @bed_file_paths;

    SAMPLE_ID:
    foreach my $sample_id ( @{$sample_ids_ref} ) {

        push @bed_file_paths, get_exome_target_bed_file(
            {
                exome_target_bed_href => $exome_target_bed_href,
                sample_id             => $sample_id,
            }
        );
    }
    @bed_file_paths = uniq(@bed_file_paths);

    ## All capture kits are the same
    return if ( @bed_file_paths == 1 );

    my $intersect_bed_file_path = catfile( $outdir_path_prefix, q{intersect.bed} );
    my $infile_path             = pop @bed_file_paths;
    BED_FILE_PATH:
    foreach my $bed_file_path (@bed_file_paths) {

        my $temp_intersect_path = catfile( $outdir_path_prefix, q{temp_intersect.bed} );
        bedtools_intersect(
            {
                filehandle         => $filehandle,
                infile_path        => $infile_path,
                intersectfile_path => $bed_file_path,
                stdoutfile_path    => $temp_intersect_path,
            }
        );
        print {$filehandle} $NEWLINE;
        gnu_mv(
            {
                filehandle   => $filehandle,
                infile_path  => $temp_intersect_path,
                outfile_path => $intersect_bed_file_path,
            }
        );
        print {$filehandle} $NEWLINE;

        ## Setup for next iteration
        $infile_path = $intersect_bed_file_path;
    }
    print {$filehandle} $NEWLINE;

    return $intersect_bed_file_path;
}

1;
