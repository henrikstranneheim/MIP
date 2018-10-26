package MIP::Recipes::Analysis::Gatk_variantrecalibration;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile splitpath };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use List::MoreUtils qw { uniq };
use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.06;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ analysis_gatk_variantrecalibration_wes analysis_gatk_variantrecalibration_wgs };

}

## Constants
Readonly my $ASTERISK   => q{*};
Readonly my $AMPERSAND  => q{&};
Readonly my $COLON      => q{:};
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

sub analysis_gatk_variantrecalibration_wes {

## Function : GATK VariantRecalibrator/ApplyRecalibration analysis recipe for wes data
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id               => Family id
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $recipe_name            => Program name
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

    use MIP::Delete::List qw{ delete_contig_elements };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_parameters get_recipe_attributes };
    use MIP::Gnu::Coreutils qw{ gnu_mv };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Variantcalling::Bcftools qw{ bcftools_norm };
    use MIP::Program::Variantcalling::Gatk
      qw{ gatk_variantrecalibrator gatk_applyvqsr gatk_selectvariants gatk_calculategenotypeposteriors };
    use MIP::QC::Record qw{ add_recipe_outfile_to_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Constants
    Readonly my $MAX_GAUSSIAN_LEVEL => 4;

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
    my $infile_path        = $io{in}{file_path};

    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};
    my $enable_indel_max_gaussians_filter =
      $active_parameter_href->{gatk_variantrecalibration_indel_max_gaussians};
    my $enable_snv_max_gaussians_filter =
      $active_parameter_href->{gatk_variantrecalibration_snv_max_gaussians};
    my $job_id_chain = get_recipe_attributes(
        {
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            attribute      => q{chain},
        }
    );
    my $recipe_mode        = $active_parameter_href->{$recipe_name};
    my $referencefile_path = $active_parameter_href->{human_genome_reference};
    my $resource_indel_href =
      $active_parameter_href->{gatk_variantrecalibration_resource_indel};
    my $resource_snv_href =
      $active_parameter_href->{gatk_variantrecalibration_resource_snv};
    my ( $core_number, $time, @source_environment_cmds ) = get_recipe_parameters(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id               => $job_id_chain,
                id                     => $case_id,
                file_info_href         => $file_info_href,
                file_name_prefixes_ref => [$infile_name_prefix],
                outdata_dir            => $active_parameter_href->{outdata_dir},
                parameter_href         => $parameter_href,
                recipe_name            => $recipe_name,
                temp_directory         => $temp_directory,
            }
        )
    );
    my $outdir_path_prefix       = $io{out}{dir_path_prefix};
    my $outfile_path_prefix      = $io{out}{file_path_prefix};
    my $outfile_suffix           = $io{out}{file_suffix};
    my $outfile_path             = $io{out}{file_path};
    my $temp_outfile_path_prefix = $io{temp}{file_path_prefix};

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

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

    ## Split to enable submission to &sample_info_qc later
    my ( $volume, $directory, $stderr_file ) =
      splitpath( $recipe_info_path . $DOT . q{stderr.txt} );

    ## Create .fam file to be used in variant calling analyses
    my $fam_file_path = catfile( $outdir_path_prefix, $case_id . $DOT . q{fam} );
    create_fam_file(
        {
            active_parameter_href => $active_parameter_href,
            execution_mode        => q{system},
            fam_file_path         => $fam_file_path,
            FILEHANDLE            => $FILEHANDLE,
            parameter_href        => $parameter_href,
            sample_info_href      => $sample_info_href,
        }
    );

    ## Check if "--pedigree" and "--pedigreeValidationType" should be included in analysis
    my %commands = gatk_pedigree_flag(
        {
            fam_file_path => $fam_file_path,
            recipe_name   => $recipe_name,
        }
    );

    ### GATK VariantRecalibrator
    ## Set mode to be used in variant recalibration
    # Exome will be processed using mode BOTH since there are to few INDELS
    # to use in the recalibration model even though using 30 exome BAMS in
    # Haplotypecaller step.
    my @modes = q{BOTH};

  MODE:
    foreach my $mode (@modes) {

        say {$FILEHANDLE} q{## GATK VariantRecalibrator};

        ## Get parameters
        my $max_gaussian_level;

        my @annotations =
          @{ $active_parameter_href->{gatk_variantrecalibration_annotations} };

        ### Special case: Not to be used with hybrid capture
        ## Removes an element from array and return new array while leaving orginal elements_ref untouched
        @annotations = delete_contig_elements(
            {
                elements_ref       => \@annotations,
                remove_contigs_ref => [qw{ DP }],
            }
        );
        my @snv_resources =
          _build_gatk_resource_command( { resources_href => $resource_snv_href, } );
        my @indel_resources =
          _build_gatk_resource_command( { resources_href => $resource_indel_href, } );

        # Create distinct set i.e. no duplicates.
        my @resources = uniq( @indel_resources, @snv_resources );

        ## Use hard filtering
        if (   $enable_snv_max_gaussians_filter
            || $enable_indel_max_gaussians_filter )
        {

            $max_gaussian_level = $MAX_GAUSSIAN_LEVEL;
        }

        my $recal_file_path = $temp_outfile_path_prefix . $DOT . q{intervals};
        gatk_variantrecalibrator(
            {
                annotations_ref      => \@annotations,
                FILEHANDLE           => $FILEHANDLE,
                infile_path          => $infile_path,
                java_use_large_pages => $active_parameter_href->{java_use_large_pages},
                max_gaussian_level   => $max_gaussian_level,
                memory_allocation    => q{Xmx10g},
                mode                 => $mode,
                outfile_path         => $recal_file_path,
                referencefile_path   => $referencefile_path,
                resources_ref        => \@resources,
                rscript_file_path    => $recal_file_path . $DOT . q{plots.R},
                temp_directory       => $temp_directory,
                tranches_file_path   => $recal_file_path . $DOT . q{tranches},
                verbosity            => $active_parameter_href->{gatk_logging_level},
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        ## GATK ApplyVQSR
        say {$FILEHANDLE} q{## GATK ApplyVQSR};

        ## Get parameters
        my $ts_filter_level;
        ## Exome analysis use combined reference for more power

        ## Infile genotypegvcfs combined vcf which used reference gVCFs to create combined vcf file
        $ts_filter_level =
          $active_parameter_href->{gatk_variantrecalibration_snv_tsfilter_level};

        gatk_applyvqsr(
            {
                FILEHANDLE           => $FILEHANDLE,
                infile_path          => $infile_path,
                java_use_large_pages => $active_parameter_href->{java_use_large_pages},
                memory_allocation    => q{Xmx10g},
                mode                 => $mode,
                outfile_path         => $outfile_path,
                recal_file_path      => $recal_file_path,
                referencefile_path   => $referencefile_path,
                temp_directory       => $temp_directory,
                tranches_file_path   => $recal_file_path . $DOT . q{tranches},
                ts_filter_level      => $ts_filter_level,
                verbosity            => $active_parameter_href->{gatk_logging_level},
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    ## BcfTools norm, Left-align and normalize indels, split multiallelics
    my $norm_outfile_path =
      $outfile_path_prefix . $UNDERSCORE . q{normalized} . $outfile_suffix;
    bcftools_norm(
        {
            FILEHANDLE      => $FILEHANDLE,
            infile_path     => $outfile_path,
            multiallelic    => q{-},
            outfile_path    => $norm_outfile_path,
            output_type     => q{v},
            reference_path  => $referencefile_path,
            stderrfile_path => $outfile_path_prefix . $UNDERSCORE . q{normalized.stderr},
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ### GATK SelectVariants

    ## Removes all genotype information for exome ref and recalulates meta-data info for remaining samples in new file.
    # Exome analysis

    say {$FILEHANDLE} q{## GATK SelectVariants};

    gatk_selectvariants(
        {
            FILEHANDLE           => $FILEHANDLE,
            exclude_non_variants => 1,
            infile_path          => $norm_outfile_path,
            java_use_large_pages => $active_parameter_href->{java_use_large_pages},
            memory_allocation    => q{Xmx2g},
            outfile_path         => $outfile_path,
            referencefile_path   => $referencefile_path,
            sample_names_ref     => \@{ $active_parameter_href->{sample_ids} },
            temp_directory       => $temp_directory,
            verbosity            => $active_parameter_href->{gatk_logging_level},
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Genotype refinement
    if ( $parameter_href->{dynamic_parameter}{trio} ) {

        say {$FILEHANDLE} q{## GATK CalculateGenotypePosteriors};

        my $calculategt_outfile_path =
          $outfile_path_prefix . $UNDERSCORE . q{refined} . $outfile_suffix;
        gatk_calculategenotypeposteriors(
            {
                FILEHANDLE           => $FILEHANDLE,
                infile_path          => $outfile_path,
                java_use_large_pages => $active_parameter_href->{java_use_large_pages},
                memory_allocation    => q{Xmx6g},
                outfile_path         => $calculategt_outfile_path,
                pedigree             => $commands{pedigree},
                referencefile_path   => $referencefile_path,
                supporting_callset_file_path =>
                  $active_parameter_href->{gatk_calculategenotypeposteriors_support_set},
                temp_directory => $temp_directory,
                verbosity      => $active_parameter_href->{gatk_logging_level},
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        ## Change name of file to accomodate downstream
        gnu_mv(
            {
                FILEHANDLE   => $FILEHANDLE,
                infile_path  => $calculategt_outfile_path,
                outfile_path => $outfile_path,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    ## BcfTools norm, Left-align and normalize indels, split multiallelics
    my $filtered_norm_outfile_path =
      $outfile_path_prefix . $UNDERSCORE . q{filtered_normalized} . $outfile_suffix;
    bcftools_norm(
        {
            FILEHANDLE     => $FILEHANDLE,
            infile_path    => $outfile_path,
            multiallelic   => q{-},
            outfile_path   => $filtered_norm_outfile_path,
            output_type    => q{v},
            reference_path => $referencefile_path,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Change name of file to accomodate downstream
    gnu_mv(
        {
            FILEHANDLE   => $FILEHANDLE,
            infile_path  => $filtered_norm_outfile_path,
            outfile_path => $outfile_path,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    close $FILEHANDLE;

    if ( $recipe_mode == 1 ) {

        ## Collect QC metadata info for later use
        add_recipe_outfile_to_sample_info(
            {
                path             => $outfile_path,
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
            }
        );

        # Used to find order of samples in qccollect downstream
        add_recipe_outfile_to_sample_info(
            {
                path             => $outfile_path,
                recipe_name      => q{pedigree_check},
                sample_info_href => $sample_info_href,
            }
        );

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

sub analysis_gatk_variantrecalibration_wgs {

## Function : GATK VariantRecalibrator/ApplyRecalibration analysis recipe for wgs data
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id               => Family id
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $recipe_name            => Program name
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

    use MIP::Delete::List qw{ delete_contig_elements };
    use MIP::File::Format::Pedigree qw{ create_fam_file gatk_pedigree_flag };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_parameters get_recipe_attributes };
    use MIP::Gnu::Coreutils qw{ gnu_mv };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Variantcalling::Bcftools qw{ bcftools_norm };
    use MIP::Program::Variantcalling::Gatk
      qw{ gatk_variantrecalibrator gatk_applyvqsr gatk_selectvariants gatk_calculategenotypeposteriors };
    use MIP::QC::Record
      qw{ add_recipe_outfile_to_sample_info add_processing_metafile_to_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Constants
    Readonly my $MAX_GAUSSIAN_LEVEL_INDEL             => 4;
    Readonly my $MAX_GAUSSIAN_LEVEL_SNV               => 6;
    Readonly my $MAX_GAUSSIAN_LEVEL_SNV_SINGLE_SAMPLE => 4;

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
    my $infile_path        = $io{in}{file_path};

    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};
    my $enable_indel_max_gaussians_filter =
      $active_parameter_href->{gatk_variantrecalibration_indel_max_gaussians};
    my $enable_snv_max_gaussians_filter =
      $active_parameter_href->{gatk_variantrecalibration_snv_max_gaussians};
    my $job_id_chain = get_recipe_attributes(
        {
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            attribute      => q{chain},
        }
    );
    my $recipe_mode        = $active_parameter_href->{$recipe_name};
    my $referencefile_path = $active_parameter_href->{human_genome_reference};
    my $resource_indel_href =
      $active_parameter_href->{gatk_variantrecalibration_resource_indel};
    my $resource_snv_href =
      $active_parameter_href->{gatk_variantrecalibration_resource_snv};
    my ( $core_number, $time, @source_environment_cmds ) = get_recipe_parameters(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id               => $job_id_chain,
                id                     => $case_id,
                file_info_href         => $file_info_href,
                file_name_prefixes_ref => [$infile_name_prefix],
                outdata_dir            => $active_parameter_href->{outdata_dir},
                parameter_href         => $parameter_href,
                recipe_name            => $recipe_name,
                temp_directory         => $temp_directory,
            }
        )
    );
    my $outdir_path_prefix       = $io{out}{dir_path_prefix};
    my $outfile_path_prefix      = $io{out}{file_path_prefix};
    my $outfile_suffix           = $io{out}{file_suffix};
    my $outfile_path             = $io{out}{file_path};
    my $temp_outfile_path_prefix = $io{temp}{file_path_prefix};

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

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

    ## Split to enable submission to &sample_info_qc later
    my ( $volume, $directory, $stderr_file ) =
      splitpath( $recipe_info_path . $DOT . q{stderr.txt} );

    ## Create .fam file to be used in variant calling analyses
    my $fam_file_path = catfile( $outdir_path_prefix, $case_id . $DOT . q{fam} );
    create_fam_file(
        {
            active_parameter_href => $active_parameter_href,
            execution_mode        => q{system},
            fam_file_path         => $fam_file_path,
            FILEHANDLE            => $FILEHANDLE,
            parameter_href        => $parameter_href,
            sample_info_href      => $sample_info_href,
        }
    );

    ## Check if "--pedigree" should be included in analysis
    my %commands = gatk_pedigree_flag(
        {
            fam_file_path => $fam_file_path,
            recipe_name   => $recipe_name,
        }
    );

    ### GATK VariantRecalibrator
    ## Set mode to be used in variant recalibration
    # SNP and INDEL will be recalibrated successively in the same file
    # because when you specify eg SNP mode, the indels are emitted
    # without modification, and vice-versa.
    my @modes = qw{ SNP INDEL };

  MODE:
    foreach my $mode (@modes) {

        say {$FILEHANDLE} q{## GATK VariantRecalibrator};

        ## Get parameters
        my $max_gaussian_level;
        my $varrecal_infile_path;

        if ( $mode eq q{SNP} ) {

            $varrecal_infile_path = $infile_path;

            ## Use hard filtering
            if ($enable_snv_max_gaussians_filter) {

                ## Use fewer Gaussians for single sample cases
                if ( scalar @{ $active_parameter_href->{sample_ids} } == 1 ) {
                    $max_gaussian_level = $MAX_GAUSSIAN_LEVEL_SNV_SINGLE_SAMPLE;
                }
                else {
                    $max_gaussian_level = $MAX_GAUSSIAN_LEVEL_SNV;
                }
            }
        }

        ## Use created recalibrated snp vcf as input
        if ( $mode eq q{INDEL} ) {

            $varrecal_infile_path =
              $temp_outfile_path_prefix . $DOT . q{SNV} . $outfile_suffix;

            ## Use hard filtering
            if ($enable_indel_max_gaussians_filter) {

                $max_gaussian_level = $MAX_GAUSSIAN_LEVEL_INDEL;
            }
        }

        my @annotations =
          @{ $active_parameter_href->{gatk_variantrecalibration_annotations} };

        ## Special case: Not to be used with hybrid capture. NOTE: Disabled when analysing wes + wgs in the same run
        if ( $consensus_analysis_type ne q{wgs} ) {

            ## Removes an element from array and return new array while leaving orginal elements_ref untouched
            @annotations = delete_contig_elements(
                {
                    elements_ref       => \@annotations,
                    remove_contigs_ref => [qw{ DP }],
                }
            );
        }

        my @resources;
        if ( $mode eq q{SNP} ) {

            @resources =
              _build_gatk_resource_command( { resources_href => $resource_snv_href, } );
        }
        if ( $mode eq q{INDEL} ) {

            @resources =
              _build_gatk_resource_command( { resources_href => $resource_indel_href, } );
        }

        my $recal_file_path = $temp_outfile_path_prefix . $DOT . q{intervals};
        gatk_variantrecalibrator(
            {
                annotations_ref       => \@annotations,
                FILEHANDLE            => $FILEHANDLE,
                infile_path           => $varrecal_infile_path,
                java_use_large_pages  => $active_parameter_href->{java_use_large_pages},
                max_gaussian_level    => $max_gaussian_level,
                memory_allocation     => q{Xmx24g},
                mode                  => $mode,
                outfile_path          => $recal_file_path,
                referencefile_path    => $referencefile_path,
                resources_ref         => \@resources,
                rscript_file_path     => $recal_file_path . $DOT . q{plots.R},
                temp_directory        => $temp_directory,
                tranches_file_path    => $recal_file_path . $DOT . q{tranches},
                trust_all_polymorphic => $active_parameter_href
                  ->{gatk_variantrecalibration_trust_all_polymorphic},
                verbosity => $active_parameter_href->{gatk_logging_level},
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        ## GATK ApplyVQSR
        say {$FILEHANDLE} q{## GATK ApplyVQSR};

        ## Get parameters
        my $applyvqsr_infile_path;
        my $applyvqsr_outfile_path;
        my $ts_filter_level;

        if ( $mode eq q{SNP} ) {

            $applyvqsr_infile_path = $varrecal_infile_path;
            $applyvqsr_outfile_path =
              $temp_outfile_path_prefix . $DOT . q{SNV} . $outfile_suffix;
            $ts_filter_level =
              $active_parameter_href->{gatk_variantrecalibration_snv_tsfilter_level};
        }

        ## Use created recalibrated snp vcf as input
        if ( $mode eq q{INDEL} ) {

            $applyvqsr_infile_path =
              $temp_outfile_path_prefix . $DOT . q{SNV} . $outfile_suffix;
            $applyvqsr_outfile_path = $outfile_path;
            $ts_filter_level =
              $active_parameter_href->{gatk_variantrecalibration_indel_tsfilter_level};
        }

        gatk_applyvqsr(
            {
                FILEHANDLE           => $FILEHANDLE,
                infile_path          => $applyvqsr_infile_path,
                java_use_large_pages => $active_parameter_href->{java_use_large_pages},
                memory_allocation    => q{Xmx10g},
                mode                 => $mode,
                outfile_path         => $applyvqsr_outfile_path,
                recal_file_path      => $recal_file_path,
                referencefile_path   => $referencefile_path,
                temp_directory       => $temp_directory,
                tranches_file_path   => $recal_file_path . $DOT . q{tranches},
                ts_filter_level      => $ts_filter_level,
                verbosity            => $active_parameter_href->{gatk_logging_level},
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    ## GenotypeRefinement
    if ( $parameter_href->{dynamic_parameter}{trio} ) {

        say {$FILEHANDLE} q{## GATK CalculateGenotypePosteriors};

        my $calculategt_outfile_path =
          $outfile_path_prefix . $UNDERSCORE . q{refined} . $outfile_suffix;
        gatk_calculategenotypeposteriors(
            {
                FILEHANDLE           => $FILEHANDLE,
                infile_path          => $outfile_path,
                java_use_large_pages => $active_parameter_href->{java_use_large_pages},
                memory_allocation    => q{Xmx6g},
                outfile_path         => $calculategt_outfile_path,
                pedigree             => $commands{pedigree},
                referencefile_path   => $referencefile_path,
                supporting_callset_file_path =>
                  $active_parameter_href->{gatk_calculategenotypeposteriors_support_set},
                temp_directory => $temp_directory,
                verbosity      => $active_parameter_href->{gatk_logging_level},
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        ## Change name of file to accomodate downstream
        gnu_mv(
            {
                FILEHANDLE   => $FILEHANDLE,
                infile_path  => $calculategt_outfile_path,
                outfile_path => $outfile_path,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    ## BcfTools norm, Left-align and normalize indels, split multiallelics
    my $bcftools_outfile_path =
      $outfile_path_prefix . $UNDERSCORE . q{normalized} . $outfile_suffix;
    bcftools_norm(
        {
            FILEHANDLE     => $FILEHANDLE,
            infile_path    => $outfile_path,
            multiallelic   => q{-},
            output_type    => q{v},
            outfile_path   => $bcftools_outfile_path,
            reference_path => $referencefile_path,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Change name of file to accomodate downstream
    gnu_mv(
        {
            FILEHANDLE   => $FILEHANDLE,
            infile_path  => $bcftools_outfile_path,
            outfile_path => $outfile_path,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    close $FILEHANDLE;

    if ( $recipe_mode == 1 ) {

        ## Collect QC metadata info for later use
        add_recipe_outfile_to_sample_info(
            {
                path             => $outfile_path,
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
            }
        );

        # Used to find order of samples in qccollect downstream
        add_recipe_outfile_to_sample_info(
            {
                path             => $outfile_path,
                recipe_name      => q{pedigree_check},
                sample_info_href => $sample_info_href,
            }
        );

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

sub _build_gatk_resource_command {

## Function : Build resources in the correct format for GATK
## Returns  :
## Arguments: $resources_href => Resources to build comand for {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $resources_href;

    my $tmpl = {
        resources_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$resources_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @built_resources;

  RESOURCE:
    while ( my ( $file, $string ) = each %{$resources_href} ) {

        ## Build resource string
        push @built_resources, $string . $COLON . $file;
    }
    return @built_resources;
}

1;
