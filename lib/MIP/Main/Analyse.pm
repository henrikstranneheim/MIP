package MIP::Main::Analyse;

#### Master script for analysing paired end reads from the Illumina plattform in fastq(.gz) format to annotated ranked disease causing variants. The program performs QC, aligns reads using BWA, performs variant discovery and annotation as well as ranking the found variants according to disease potential.

#### Copyright 2011 Henrik Stranneheim

use 5.018;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use Cwd qw{ abs_path };
use English qw{ -no_match_vars };
use File::Basename qw{ basename dirname fileparse };
use File::Copy qw{ copy };
use File::Spec::Functions qw{ catdir catfile devnull splitpath };
use FindBin qw{ $Bin };
use Getopt::Long;
use IPC::Cmd qw{ can_run run};
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use POSIX;
use Time::Piece;
use utf8;
use warnings qw{ FATAL utf8 };

## Third party module(s)
use autodie qw{ open close :all };
use IPC::System::Simple;
use List::MoreUtils qw { any uniq all };
use Modern::Perl qw{ 2014 };
use Path::Iterator::Rule;
use Readonly;

## MIPs lib/
# Add MIPs internal lib
use MIP::Check::Cluster qw{ check_max_core_number };
use MIP::Check::Modules qw{ check_perl_modules };
use MIP::Check::Parameter qw{ check_allowed_temp_directory
  check_aligner
  check_cmd_config_vs_definition_file
  check_email_address
  check_parameter_hash
  check_pprogram_exists_in_hash
  check_program_mode
  check_sample_ids
  check_sample_id_in_hash_parameter
  check_sample_id_in_hash_parameter_path
  check_snpsift_keys
  check_vep_directories
};
use MIP::Check::Path
  qw{ check_command_in_path check_parameter_files check_target_bed_file_suffix check_vcfanno_toml };
use MIP::Check::Reference
  qw{ check_human_genome_file_endings check_parameter_metafiles };
use MIP::File::Format::Config qw{ write_mip_config };
use MIP::File::Format::Pedigree
  qw{ create_fam_file detect_founders detect_sample_id_gender detect_trio parse_yaml_pedigree_file reload_previous_pedigree_info };
use MIP::File::Format::Yaml qw{ load_yaml write_yaml order_parameter_names };
use MIP::Get::Analysis qw{ get_overall_analysis_type };
use MIP::Get::File qw{ get_select_file_contigs };
use MIP::Log::MIP_log4perl qw{ initiate_logger set_default_log4perl_file };
use MIP::Parse::Parameter
  qw{ parse_prioritize_variant_callers parse_start_with_program };
use MIP::Script::Utils qw{ help };
use MIP::Set::Contigs qw{ set_contigs };
use MIP::Set::Parameter
  qw{ set_config_to_active_parameters set_custom_default_to_active_parameter set_default_config_dynamic_parameters set_default_to_active_parameter set_dynamic_parameter set_human_genome_reference_features set_parameter_reference_dir_path set_parameter_to_broadcast };
use MIP::Update::Contigs qw{ update_contigs_for_run };
use MIP::Update::Parameters
  qw{ update_dynamic_config_parameters update_exome_target_bed update_reference_parameters update_vcfparser_outfile_counter };
use MIP::Update::Path qw{ update_to_absolute_path };
use MIP::Update::Programs
  qw{ update_prioritize_flag update_program_mode_for_analysis_type update_program_mode_with_dry_run_all };

## Recipes
use MIP::Recipes::Analysis::Gzip_fastq qw{ analysis_gzip_fastq };
use MIP::Recipes::Analysis::Split_fastq_file qw{ analysis_split_fastq_file };
use MIP::Recipes::Analysis::Vt_core qw{ analysis_vt_core };
use MIP::Recipes::Pipeline::Rare_disease qw{ pipeline_rare_disease };
use MIP::Recipes::Pipeline::Rna qw{ pipeline_rna };
use MIP::Recipes::Pipeline::Cancer qw{ pipeline_cancer };

BEGIN {

    use base qw{ Exporter };
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.06;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ mip_analyse };
}

## Constants
Readonly my $DOT          => q{.};
Readonly my $EMPTY_STR    => q{};
Readonly my $NEWLINE      => qq{\n};
Readonly my $SINGLE_QUOTE => q{'};
Readonly my $SPACE        => q{ };
Readonly my $TAB          => qq{\t};

sub mip_analyse {

## Function : Creates program directories (info & programData & programScript), program script filenames and writes sbatch header.
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $file_info_href                       => File info hash {REF}
##          : $parameter_href        => Parameter hash {REF}
#           : $order_parameters_ref  => Order of addition to parameter array {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $parameter_href;
    my $file_info_href;
    my $order_parameters_ref;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        order_parameters_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$order_parameters_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Transfer to lexical variables
    my %active_parameter = %{$active_parameter_href};
    my %file_info        = %{$file_info_href};
    my @order_parameters = @{$order_parameters_ref};
    my %parameter        = %{$parameter_href};

#### Script parameters

## Add date_time_stamp for later use in log and qc_metrics yaml file
    my $date_time       = localtime;
    my $date_time_stamp = $date_time->datetime;
    my $date            = $date_time->ymd;

    # Catches script name and removes ending
    my $script = fileparse( basename( $PROGRAM_NAME, $DOT . q{pl} ) );
    chomp( $date_time_stamp, $date, $script );

#### Set program parameters

## Set MIP version
    our $VERSION = 'v7.0.1';

    if ( $active_parameter{version} ) {

        say STDOUT $NEWLINE . basename($PROGRAM_NAME) . $SPACE . $VERSION,
          $NEWLINE;
        exit;
    }

## Directories, files, job_ids and sample_info
    my ( %infile, %indir_path, %infile_lane_prefix, %lane,
        %infile_both_strands_prefix, %job_id, %sample_info );

#### Staging Area
### Get and/or set input parameters

## Special case for boolean flag that will be removed from
## config upon loading
    my @boolean_parameter = qw{dry_run_all};
    foreach my $parameter (@boolean_parameter) {

        if ( not defined $active_parameter{$parameter} ) {

            delete $active_parameter{$parameter};
        }
    }

## Change relative path to absolute path for parameter with "update_path: absolute_path" in config
    update_to_absolute_path(
        {
            active_parameter_href => \%active_parameter,
            parameter_href        => \%parameter,
        }
    );

### Config file
## If config from cmd
    if ( exists $active_parameter{config_file}
        && defined $active_parameter{config_file} )
    {

        ## Loads a YAML file into an arbitrary hash and returns it.
        my %config_parameter =
          load_yaml( { yaml_file => $active_parameter{config_file}, } );

        ## Remove previous analysis specific info not relevant for current run e.g. log file, which is read from pedigree or cmd
        my @remove_keys = (qw{ log_file dry_run_all });

      KEY:
        foreach my $key (@remove_keys) {

            delete $config_parameter{$key};
        }

## Set config parameters into %active_parameter unless $parameter
## has been supplied on the command line
        set_config_to_active_parameters(
            {
                active_parameter_href => \%active_parameter,
                config_parameter_href => \%config_parameter,
            }
        );

        ## Compare keys from config and cmd (%active_parameter) with definitions file (%parameter)
        check_cmd_config_vs_definition_file(
            {
                active_parameter_href => \%active_parameter,
                parameter_href        => \%parameter,
            }
        );

        my @config_dynamic_parameters =
          qw{ analysis_constant_path outaligner_dir };

        ## Replace config parameter with cmd info for config dynamic parameter
        set_default_config_dynamic_parameters(
            {
                active_parameter_href => \%active_parameter,
                parameter_href        => \%parameter,
                parameter_names_ref   => \@config_dynamic_parameters,
            }
        );

        ## Loop through all parameters and update info
      PARAMETER:
        foreach my $parameter_name (@order_parameters) {

            ## Updates the active parameters to particular user/cluster for dynamic config parameters following specifications. Leaves other entries untouched.
            update_dynamic_config_parameters(
                {
                    active_parameter_href => \%active_parameter,
                    parameter_name        => $parameter_name,
                }
            );
        }
    }

## Set the default Log4perl file using supplied dynamic parameters.
    $active_parameter{log_file} = set_default_log4perl_file(
        {
            active_parameter_href => \%active_parameter,
            cmd_input             => $active_parameter{log_file},
            date                  => $date,
            date_time_stamp       => $date_time_stamp,
            script                => $script,
        }
    );

## Creates log object
    my $log = initiate_logger(
        {
            file_path => $active_parameter{log_file},
            log_name  => q{MIP},
        }
    );

## Parse pedigree file
## Reads family_id_pedigree file in YAML format. Checks for pedigree data for allowed entries and correct format. Add data to sample_info depending on user info.
    # Meta data in YAML format
    if ( defined $active_parameter{pedigree_file} ) {

        ## Loads a YAML file into an arbitrary hash and returns it. Load parameters from previous run from sample_info_file
        my %pedigree =
          load_yaml( { yaml_file => $active_parameter{pedigree_file}, } );

        $log->info( q{Loaded: } . $active_parameter{pedigree_file} );

        parse_yaml_pedigree_file(
            {
                active_parameter_href => \%active_parameter,
                file_path             => $active_parameter{pedigree_file},
                parameter_href        => \%parameter,
                pedigree_href         => \%pedigree,
                sample_info_href      => \%sample_info,
            }
        );
    }

# Detect if all samples has the same sequencing type and return consensus if reached
    $parameter{dynamic_parameter}{consensus_analysis_type} =
      get_overall_analysis_type(
        { analysis_type_href => \%{ $active_parameter{analysis_type} }, } );

### Populate uninitilized active_parameters{parameter_name} with default from parameter
  PARAMETER:
    foreach my $parameter_name (@order_parameters) {

        ## If hash and set - skip
        next PARAMETER
          if ( ref $active_parameter{$parameter_name} eq qw{HASH}
            && keys %{ $active_parameter{$parameter_name} } );

        ## If array and set - skip
        next PARAMETER
          if ( ref $active_parameter{$parameter_name} eq qw{ARRAY}
            && @{ $active_parameter{$parameter_name} } );

        ## If scalar and set - skip
        next PARAMETER
          if ( defined $active_parameter{$parameter_name}
            and not ref $active_parameter{$parameter_name} );

        ### Special case for parameters that are dependent on other parameters values
        my @custom_default_parameters = qw{
          analysis_type
          bwa_build_reference
          exome_target_bed
          expansionhunter_repeat_specs_dir
          gatk_path
          infile_dirs
          picardtools_path
          sample_info_file
          snpeff_path
          rtg_vcfeval_reference_genome
          vep_directory_path
        };

        if ( any { $_ eq $parameter_name } @custom_default_parameters ) {

            set_custom_default_to_active_parameter(
                {
                    active_parameter_href => \%active_parameter,
                    parameter_href        => \%parameter,
                    parameter_name        => $parameter_name,
                }
            );
            next PARAMETER;
        }

        ## Checks and sets user input or default values to active_parameters
        set_default_to_active_parameter(
            {
                active_parameter_href => \%active_parameter,
                associated_programs_ref =>
                  \@{ $parameter{$parameter_name}{associated_program} },
                log            => $log,
                parameter_href => \%parameter,
                parameter_name => $parameter_name,
            }
        );
    }

## Update path for supplied reference(s) associated with parameter that should reside in the mip reference directory to full path
    set_parameter_reference_dir_path(
        {
            active_parameter_href => \%active_parameter,
            parameter_name        => q{human_genome_reference},
        }
    );

## Detect version and source of the human_genome_reference: Source (hg19 or GRCh).
    set_human_genome_reference_features(
        {
            file_info_href => \%file_info,
            human_genome_reference =>
              basename( $active_parameter{human_genome_reference} ),
            log => $log,
        }
    );

## Update exome_target_bed files with human_genome_reference_source and human_genome_reference_version
    update_exome_target_bed(
        {
            exome_target_bed_file_href => $active_parameter{exome_target_bed},
            human_genome_reference_source =>
              $file_info{human_genome_reference_source},
            human_genome_reference_version =>
              $file_info{human_genome_reference_version},
        }
    );

    # Holds all active parameters values for broadcasting
    my @broadcasts;

    if ( $active_parameter{verbose} ) {

        set_parameter_to_broadcast(
            {
                parameter_href        => \%parameter,
                active_parameter_href => \%active_parameter,
                order_parameters_ref  => \@order_parameters,
                broadcasts_ref        => \@broadcasts,
            }
        );
    }

## Reference in MIP reference directory
  PARAMETER:
    foreach my $parameter_name ( keys %parameter ) {

        ## Expect file to be in reference directory
        if ( exists $parameter{$parameter_name}{reference} ) {

            update_reference_parameters(
                {
                    active_parameter_href => \%active_parameter,
                    associated_programs_ref =>
                      \@{ $parameter{$parameter_name}{associated_program} },
                    parameter_name => $parameter_name,
                }
            );
        }
    }

### Checks

## Check existence of files and directories
  PARAMETER:
    foreach my $parameter_name ( keys %parameter ) {

        if ( exists $parameter{$parameter_name}{exists_check} ) {

            check_parameter_files(
                {
                    active_parameter_href => \%active_parameter,
                    associated_programs_ref =>
                      \@{ $parameter{$parameter_name}{associated_program} },
                    log => $log,
                    parameter_exists_check =>
                      $parameter{$parameter_name}{exists_check},
                    parameter_href => \%parameter,
                    parameter_name => $parameter_name,
                }
            );
        }
    }

## Updates sample_info hash with previous run pedigree info
    reload_previous_pedigree_info(
        {
            log                   => $log,
            sample_info_href      => \%sample_info,
            sample_info_file_path => $active_parameter{sample_info_file},
        }
    );

## Special case since dict is created with .fastq removed
## Check the existance of associated human genome files
    check_human_genome_file_endings(
        {
            active_parameter_href => \%active_parameter,
            file_info_href        => \%file_info,
            log                   => $log,
            parameter_href        => \%parameter,
            parameter_name        => q{human_genome_reference},
        }
    );

## Check that supplied target file ends with ".bed" and otherwise croaks
  TARGET_FILE:
    foreach
      my $target_bed_file ( keys %{ $active_parameter{exome_target_bed} } )
    {

        check_target_bed_file_suffix(
            {
                parameter_name => q{exome_target_bed},
                path           => $target_bed_file,
            }
        );
    }

## Checks parameter metafile exists and set build_file parameter
    check_parameter_metafiles(
        {
            parameter_href        => \%parameter,
            active_parameter_href => \%active_parameter,
            file_info_href        => \%file_info,
        }
    );

## Update the expected number of outfile after vcfparser
    update_vcfparser_outfile_counter(
        { active_parameter_href => \%active_parameter, } );

## Collect select file contigs to loop over downstream
    if ( $active_parameter{vcfparser_select_file} ) {

## Collects sequences contigs used in select file
        @{ $file_info{select_file_contigs} } = get_select_file_contigs(
            {
                select_file_path =>
                  catfile( $active_parameter{vcfparser_select_file} ),
                log => $log,
            }
        );
    }

## Detect family constellation based on pedigree file
    $parameter{dynamic_parameter}{trio} = detect_trio(
        {
            active_parameter_href => \%active_parameter,
            log                   => $log,
            sample_info_href      => \%sample_info,
        }
    );

## Detect number of founders (i.e. parents ) based on pedigree file
    detect_founders(
        {
            active_parameter_href => \%active_parameter,
            sample_info_href      => \%sample_info,
        }
    );

## Check email adress syntax and mail host
    if ( defined $active_parameter{email} ) {

        check_email_address(
            {
                email => $active_parameter{email},
                log   => $log,
            }
        );
    }

## Check that the temp directory value is allowed
    check_allowed_temp_directory(
        {
            log            => $log,
            temp_directory => $active_parameter{temp_directory},
        }
    );

## Parameters that have keys as MIP program names
    my @parameter_keys_to_check =
      (qw{ module_time module_core_number module_source_environment_command });
  PARAMETER_NAME:
    foreach my $parameter_name (@parameter_keys_to_check) {

        ## Test if key from query hash exists truth hash
        check_pprogram_exists_in_hash(
            {
                log            => $log,
                parameter_name => $parameter_name,
                query_ref      => \%{ $active_parameter{$parameter_name} },
                truth_href     => \%parameter,
            }
        );
    }

## Parameters with key(s) that have elements as MIP program names
    my @parameter_element_to_check = qw(associated_program);
  PARAMETER:
    foreach my $parameter ( keys %parameter ) {

      KEY:
        foreach my $parameter_name (@parameter_element_to_check) {

            next KEY if ( not exists $parameter{$parameter}{$parameter_name} );

            ## Test if element from query array exists truth hash
            check_pprogram_exists_in_hash(
                {
                    log            => $log,
                    parameter_name => $parameter_name,
                    query_ref  => \@{ $parameter{$parameter}{$parameter_name} },
                    truth_href => \%parameter,
                }
            );
        }
    }

## Parameters that have elements as MIP program names
    my @parameter_elements_to_check =
      (qw(associated_program decompose_normalize_references));
    foreach my $parameter_name (@parameter_elements_to_check) {

        ## Test if element from query array exists truth hash
        check_pprogram_exists_in_hash(
            {
                log            => $log,
                parameter_name => $parameter_name,
                query_ref      => \@{ $active_parameter{$parameter_name} },
                truth_href     => \%parameter,
            }
        );
    }

## Check that the module core number do not exceed the maximum per node
    foreach my $program_name ( keys %{ $active_parameter{module_core_number} } )
    {

        ## Limit number of cores requested to the maximum number of cores available per node
        $active_parameter{module_core_number}{$program_name} =
          check_max_core_number(
            {
                max_cores_per_node => $active_parameter{max_cores_per_node},
                core_number_requested =>
                  $active_parameter{module_core_number}{$program_name},
            }
          );
    }

## Check programs in path, and executable
    check_command_in_path(
        {
            active_parameter_href => \%active_parameter,
            log                   => $log,
            parameter_href        => \%parameter,
        }
    );

## Test that the family_id and the sample_id(s) exists and are unique. Check if id sample_id contains "_".
    check_sample_ids(
        {
            family_id      => $active_parameter{family_id},
            log            => $log,
            sample_ids_ref => \@{ $active_parameter{sample_ids} },
        }
    );

## Check sample_id provided in hash parameter is included in the analysis and only represented once
    check_sample_id_in_hash_parameter(
        {
            active_parameter_href => \%active_parameter,
            log                   => $log,
            parameter_names_ref =>
              [qw{ analysis_type expected_coverage sample_origin }],
            parameter_href => \%parameter,
            sample_ids_ref => \@{ $active_parameter{sample_ids} },
        }
    );

## Check sample_id provided in hash path parameter is included in the analysis and only represented once
    check_sample_id_in_hash_parameter_path(
        {
            active_parameter_href => \%active_parameter,
            log                   => $log,
            parameter_names_ref   => [qw{ infile_dirs exome_target_bed }],
            sample_ids_ref        => \@{ $active_parameter{sample_ids} },
        }
    );

## Check that VEP directory and VEP cache match
    if ( exists $active_parameter{pvarianteffectpredictor} ) {
        check_vep_directories(
            {
                log                 => $log,
                vep_directory_cache => $active_parameter{vep_directory_cache},
                vep_directory_path  => $active_parameter{vep_directory_path},
            }
        );
    }

## Check that the supplied vcfanno toml frequency file match record 'file=' within toml config file
    if (    exists $active_parameter{psv_combinevariantcallsets}
        and $active_parameter{psv_combinevariantcallsets} > 0
        and $active_parameter{sv_vcfanno} > 0 )
    {

        check_vcfanno_toml(
            {
                log               => $log,
                vcfanno_file_freq => $active_parameter{sv_vcfanno_config_file},
                vcfanno_file_toml => $active_parameter{sv_vcfanno_config},
            }
        );
    }

    check_snpsift_keys(
        {
            log => $log,
            snpsift_annotation_files_href =>
              \%{ $active_parameter{snpsift_annotation_files} },
            snpsift_annotation_outinfo_key_href =>
              \%{ $active_parameter{snpsift_annotation_outinfo_key} },
        }
    );

## Adds dynamic aggregate information from definitions to parameter hash
    set_dynamic_parameter(
        {
            aggregates_ref => [
                ## Collects all programs that MIP can handle
                q{type:program},
                ## Collects all variant_callers
                q{program_type:variant_callers},
                ## Collects all structural variant_callers
                q{program_type:structural_variant_callers},
                ## Collect all aligners
                q{program_type:aligners},
                ## Collects all references in that are supposed to be in reference directory
                q{reference:reference_dir},
            ],
            parameter_href => \%parameter,
        }
    );

## Check correct value for program mode in MIP
    check_program_mode(
        {
            active_parameter_href => \%active_parameter,
            log                   => $log,
            parameter_href        => \%parameter,
        }
    );

## Get initiation program, downstream dependencies and update program modes for start_with_program parameter
    parse_start_with_program(
        {
            active_parameter_href => \%active_parameter,
            initiation_file =>
              catfile( $Bin, qw{ definitions rare_disease_initiation.yaml } ),
            parameter_href => \%parameter,
        },
    );

## Update program mode depending on dry_run_all flag
    update_program_mode_with_dry_run_all(
        {
            active_parameter_href => \%active_parameter,
            dry_run_all           => $active_parameter{dry_run_all},
            programs_ref => \@{ $parameter{dynamic_parameter}{program} },
        }
    );

## Check that the correct number of aligners is used in MIP and sets the aligner flag accordingly
    check_aligner(
        {
            active_parameter_href => \%active_parameter,
            broadcasts_ref        => \@broadcasts,
            log                   => $log,
            parameter_href        => \%parameter,
            verbose               => $active_parameter{verbose},
        }
    );

## Check that all active variant callers have a prioritization order and that the prioritization elements match a supported variant caller
    parse_prioritize_variant_callers(
        {
            active_parameter_href => \%active_parameter,
            log                   => $log,
            parameter_href        => \%parameter,
        }
    );

## Broadcast set parameters info
    foreach my $parameter_info (@broadcasts) {

        $log->info($parameter_info);
    }

## Update program mode depending on analysis run value as some programs are not applicable for e.g. wes
    update_program_mode_for_analysis_type(
        {
            active_parameter_href => \%active_parameter,
            consensus_analysis_type =>
              $parameter{dynamic_parameter}{consensus_analysis_type},
            log          => $log,
            programs_ref => [
                qw{ cnvnator delly_call delly_reformat expansionhunter tiddit samtools_subsample_mt }
            ],
        }
    );

## Update prioritize flag depending on analysis run value as some programs are not applicable for e.g. wes
    $active_parameter{sv_svdb_merge_prioritize} = update_prioritize_flag(
        {
            consensus_analysis_type =>
              $parameter{dynamic_parameter}{consensus_analysis_type},
            prioritize_key => $active_parameter{sv_svdb_merge_prioritize},
            programs_ref   => [qw{ cnvnator delly_call delly_reformat tiddit }],
        }
    );

## Write config file for family
    write_mip_config(
        {
            active_parameter_href => \%active_parameter,
            log                   => $log,
            remove_keys_ref       => [qw{ associated_program }],
            sample_info_href      => \%sample_info,
        }
    );

## Detect the gender(s) included in current analysis
    (

        $active_parameter{found_male},
        $active_parameter{found_female},
        $active_parameter{found_other},
        $active_parameter{found_other_count},
      )
      = detect_sample_id_gender(
        {
            active_parameter_href => \%active_parameter,
            sample_info_href      => \%sample_info,
        }
      );

### Contigs
## Set contig prefix and contig names depending on reference used
    set_contigs(
        {
            file_info_href         => \%file_info,
            human_genome_reference => $active_parameter{human_genome_reference},
        }
    );

## Update contigs depending on settings in run (wes or if only male samples)
    update_contigs_for_run(
        {
            file_info_href     => \%file_info,
            analysis_type_href => \%{ $active_parameter{analysis_type} },
            found_male         => $active_parameter{found_male},
        }
    );

## Sorts array depending on reference array. NOTE: Only entries present in reference array will survive in sorted array.
    @{ $file_info{sorted_select_file_contigs} } = size_sort_select_file_contigs(
        {
            file_info_href => \%file_info,
            consensus_analysis_type_ref =>
              \$parameter{dynamic_parameter}{consensus_analysis_type},
            hash_key_to_sort        => 'select_file_contigs',
            hash_key_sort_reference => 'contigs_size_ordered',
        }
    );

    if ( $active_parameter{verbose} ) {
## Write CMD to MIP log file
        write_cmd_mip_log(
            {
                parameter_href        => \%parameter,
                active_parameter_href => \%active_parameter,
                order_parameters_ref  => \@order_parameters,
                script_ref            => \$script,
                log_file_ref          => \$active_parameter{log_file},
                mip_version_ref       => \$VERSION,
            }
        );
    }

## Collects the ".fastq(.gz)" files from the supplied infiles directory. Checks if any of the files exist
    collect_infiles(
        {
            active_parameter_href => \%active_parameter,
            indir_path_href       => \%indir_path,
            infile_href           => \%infile,
        }
    );

## Reformat files for MIP output, which have not yet been created into, correct format so that a sbatch script can be generated with the correct filenames
    my $uncompressed_file_switch = infiles_reformat(
        {
            active_parameter_href           => \%active_parameter,
            sample_info_href                => \%sample_info,
            file_info_href                  => \%file_info,
            infile_href                     => \%infile,
            indir_path_href                 => \%indir_path,
            infile_lane_prefix_href         => \%infile_lane_prefix,
            infile_both_strands_prefix_href => \%infile_both_strands_prefix,
            lane_href                       => \%lane,
            job_id_href                     => \%job_id,
            outaligner_dir_ref => \$active_parameter{outaligner_dir},
            program_name       => 'infiles_reformat',
        }
    );

## Creates all fileendings as the samples is processed depending on the chain of modules activated
    create_file_endings(
        {
            parameter_href          => \%parameter,
            active_parameter_href   => \%active_parameter,
            file_info_href          => \%file_info,
            infile_lane_prefix_href => \%infile_lane_prefix,
            order_parameters_ref    => \@order_parameters,
        }
    );

## Create .fam file to be used in variant calling analyses
    create_fam_file(
        {
            parameter_href        => \%parameter,
            active_parameter_href => \%active_parameter,
            sample_info_href      => \%sample_info,
            execution_mode        => 'system',
            fam_file_path         => catfile(
                $active_parameter{outdata_dir},
                $active_parameter{family_id},
                $active_parameter{family_id} . '.fam'
            ),
        }
    );

## Add to SampleInfo
    add_to_sample_info(
        {
            active_parameter_href => \%active_parameter,
            sample_info_href      => \%sample_info,
            file_info_href        => \%file_info,
        }
    );

############
####MAIN####
############

    if ( not $active_parameter{dry_run_all} ) {

        my %no_dry_run_info = (
            analysisrunstatus => q{not_finished},
            analysis_date     => $date_time_stamp,
            mip_version       => $VERSION,
        );

      KEY_VALUE_PAIR:
        while ( my ( $key, $value ) = each %no_dry_run_info ) {

            $sample_info{$key} = $value;
        }
    }

    my $consensus_analysis_type =
      $parameter{dynamic_parameter}{consensus_analysis_type};

## Split of fastq files in batches
    if ( $active_parameter{psplit_fastq_file} ) {

        $log->info(q{[Split fastq files in batches]});

      SAMPLE_ID:
        foreach my $sample_id ( @{ $active_parameter{sample_ids} } ) {

            ## Split input fastq files into batches of reads, versions and compress. Moves original file to subdirectory
            analysis_split_fastq_file(
                {
                    parameter_href        => \%parameter,
                    active_parameter_href => \%active_parameter,
                    infile_href           => \%infile,
                    job_id_href           => \%job_id,
                    insample_directory    => $indir_path{$sample_id},
                    outsample_directory   => $indir_path{$sample_id},
                    sample_id             => $sample_id,
                    program_name          => q{split_fastq_file},
                    sequence_read_batch =>
                      $active_parameter{split_fastq_file_read_batch},
                }
            );
        }

        ## End here if this module is turned on
        exit;
    }

## GZip of fastq files
    if (   $active_parameter{pgzip_fastq}
        && $uncompressed_file_switch eq q{uncompressed} )
    {

        $log->info(q{[Gzip for fastq files]});

      SAMPLES:
        foreach my $sample_id ( @{ $active_parameter{sample_ids} } ) {

            ## Determine which sample id had the uncompressed files
          INFILES:
            foreach my $infile ( @{ $infile{$sample_id} } ) {

                my $infile_suffix = $parameter{pgzip_fastq}{infile_suffix};

                if ( $infile =~ /$infile_suffix$/ ) {

                    ## Automatically gzips fastq files
                    analysis_gzip_fastq(
                        {
                            parameter_href          => \%parameter,
                            active_parameter_href   => \%active_parameter,
                            sample_info_href        => \%sample_info,
                            infile_href             => \%infile,
                            infile_lane_prefix_href => \%infile_lane_prefix,
                            job_id_href             => \%job_id,
                            insample_directory      => $indir_path{$sample_id},
                            sample_id               => $sample_id,
                            program_name            => q{gzip_fastq},
                        }
                    );

                    # Call once per sample_id
                    last INFILES;
                }
            }
        }
    }

### Cancer
    if ( $consensus_analysis_type eq q{cancer} )

    {

        $log->info( q{Pipeline analysis type: } . $consensus_analysis_type );

        ## Pipeline recipe for cancer data
        pipeline_cancer(
            {
                parameter_href          => \%parameter,
                active_parameter_href   => \%active_parameter,
                sample_info_href        => \%sample_info,
                file_info_href          => \%file_info,
                indir_path_href         => \%indir_path,
                infile_href             => \%infile,
                infile_lane_prefix_href => \%infile_lane_prefix,
                lane_href               => \%lane,
                job_id_href             => \%job_id,
                outaligner_dir          => $active_parameter{outaligner_dir},
                log                     => $log,
            }
        );
    }

### RNA
    if ( $consensus_analysis_type eq q{wts} ) {

        $log->info( q{Pipeline analysis type: } . $consensus_analysis_type );

        ## Pipeline recipe for rna data
        pipeline_rna(
            {
                parameter_href          => \%parameter,
                active_parameter_href   => \%active_parameter,
                sample_info_href        => \%sample_info,
                file_info_href          => \%file_info,
                indir_path_href         => \%indir_path,
                infile_href             => \%infile,
                infile_lane_prefix_href => \%infile_lane_prefix,
                lane_href               => \%lane,
                job_id_href             => \%job_id,
                outaligner_dir          => $active_parameter{outaligner_dir},
                log                     => $log,
            }
        );
    }

### WES|WGS
    if (   $consensus_analysis_type eq q{wgs}
        || $consensus_analysis_type eq q{wes}
        || $consensus_analysis_type eq q{mixed} )
    {

        $log->info( q{Pipeline analysis type: } . $consensus_analysis_type );

        ## Pipeline recipe for rna data
        pipeline_rare_disease(
            {
                parameter_href          => \%parameter,
                active_parameter_href   => \%active_parameter,
                sample_info_href        => \%sample_info,
                file_info_href          => \%file_info,
                indir_path_href         => \%indir_path,
                infile_href             => \%infile,
                infile_lane_prefix_href => \%infile_lane_prefix,
                lane_href               => \%lane,
                job_id_href             => \%job_id,
                outaligner_dir          => $active_parameter{outaligner_dir},
                log                     => $log,
            }
        );
    }

## Write QC for programs used in analysis
    # Write SampleInfo to yaml file
    if ( $active_parameter{sample_info_file} ) {

        ## Writes a YAML hash to file
        write_yaml(
            {
                yaml_href      => \%sample_info,
                yaml_file_path => $active_parameter{sample_info_file},
            }
        );
        $log->info( q{Wrote: } . $active_parameter{sample_info_file} );
    }

}
######################
####Sub routines######
######################

sub collect_infiles {

##collect_infiles

##Function : Collects the ".fastq(.gz)" files from the supplied infiles directory. Checks if any files exist.
##Returns  : ""
##Arguments: $active_parameter_href, $indir_path_href, $infile_href
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $indir_path_href       => Indirectories path(s) hash {REF}
##         : $infile_href           => Infiles hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $indir_path_href;
    my $infile_href;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        indir_path_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$indir_path_href
        },
        infile_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_href
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    $log->info("Reads from platform:\n");

    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } )
    {    #Collects inputfiles govern by sample_ids

        ## Return the key if the hash value and query match
        my $infile_directory_ref = \get_matching_values_key(
            {
                active_parameter_href => $active_parameter_href,
                query_value_ref       => \$sample_id,
                parameter_name        => "infile_dirs",
            }
        );

        my @infiles;

        ## Collect all fastq files in supplied indirectories
        my $rule = Path::Iterator::Rule->new;
        $rule->skip_subdirs("original_fastq_files")
          ;    #Ignore if original fastq files sub directory
        $rule->name("*.fastq*");    #Only look for fastq or fastq.gz files
        my $it = $rule->iter($$infile_directory_ref);

        while ( my $file = $it->() ) {    #Iterate over directory

            my ( $volume, $directory, $fastq_file ) = splitpath($file);
            push( @infiles, $fastq_file );
        }
        chomp(@infiles);    #Remove newline from every entry in array

        if ( !@infiles ) {  #No "*.fastq*" infiles

            $log->fatal(
"Could not find any '.fastq' files in supplied infiles directory "
                  . $$infile_directory_ref,
                "\n"
            );
            exit 1;
        }
        foreach my $infile (@infiles)
        {    #Check that inFileDirs/infile contains sample_id in filename

            unless ( $infile =~ /$sample_id/ ) {

                $log->fatal(
                    "Could not detect sample_id: "
                      . $sample_id
                      . " in supplied infile: "
                      . $$infile_directory_ref . "/"
                      . $infile,
                    "\n"
                );
                $log->fatal(
"Check that: '--sample_ids' and '--inFileDirs' contain the same sample_id and that the filename of the infile contains the sample_id.",
                    "\n"
                );
                exit 1;
            }
        }
        $log->info( "Sample id: " . $sample_id . "\n" );
        $log->info("\tInputfiles:\n");

        ## Log each file from platform
        foreach my $file (@infiles) {

            $log->info( "\t\t", $file, "\n" );    #Indent for visability
        }
        $indir_path_href->{$sample_id} =
          $$infile_directory_ref;                 #Catch inputdir path
        $infile_href->{$sample_id} = [@infiles];  #Reload files into hash
    }
}

sub infiles_reformat {

## Function : Reformat files for MIP output, which have not yet been created into, correct format so that a sbatch script can be generated with the correct filenames.
## Returns  : "$uncompressed_file_counter"
## Arguments: $active_parameter_href           => Active parameters for this analysis hash {REF}
##          : $file_info_href                  => File info hash {REF}
##          : $indir_path_href                 => Indirectories path(s) hash {REF}
##          : $infile_both_strands_prefix_href => The infile(s) without the ".ending" and strand info {REF}
##          : $infile_href                     => Infiles hash {REF}
##          : $infile_lane_prefix_href         => Infile(s) without the ".ending" {REF}
##          : $job_id_href                     => Job id hash {REF}
##          : $lane_href                       => The lane info hash {REF}
##          : $outaligner_dir_ref              => Outaligner_dir used in the analysis {REF}
##          : $program_name                    => Program name {REF}
##          : $sample_info_href                => Info on samples and family hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $indir_path_href;
    my $infile_both_strands_prefix_href;
    my $infile_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $lane_href;
    my $program_name;
    my $sample_info_href;

    ## Default(s)
    my $outaligner_dir_ref;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        indir_path_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$indir_path_href,
            strict_type => 1,
        },
        infile_both_strands_prefix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_both_strands_prefix_href,
            strict_type => 1,
        },
        infile_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_href,
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
        lane_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$lane_href,
            strict_type => 1,
        },
        outaligner_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{outaligner_dir},
            store       => \$outaligner_dir_ref,
            strict_type => 1,
        },
        program_name => {
            defined     => 1,
            required    => 1,
            store       => \$program_name,
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

    use MIP::Check::Parameter qw{ check_gzipped };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

# Used to decide later if any inputfiles needs to be compressed before starting analysis
    my $uncompressed_file_counter = 0;

  SAMPLE_ID:
    for my $sample_id ( keys %{$infile_href} ) {

        # Needed to be able to track when lanes are finished
        my $lane_tracker = 0;

      INFILE:
        while ( my ( $file_index, $file_name ) =
            each( @{ $infile_href->{$sample_id} } ) )
        {

            ## Check if a file is gzipped.
            my $compressed_switch =
              check_gzipped( { file_name => $file_name, } );
            my $read_file_command = q{zcat};

            ## Not compressed
            if ( not $compressed_switch ) {

                ## File needs compression before starting analysis. Note: All files are rechecked downstream and uncompressed ones are gzipped automatically
                $uncompressed_file_counter = q{uncompressed};
                $read_file_command         = q{cat};
            }

            ## Parse 'new' no "index" format $1=lane, $2=date,
            ## $3=Flow-cell, $4=Sample_id, $5=index,$6=direction
            if (
                $file_name =~ /(\d+)_(\d+)_([^_]+)_([^_]+)_([^_]+)_(\d).fastq/ )
            {

                ## Check that the sample_id provided and sample_id in infile name match.
                check_sample_id_match(
                    {
                        active_parameter_href => $active_parameter_href,
                        file_index            => $file_index,
                        infile_href           => $infile_href,
                        infile_sample_id => $4,    #$4 = Sample_id from filename
                        sample_id => $sample_id,
                    }
                );

                ## Adds information derived from infile name to sample_info hash. Tracks the number of lanes sequenced and checks unique array elementents.
                add_infile_info(
                    {
                        active_parameter_href => $active_parameter_href,
                        compressed_switch     => $compressed_switch,
                        date                  => $2,
                        direction             => $6,
                        file_info_href        => $file_info_href,
                        file_index            => $file_index,
                        flowcell              => $3,
                        index                 => $5,
                        indir_path_href       => $indir_path_href,
                        infile_both_strands_prefix_href =>
                          $infile_both_strands_prefix_href,
                        infile_href             => $infile_href,
                        infile_lane_prefix_href => $infile_lane_prefix_href,
                        lane                    => $1,
                        lane_href               => $lane_href,
                        lane_tracker_ref        => \$lane_tracker,
                        sample_id               => $4,
                        sample_info_href        => $sample_info_href,
                    }
                );
            }
            else
            {    #No regexp match i.e. file does not follow filename convention

                $log->warn(
                        q{Could not detect MIP file name convention for file: }
                      . $file_name
                      . q{.} );
                $log->warn(
                    q{Will try to find mandatory information in fastq header.});

                ## Check that file name at least contains sample_id
                if ( $file_name !~ /$sample_id/ ) {

                    $log->fatal(
q{Please check that the file name contains the sample_id.}
                    );
                }

                ## Get run info from fastq file header
                my @fastq_info_headers = get_run_info(
                    {
                        directory         => $indir_path_href->{$sample_id},
                        file              => $file_name,
                        read_file_command => $read_file_command,
                    }
                );

                ## Adds information derived from infile name to sample_info hash. Tracks the number of lanes sequenced and checks unique array elementents.
                add_infile_info(
                    {
                        active_parameter_href => $active_parameter_href,
                        compressed_switch     => $compressed_switch,
                        ## fastq format does not contain a date of the run,
                        ## so fake it with constant impossible date
                        date            => q{000101},
                        direction       => $fastq_info_headers[4],
                        file_index      => $file_index,
                        file_info_href  => $file_info_href,
                        flowcell        => $fastq_info_headers[2],
                        index           => $fastq_info_headers[5],
                        indir_path_href => $indir_path_href,
                        infile_both_strands_prefix_href =>
                          $infile_both_strands_prefix_href,
                        infile_href             => $infile_href,
                        infile_lane_prefix_href => $infile_lane_prefix_href,
                        lane                    => $fastq_info_headers[3],
                        lane_href               => $lane_href,
                        lane_tracker_ref        => \$lane_tracker,
                        sample_id               => $sample_id,
                        sample_info_href        => $sample_info_href,
                    }
                );

                $log->info(
                        q{Found following information from fastq header: lane=}
                      . $fastq_info_headers[3]
                      . q{ flow-cell=}
                      . $fastq_info_headers[2]
                      . q{ index=}
                      . $fastq_info_headers[5]
                      . q{ direction=}
                      . $fastq_info_headers[4],
                );
                $log->warn(
q{Will add fake date '20010101' to follow file convention since this is not recorded in fastq header}
                );
            }
        }
    }
    return $uncompressed_file_counter;
}

sub check_sample_id_match {

##check_sample_id_match

##Function : Check that the sample_id provided and sample_id in infile name match.
##Returns  : ""
##Arguments: $active_parameter_href, $infile_href, $sample_id, $infile_sample_id, $file_index
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $infile_href           => Infiles hash {REF}
##         : $sample_id             => Sample id from user
##         : $infile_sample_id      => Sample_id collect with regexp from infile
##         : $file_index            => Counts the number of infiles

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $infile_href;
    my $sample_id;
    my $infile_sample_id;
    my $file_index;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        infile_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_href
        },
        sample_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sample_id,
        },
        infile_sample_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_sample_id
        },
        file_index => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$file_index
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my %seen = ( $infile_sample_id => 1 );    #Add input as first increment

    foreach my $sample_id_supplied ( @{ $active_parameter_href->{sample_ids} } )
    {

        $seen{$sample_id_supplied}++;
    }
    unless ( $seen{$infile_sample_id} > 1 ) {

        $log->fatal( $sample_id
              . " supplied and sample_id "
              . $infile_sample_id
              . " found in file : "
              . $infile_href->{$sample_id}[$file_index]
              . " does not match. Please rename file to match sample_id: "
              . $sample_id
              . "\n" );
        exit 1;
    }
}

sub get_run_info {

##get_run_info

##Function : Get run info from fastq file header
##Returns  : ""
##Arguments: $directory, $read_file, $file
##         : $directory       => Directory of file
##         : $read_file_command => Command used to read file
##         : $file            => File to parse

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $directory;
    my $read_file_command;
    my $file;

    my $tmpl = {
        directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$directory
        },
        read_file_command => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$read_file_command
        },
        file =>
          { required => 1, defined => 1, strict_type => 1, store => \$file },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my $fastq_header_regexp =
q?perl -nae 'chomp($_); if($_=~/^(@\w+):(\w+):(\w+):(\w+)\S+\s(\w+):\w+:\w+:(\w+)/) {print $1." ".$2." ".$3." ".$4." ".$5." ".$6."\n";} if($.=1) {last;}' ?;

    my $pwd = cwd();      #Save current direcory
    chdir($directory);    #Move to sample_id infile directory

    my $fastq_info_headers = `$read_file_command $file | $fastq_header_regexp;`
      ;                   #Collect fastq header info
    my @fastq_info_headers = split( " ", $fastq_info_headers );

    chdir($pwd);          #Move back to original directory

    unless ( scalar(@fastq_info_headers) eq 6 ) {

        $log->fatal(
"Could not detect reuired sample sequencing run info from fastq file header - PLease proved MIP file in MIP file convention format to proceed\n"
        );
        exit 1;
    }

    return @fastq_info_headers;
}

sub add_infile_info {

##add_infile_info

##Function : Adds information derived from infile name to sample_info hash. Tracks the number of lanes sequenced and checks unique array elementents.
##Returns  : ""
##Arguments: $active_parameter_href, $sample_info_href, $file_info_href, $infile_href, $infile_lane_prefix_href, $infile_both_strands_prefix_href, $indir_path_href, $lane_href, $lane, $date, $flowcell, $sample_id, $index, $direction, $lane_tracker_ref, $file_index, $compressed_switch
##         : $active_parameter_href              => Active parameters for this analysis hash {REF}
##         : $sample_info_href                   => Info on samples and family hash {REF}
##         : $file_info_href                     => File info hash {REF}
##         : $infile_href                        => Infiles hash {REF}
##         : $infile_lane_prefix_href         => Infile(s) without the ".ending" {REF}
##         : $infile_both_strands_prefix_href => The infile(s) without the ".ending" and strand info {REF}
##         : $indir_path_href                    => Indirectories path(s) hash {REF}
##         : $lane_href                          => The lane info hash {REF}
##         : $lane                               => Flow-cell lane
##         : $date                               => Flow-cell sequencing date
##         : $flowcell                           => Flow-cell id
##         : $sample_id                          => Sample id
##         : $index                              => The DNA library preparation molecular barcode
##         : $direction                          => Sequencing read direction
##         : $lane_tracker_ref                   => Counts the number of lanes sequenced {REF}
##         : $file_index                         => Index of file
##         : $compressed_switch                  => ".fastq.gz" or ".fastq" info governs zcat or cat downstream

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_href;
    my $indir_path_href;
    my $infile_lane_prefix_href;
    my $infile_both_strands_prefix_href;
    my $lane_href;
    my $lane_tracker_ref;
    my $sample_id;
    my $lane;
    my $date;
    my $flowcell;
    my $index;
    my $direction;
    my $file_index;
    my $compressed_switch;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        infile_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_href
        },
        indir_path_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$indir_path_href
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        infile_both_strands_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_both_strands_prefix_href
        },
        lane_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$lane_href
        },
        sample_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sample_id,
        },
        lane => {
            required    => 1,
            defined     => 1,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$lane
        },
        lane_tracker_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$lane_tracker_ref
        },
        date =>
          { required => 1, defined => 1, strict_type => 1, store => \$date },
        flowcell => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$flowcell
        },
        index =>
          { required => 1, defined => 1, strict_type => 1, store => \$index },
        direction => {
            required    => 1,
            defined     => 1,
            allow       => [ 1, 2 ],
            strict_type => 1,
            store       => \$direction
        },
        file_index => {
            required    => 1,
            defined     => 1,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$file_index
        },
        compressed_switch => {
            required    => 1,
            defined     => 1,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$compressed_switch
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $read_file;
    my $file_at_lane_level_ref;
    my $file_at_direction_level_ref;

    my $parsed_date = Time::Piece->strptime( $date, "%y%m%d" );
    $parsed_date = $parsed_date->ymd;

    if ($compressed_switch) {

        $read_file = "zcat";    #Read file in compressed format
    }
    else {

        $read_file = "cat";     #Read file in uncompressed format
    }

    if ( $direction == 1 ) {    #Read 1

        push( @{ $lane_href->{$sample_id} }, $lane );    #Lane
        $infile_lane_prefix_href->{$sample_id}[$$lane_tracker_ref] =
            $sample_id . "."
          . $date . "_"
          . $flowcell . "_"
          . $index . ".lane"
          . $lane
          ; #Save new format (sample_id_date_flow-cell_index_lane) in hash with samplid as keys and inputfiles in array. Note: These files have not been created yet and there is one entry into hash for both strands and .ending is removed (.fastq).

        $file_at_lane_level_ref =
          \$infile_lane_prefix_href->{$sample_id}[$$lane_tracker_ref];    #Alias
        $sample_info_href->{sample}{$sample_id}{file}{$$file_at_lane_level_ref}
          {sequence_run_type} = "single_end"; #Single_end until proven otherwise

        ## Collect read length from an infile
        $sample_info_href->{sample}{$sample_id}{file}{$$file_at_lane_level_ref}
          {sequence_length} = collect_read_length(
            {
                directory         => $indir_path_href->{$sample_id},
                read_file_command => $read_file,
                file              => $infile_href->{$sample_id}[$file_index],
            }
          );

        ## Check if fastq file is interleaved
        $sample_info_href->{sample}{$sample_id}{file}{$$file_at_lane_level_ref}
          {interleaved} = detect_interleaved(
            {
                directory         => $indir_path_href->{$sample_id},
                read_file_command => $read_file,
                file              => $infile_href->{$sample_id}[$file_index],
            }
          );

        ## Detect "regexp" in string
        $file_info_href->{undetermined_in_file_name}
          { $infile_lane_prefix_href->{$sample_id}[$$lane_tracker_ref] } =
          check_string(
            {
                string => $flowcell,
                regexp => "Undetermined",
            }
          );
        $$lane_tracker_ref++;
    }
    if ( $direction == 2 ) {    #2nd read direction

        $file_at_lane_level_ref =
          \$infile_lane_prefix_href->{$sample_id}[ $$lane_tracker_ref - 1 ]
          ;                     #Alias
        $sample_info_href->{sample}{$sample_id}{file}{$$file_at_lane_level_ref}
          {sequence_run_type} = 'paired-end'
          ;    #$lane_tracker -1 since it gets incremented after direction eq 1.
    }

    $infile_both_strands_prefix_href->{$sample_id}[$file_index] =
        $sample_id . "."
      . $date . "_"
      . $flowcell . "_"
      . $index . ".lane"
      . $lane . "_"
      . $direction
      ; #Save new format in hash with samplid as keys and inputfiles in array. Note: These files have not been created yet and there is one entry per strand and .ending is removed (.fastq).

    $file_at_direction_level_ref =
      \$infile_both_strands_prefix_href->{$sample_id}[$file_index];    #Alias
    $sample_info_href->{sample}{$sample_id}{file}{$$file_at_lane_level_ref}
      {read_direction_file}{$$file_at_direction_level_ref}{original_file_name}
      = $infile_href->{$sample_id}[$file_index];    #Original file_name

    $sample_info_href->{sample}{$sample_id}{file}{$$file_at_lane_level_ref}
      {read_direction_file}{$$file_at_direction_level_ref}
      {original_file_name_prefix} =
        $lane . "_"
      . $date . "_"
      . $flowcell . "_"
      . $sample_id . "_"
      . $index . "_"
      . $direction;    #Original file_name, but no ending

    $sample_info_href->{sample}{$sample_id}{file}{$$file_at_lane_level_ref}
      {read_direction_file}{$$file_at_direction_level_ref}{lane} =
      $lane;           #Save sample lane

    $sample_info_href->{sample}{$sample_id}{file}{$$file_at_lane_level_ref}
      {read_direction_file}{$$file_at_direction_level_ref}{date} =
      $parsed_date;    #Save Sequence run date

    $sample_info_href->{sample}{$sample_id}{file}{$$file_at_lane_level_ref}
      {read_direction_file}{$$file_at_direction_level_ref}{flowcell} =
      $flowcell;       #Save Sequence flow-cell

    $sample_info_href->{sample}{$sample_id}{file}{$$file_at_lane_level_ref}
      {read_direction_file}{$$file_at_direction_level_ref}{sample_barcode} =
      $index;          #Save sample barcode

    $sample_info_href->{sample}{$sample_id}{file}{$$file_at_lane_level_ref}
      {read_direction_file}{$$file_at_direction_level_ref}{run_barcode} =
      $date . "_" . $flowcell . "_" . $lane . "_" . $index;    #Save run barcode

    $sample_info_href->{sample}{$sample_id}{file}{$$file_at_lane_level_ref}
      {read_direction_file}{$$file_at_direction_level_ref}{read_direction} =
      $direction;
}

sub detect_interleaved {

##detect_interleaved

##Function : Detect if fastq file is interleaved
##Returns  : "1(=interleaved)"
##Arguments: $directory, $read_file, $file
##         : $directory         => Directory of file
##         : $read_file_command => Command used to read file
##         : $file              => File to parse

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $directory;
    my $read_file_command;
    my $file;

    my $tmpl = {
        directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$directory
        },
        read_file_command => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$read_file_command
        },
        file =>
          { required => 1, defined => 1, strict_type => 1, store => \$file },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my $interleaved_regexp =
q?perl -nae 'chomp($_); if( ($_=~/^@\S+:\w+:\w+:\w+\S+\s(\w+):\w+:\w+:\w+/) && ($.==5) ) {print $1."\n";last;} elsif ($.==6) {last;}' ?;

    my $pwd = cwd();      #Save current direcory
    chdir($directory);    #Move to sample_id infile directory

    my $fastq_info_headers = `$read_file_command $file | $interleaved_regexp;`
      ;                   #Collect interleaved info

    if ( !$fastq_info_headers ) {

        my $interleaved_regexp =
q?perl -nae 'chomp($_); if( ($_=~/^@\w+-\w+:\w+:\w+:\w+:\w+:\w+:\w+\/(\w+)/) && ($.==5) ) {print $1."\n";last;} elsif ($.==6) {last;}' ?;
        $fastq_info_headers = `$read_file_command $file | $interleaved_regexp;`
          ;               #Collect interleaved info
    }

    chdir($pwd);          #Move back to original directory

    unless ( $fastq_info_headers =~ /[1, 2, 3]/ ) {

        $log->fatal("Malformed fastq file!\n");
        $log->fatal( "Read direction is: "
              . $fastq_info_headers
              . " allowed entries are '1', '2', '3'. Please check fastq file\n"
        );
        exit 1;
    }
    if ( $fastq_info_headers > 1 ) {

        $log->info( "Found interleaved fastq file: " . $file, "\n" );
        return 1;
    }
    return;
}

sub create_file_endings {

## Function : Creates the file_tags depending on which modules are used by the user to relevant chain.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $family_id_ref           => Family id {REF}
##          : $file_info_href          => Info on files hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $order_parameters_ref    => Order of addition to parameter array {REF}
##          : $parameter_href          => Parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $order_parameters_ref;
    my $parameter_href;

    ## Default(s)
    my $family_id_ref;

    my $tmpl = {
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
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
        order_parameters_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$order_parameters_ref,
            strict_type => 1,
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            store       => \$family_id_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};

    ## Used to enable seqential build-up of file_tags between modules
    my %temp_file_ending;

  PARAMETER:
    foreach my $order_parameter_element (@$order_parameters_ref) {

        ## Only active parameters
        if ( defined $active_parameter_href->{$order_parameter_element} ) {

            ## Only process programs
            if (
                any { $_ eq $order_parameter_element }
                @{ $parameter_href->{dynamic_parameter}{program} }
              )
            {

                ## MAIN chain
                if ( $parameter_href->{$order_parameter_element}{chain} eq
                    q{MAIN} )
                {

                    ##  File_tag exist
                    if ( $parameter_href->{$order_parameter_element}{file_tag}
                        ne q{nofile_tag} )
                    {

                        ## Alias
                        my $file_ending_ref =
                          \$parameter_href->{$order_parameter_element}
                          {file_tag};

                        ###MAIN/Per sample_id
                      SAMPLE_ID:
                        foreach my $sample_id (
                            @{ $active_parameter_href->{sample_ids} } )
                        {

                            ## File_ending should be added
                            if ( $active_parameter_href
                                ->{$order_parameter_element} > 0 )
                            {

                                ## Special case
                                if ( $order_parameter_element eq
                                    q{ppicardtools_mergesamfiles} )
                                {

                                    $file_info_href->{$sample_id}
                                      {ppicardtools_mergesamfiles}{file_tag} =
                                      $temp_file_ending{$sample_id} . "";
                                }
                                else {

                                    if ( defined $temp_file_ending{$sample_id} )
                                    {

                                        $file_info_href->{$sample_id}
                                          {$order_parameter_element}{file_tag}
                                          = $temp_file_ending{$sample_id}
                                          . $$file_ending_ref;
                                    }
                                    else {
                                        ## First module that should add filending

                                        $file_info_href->{$sample_id}
                                          {$order_parameter_element}{file_tag}
                                          = $$file_ending_ref;
                                    }
                                }
                            }
                            else {
                                ## Do not add new module file_tag

                                $file_info_href->{$sample_id}
                                  {$order_parameter_element}{file_tag} =
                                  $temp_file_ending{$sample_id};
                            }

                            ## To enable sequential build-up of fileending
                            $temp_file_ending{$sample_id} =
                              $file_info_href->{$sample_id}
                              {$order_parameter_element}{file_tag};
                        }

                        ###MAIN/Per family_id
                        ## File_ending should be added
                        if ( $active_parameter_href->{$order_parameter_element}
                            > 0 )
                        {

                            ## Special case - do nothing
                            if ( $order_parameter_element eq
                                q{ppicardtools_mergesamfiles} )
                            {
                            }
                            else {

                                if (
                                    defined $temp_file_ending{$$family_id_ref} )
                                {

                                    $file_info_href->{$$family_id_ref}
                                      {$order_parameter_element}{file_tag} =
                                        $temp_file_ending{$$family_id_ref}
                                      . $$file_ending_ref;
                                }
                                else {
                                    ## First module that should add filending

                                    $file_info_href->{$$family_id_ref}
                                      {$order_parameter_element}{file_tag} =
                                      $$file_ending_ref;
                                }

                                ## To enable sequential build-up of fileending
                                $temp_file_ending{$$family_id_ref} =
                                  $file_info_href->{$$family_id_ref}
                                  {$order_parameter_element}{file_tag};
                            }
                        }
                        else {
                            ## Do not add new module file_tag

                            $file_info_href->{$$family_id_ref}
                              {$order_parameter_element}{file_tag} =
                              $temp_file_ending{$$family_id_ref};
                        }
                    }
                }

                ## Other chain(s)
                if ( $parameter_href->{$order_parameter_element}{chain} ne
                    q{MAIN} )
                {

                    ## Alias
                    my $chain_fork =
                      $parameter_href->{$order_parameter_element}{chain};

                    ## File_tag exist
                    if ( $parameter_href->{$order_parameter_element}{file_tag}
                        ne q{nofile_tag} )
                    {

                        ## Alias
                        my $file_ending_ref =
                          \$parameter_href->{$order_parameter_element}
                          {file_tag};

                        ###OTHER/Per sample_id
                      SAMPLE_ID:
                        foreach my $sample_id (
                            @{ $active_parameter_href->{sample_ids} } )
                        {

                            ## File_ending should be added
                            if ( $active_parameter_href
                                ->{$order_parameter_element} > 0 )
                            {

                                if (
                                    not
                                    defined $temp_file_ending{$chain_fork}
                                    {$sample_id} )
                                {

                                    ## Inherit current MAIN chain.
                                    $temp_file_ending{$chain_fork}{$sample_id}
                                      = $temp_file_ending{$sample_id};
                                }
                                if (
                                    defined $temp_file_ending{$chain_fork}
                                    {$sample_id} )
                                {

                                    $file_info_href->{$sample_id}
                                      {$order_parameter_element}{file_tag} =
                                      $temp_file_ending{$chain_fork}{$sample_id}
                                      . $$file_ending_ref;
                                }
                                else {
                                    ## First module that should add filending

                                    $file_info_href->{$sample_id}
                                      {$order_parameter_element}{file_tag} =
                                      $$file_ending_ref;
                                }
                            }
                            else {
                                ## Do not add new module file_tag

                                $file_info_href->{$sample_id}
                                  {$order_parameter_element}{file_tag} =
                                  $temp_file_ending{$chain_fork}{$sample_id};
                            }

                            ## To enable sequential build-up of fileending
                            $temp_file_ending{$chain_fork}{$sample_id} =
                              $file_info_href->{$sample_id}
                              {$order_parameter_element}{file_tag};
                        }
                        ###Other/Per family_id

                        ## File ending should be added
                        if ( $active_parameter_href->{$order_parameter_element}
                            > 0 )
                        {

                            if (
                                not defined $temp_file_ending{$chain_fork}
                                {$$family_id_ref} )
                            {

                                ## Inherit current MAIN chain.
                                $temp_file_ending{$chain_fork}{$$family_id_ref}
                                  = $temp_file_ending{$$family_id_ref};
                            }
                            if (
                                defined $temp_file_ending{$chain_fork}
                                {$$family_id_ref} )
                            {

                                $file_info_href->{$$family_id_ref}
                                  {$order_parameter_element}{file_tag} =
                                  $temp_file_ending{$chain_fork}
                                  {$$family_id_ref} . $$file_ending_ref;
                            }
                            else {
                                ## First module that should add filending

                                $file_info_href->{$$family_id_ref}
                                  {$order_parameter_element}{file_tag} =
                                  $$file_ending_ref;
                            }

                            ## To enable sequential build-up of fileending
                            $temp_file_ending{$chain_fork}{$$family_id_ref} =
                              $file_info_href->{$$family_id_ref}
                              {$order_parameter_element}{file_tag};
                        }
                        else {
                            ## Do not add new module file_tag

                            $file_info_href->{$$family_id_ref}
                              {$order_parameter_element}{file_tag} =
                              $temp_file_ending{$chain_fork}{$$family_id_ref};
                        }
                    }
                }
            }
        }
    }
    return;
}

sub write_cmd_mip_log {

##write_cmd_mip_log

##Function : Write CMD to MIP log file
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, $order_parameters_ref, $script_ref, $log_file_ref
##         : $parameter_href        => Parameter hash {REF}
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $order_parameters_ref  => Order of addition to parameter array {REF}
##         : $script_ref            => The script that is being executed {REF}
##         : $log_file_ref          => The log file {REF}
##         : $mip_version_ref       => The MIP version

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $order_parameters_ref;
    my $script_ref;
    my $log_file_ref;
    my $mip_version_ref;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        order_parameters_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$order_parameters_ref
        },
        script_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$script_ref
        },
        log_file_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$log_file_ref
        },
        mip_version_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$mip_version_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my $cmd_line = $$script_ref . " ";

    my @nowrite = (
        "mip",                  "bwa_build_reference",
        "pbamcalibrationblock", "pvariantannotationblock",
        q{associated_program},  q{rtg_build_reference},
    );

  PARAMETER_KEY:
    foreach my $order_parameter_element ( @{$order_parameters_ref} ) {

        if ( defined $active_parameter_href->{$order_parameter_element} ) {

            ## If no config file do not print
            if (   $order_parameter_element eq q{config_file}
                && $active_parameter_href->{config_file} eq 0 )
            {
            }
            else {

                ## If element is part of array - do nothing
                if ( any { $_ eq $order_parameter_element } @nowrite ) {
                }
                elsif (
                    ## Array reference
                    (
                        exists $parameter_href->{$order_parameter_element}
                        {data_type}
                    )
                    && ( $parameter_href->{$order_parameter_element}{data_type}
                        eq q{ARRAY} )
                  )
                {

                    my $separator = $parameter_href->{$order_parameter_element}
                      {element_separator};
                    $cmd_line .= "-"
                      . $order_parameter_element . " "
                      . join(
                        $separator,
                        @{
                            $active_parameter_href->{$order_parameter_element}
                        }
                      ) . " ";
                }
                elsif (
                    ## HASH reference
                    (
                        exists $parameter_href->{$order_parameter_element}
                        {data_type}
                    )
                    && ( $parameter_href->{$order_parameter_element}{data_type}
                        eq q{HASH} )
                  )
                {

                    # First key
                    $cmd_line .= "-" . $order_parameter_element . " ";
                    $cmd_line .= join(
                        "-" . $order_parameter_element . " ",
                        map {
"$_=$active_parameter_href->{$order_parameter_element}{$_} "
                        } (
                            keys %{
                                $active_parameter_href
                                  ->{$order_parameter_element}
                            }
                        )
                    );
                }
                else {

                    $cmd_line .= "-"
                      . $order_parameter_element . " "
                      . $active_parameter_href->{$order_parameter_element}
                      . " ";
                }
            }
        }
    }
    $log->info( $cmd_line,                            "\n" );
    $log->info( q{MIP Version: } . $$mip_version_ref, "\n" );
    $log->info(
        q{Script parameters and info from }
          . $$script_ref
          . q{ are saved in file: }
          . $$log_file_ref,
        "\n"
    );
    return;
}

sub size_sort_select_file_contigs {

## Function : Sorts array depending on reference array. NOTE: Only entries present in reference array will survive in sorted array.
## Returns  : @sorted_contigs
## Arguments: $consensus_analysis_type_ref => Consensus analysis_type {REF}
##          : $file_info_href              => File info hash {REF}
##          : $hash_key_sort_reference     => The hash keys sort reference
##          : $hash_key_to_sort            => The keys to sort

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $consensus_analysis_type_ref;
    my $file_info_href;
    my $hash_key_sort_reference;
    my $hash_key_to_sort;

    my $tmpl = {
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        consensus_analysis_type_ref => {
            default     => \$$,
            defined     => 1,
            required    => 1,
            store       => \$consensus_analysis_type_ref,
            strict_type => 1,
        },
        hash_key_to_sort => {
            defined     => 1,
            required    => 1,
            store       => \$hash_key_to_sort,
            strict_type => 1,
        },
        hash_key_sort_reference => {
            defined     => 1,
            required    => 1,
            store       => \$hash_key_sort_reference,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::Hash qw{ check_element_exist_hash_of_array };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my @sorted_contigs;

    ## Sort the contigs depending on reference array
    if ( $file_info_href->{$hash_key_to_sort} ) {

        foreach my $element ( @{ $file_info_href->{$hash_key_sort_reference} } )
        {

            if (
                not check_element_exist_hash_of_array(
                    {
                        element  => $element,
                        hash_ref => $file_info_href,
                        key      => $hash_key_to_sort,
                    }
                )
              )
            {

                push @sorted_contigs, $element;
            }
        }
    }

    ## Test if all contigs collected from select file was sorted by reference contig array
    if ( @sorted_contigs
        && scalar @{ $file_info_href->{$hash_key_to_sort} } !=
        scalar @sorted_contigs )
    {

        foreach my $element ( @{ $file_info_href->{$hash_key_to_sort} } ) {

            ## If element is not part of array
            if ( not any { $_ eq $element } @sorted_contigs ) {

                ## Special case when analysing wes since Mitochondrial contigs have no baits in exome capture kits
                unless ( $$consensus_analysis_type_ref eq q{wes}
                    && $element =~ /MT$|M$/ )
                {

                    $log->fatal( q{Could not detect '##contig'= }
                          . $element
                          . q{ from meta data header in '-vcfparser_select_file' in reference contigs collected from '-human_genome_reference'}
                    );
                    exit 1;
                }
            }
        }
    }
    return @sorted_contigs;
}

sub add_to_sample_info {

##add_to_sample_info

##Function : Adds parameter info to sample_info
##Returns  : ""
##Arguments: $active_parameter_href, $sample_info_href, $file_info_href, $family_id_ref
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $sample_info_href      => Info on samples and family hash {REF}
##         : $file_info_href        => File info hash {REF}
##         : $family_id_ref         => The family_id_ref {REF}

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;
    my $human_genome_reference_ref;
    my $outdata_dir;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
        human_genome_reference_ref => {
            default =>
              \$arg_href->{active_parameter_href}{human_genome_reference},
            strict_type => 1,
            store       => \$human_genome_reference_ref
        },
        outdata_dir => {
            default     => $arg_href->{active_parameter_href}{outdata_dir},
            strict_type => 1,
            store       => \$outdata_dir
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::QC::Record qw(add_program_outfile_to_sample_info);

    if ( exists( $active_parameter_href->{analysis_type} ) ) {

        $sample_info_href->{analysis_type} =
          $active_parameter_href->{analysis_type};
    }
    if ( exists( $active_parameter_href->{expected_coverage} ) ) {

        $sample_info_href->{expected_coverage} =
          $active_parameter_href->{expected_coverage};
    }
    if ( exists $active_parameter_href->{sample_origin} ) {

        $sample_info_href->{sample_origin} =
          $active_parameter_href->{sample_origin};
    }
    if ( exists $active_parameter_href->{gatk_path}
        && $active_parameter_href->{gatk_path} )
    {

        my $gatk_version;
        if ( $active_parameter_href->{gatk_path} =~ /GenomeAnalysisTK-([^,]+)/ )
        {

            $gatk_version = $1;
        }
        else {
            # Fall back on actually calling program

            my $jar_path = catfile( $active_parameter_href->{gatk_path},
                "GenomeAnalysisTK.jar" );
            $gatk_version = (`java -jar $jar_path --version 2>&1`);
            chomp $gatk_version;
        }
        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                program_name     => 'gatk',
                version          => $gatk_version,
            }
        );
    }
    if ( exists $active_parameter_href->{picardtools_path}
        && $active_parameter_href->{picardtools_path} )
    {
        ## To enable addition of version to sample_info
        my $picardtools_version;
        if ( $active_parameter_href->{picardtools_path} =~
            /picard-tools-([^,]+)/ )
        {

            $picardtools_version = $1;
        }
        else {    #Fall back on actually calling program

            my $jar_path = catfile( $active_parameter_href->{picardtools_path},
                q{picard.jar} );
            $picardtools_version =
              (`java -jar $jar_path CreateSequenceDictionary --version 2>&1`);
            chomp $picardtools_version;
        }

        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                program_name     => 'picardtools',
                version          => $picardtools_version,
            }
        );
    }
    my @sambamba_programs =
      ( "pbwa_mem", "psambamba_depth", "markduplicates_sambamba_markdup" );
    foreach my $program (@sambamba_programs) {

        if (   ( defined $active_parameter_href->{$program} )
            && ( $active_parameter_href->{$program} == 1 ) )
        {

            if ( !$active_parameter_href->{dry_run_all} ) {

                my $regexp =
                  q?perl -nae 'if($_=~/sambamba\s(\S+)/) {print $1;last;}'?;
                my $sambamba_version = (`sambamba 2>&1 | $regexp`);
                chomp $sambamba_version;
                add_program_outfile_to_sample_info(
                    {
                        sample_info_href => $sample_info_href,
                        program_name     => 'sambamba',
                        version          => $sambamba_version,
                    }
                );
                last;    #Only need to check once
            }
        }
    }
    if ( exists $active_parameter_href->{pcnvnator} )
    {                    #To enable addition of version to sample_info

        if (   ( $active_parameter_href->{pcnvnator} == 1 )
            && ( !$active_parameter_href->{dry_run_all} ) )
        {

            my $regexp =
              q?perl -nae 'if($_=~/CNVnator\s+(\S+)/) {print $1;last;}'?;
            my $cnvnator_version = (`cnvnator 2>&1 | $regexp`);
            chomp $cnvnator_version;
            add_program_outfile_to_sample_info(
                {
                    sample_info_href => $sample_info_href,
                    program_name     => 'cnvnator',
                    version          => $cnvnator_version,
                }
            );
        }
    }
    if ( defined($$human_genome_reference_ref) )
    {    #To enable addition of version to sample_info

        $sample_info_href->{human_genome_build}{path} =
          $$human_genome_reference_ref;
        $sample_info_href->{human_genome_build}{source} =
          $file_info_href->{human_genome_reference_source};
        $sample_info_href->{human_genome_build}{version} =
          $file_info_href->{human_genome_reference_version};
    }
    if ( exists( $active_parameter_href->{pedigree_file} ) ) {

        ## Add pedigree_file to sample_info
        $sample_info_href->{pedigree_file}{path} =
          $active_parameter_href->{pedigree_file};
    }
    if ( exists( $active_parameter_href->{log_file} ) ) {

        my $path = dirname( dirname( $active_parameter_href->{log_file} ) );
        $sample_info_href->{log_file_dir} =
          $path;    #Add log_file_dir to SampleInfoFile
        $sample_info_href->{last_log_file_path} =
          $active_parameter_href->{log_file};
    }
}

sub check_string {

##check_string

##Function : Detect "regexp" in string
##Returns  : ""|1
##Arguments: $string
##         : $string => String to be searched
##         : $regexp => regexp to use on string

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $string;
    my $regexp;

    my $tmpl = {
        string =>
          { required => 1, defined => 1, strict_type => 1, store => \$string },
        regexp =>
          { required => 1, defined => 1, strict_type => 1, store => \$regexp },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    if ( $string =~ /$regexp/ ) {

        return 1;
    }
}

sub collect_read_length {

## Function : Collect read length from an infile
## Returns  : "readLength"
## Arguments: $directory => Directory of file
##          : $file      => File to parse
##          : $read_file => Command used to read file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $directory;
    my $file;
    my $read_file_command;

    my $tmpl = {
        directory => {
            defined     => 1,
            required    => 1,
            store       => \$directory,
            strict_type => 1,
        },
        read_file_command => {
            defined     => 1,
            required    => 1,
            store       => \$read_file_command,
            strict_type => 1,
        },
        file =>
          { defined => 1, required => 1, store => \$file, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Prints sequence length and exits
    my $seq_length_regexp =
q?perl -ne 'if ($_!~/@/) {chomp($_);my $seq_length = length($_);print $seq_length;last;}' ?;

    ## Save current direcory
    my $pwd = cwd();

    ## Move to sample_id infile directory
    chdir($directory);

    ## Collect sequence length
    my $ret = `$read_file_command $file | $seq_length_regexp;`;

    ## Move to original directory
    chdir($pwd);

    return $ret;
}

sub get_matching_values_key {

##get_matching_values_key

##Function : Return the key if the hash value and query match
##Returns  : "key pointing to matched value"
##Arguments: $active_parameter_href, $query_value_ref, $parameter_name
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $query_value_ref       => The value to query in the hash {REF}
##         : $parameter_name        => MIP parameter name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $query_value_ref;
    my $parameter_name;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        query_value_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$query_value_ref
        },
        parameter_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$parameter_name,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my %reversed = reverse %{ $active_parameter_href->{$parameter_name} }
      ;    #Values are now keys and vice versa

    if ( exists $reversed{$$query_value_ref} ) {

        return $reversed{$$query_value_ref};
    }
}

##Investigate potential autodie error
if ( $@ and $@->isa("autodie::exception") ) {

    if ( $@->matches("default") ) {

        say "Not an autodie error at all";
    }
    if ( $@->matches("open") ) {

        say "Error from open";
    }
    if ( $@->matches(":io") ) {

        say "Non-open, IO error.\n";
    }
}
elsif ($@) {

    say "A non-autodie exception.";
}

1;
