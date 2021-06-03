package MIP::Recipes::Analysis::Rhocall;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile devnull splitpath };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

# MIPs lib/
use MIP::Constants qw{ $ASTERISK $DOT $LOG_NAME $NEWLINE $PIPE $SEMICOLON $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_rhocall_annotate analysis_rhocall_viz };

}

sub analysis_rhocall_annotate {

## Function : Rhocall performs annotation of autozygosity regions
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $file_path               => File path
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Program name
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory {REF}
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $file_path;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_info_href;

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
        file_path   => { store => \$file_path, strict_type => 1, },
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
        xargs_file_counter => {
            allow       => qr/ ^\d+$ /xsm,
            default     => 0,
            store       => \$xargs_file_counter,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Cluster qw{ get_core_number update_memory_allocation };
    use MIP::File_info qw{ get_io_files parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Bcftools qw{ bcftools_roh };
    use MIP::Program::Rhocall qw{ rhocall_annotate };
    use MIP::Sample_info qw{ set_recipe_outfile_in_sample_info };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
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
    my %infile_path        = %{ $io{in}{file_path_href} };

    my $consensus_analysis_type = $parameter_href->{cache}{consensus_analysis_type};
    my @contigs_size_ordered    = @{ $file_info_href->{contigs_size_ordered} };
    my %recipe                  = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            recipe_name           => $recipe_name,
        }
    );
    my $core_number = $recipe{core_number};

    ## Set and get the io files per chain, id and stream
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id         => $recipe{job_id_chain},
                id               => $case_id,
                file_info_href   => $file_info_href,
                file_name_prefix => $infile_name_prefix,
                iterators_ref    => \@contigs_size_ordered,
                outdata_dir      => $active_parameter_href->{outdata_dir},
                parameter_href   => $parameter_href,
                recipe_name      => $recipe_name,
                temp_directory   => $temp_directory,
            }
        )
    );

    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my @outfile_paths       = @{ $io{out}{file_paths} };
    my %outfile_path        = %{ $io{out}{file_path_href} };

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle      = IO::Handle->new();
    my $xargsfilehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $core_number,
            directory_id          => $case_id,
            filehandle            => $filehandle,
            job_id_href           => $job_id_href,
            memory_allocation     => $recipe{memory},
            process_time          => $recipe{time},
            recipe_directory      => $recipe_name,
            recipe_name           => $recipe_name,
            temp_directory        => $temp_directory,
        }
    );

    ### SHELL:

    say {$filehandle} q{## bcftools rho calculation};

    my $xargs_file_path_prefix;

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            core_number        => $core_number,
            filehandle         => $filehandle,
            file_path          => $recipe_file_path,
            recipe_info_path   => $recipe_info_path,
            xargsfilehandle    => $xargsfilehandle,
            xargs_file_counter => $xargs_file_counter,
        }
    );

  CONTIG:
    foreach my $contig (@contigs_size_ordered) {

        ## Get parameters
        my @sample_ids;
        if ( defined $parameter_href->{cache}{affected}
            && @{ $parameter_href->{cache}{affected} } )
        {

            push @sample_ids, $parameter_href->{cache}{affected}[0];
        }
        else {
            # No affected - pick any sample_id
            push @sample_ids, $active_parameter_href->{sample_ids}[0];
        }

        my $roh_outfile_path = $outfile_path_prefix . $DOT . $contig . $DOT . q{roh};
        bcftools_roh(
            {
                af_file_path => $active_parameter_href->{rhocall_frequency_file},
                filehandle   => $xargsfilehandle,
                infile_path  => $infile_path{$contig},
                outfile_path => $roh_outfile_path,
                samples_ref  => \@sample_ids,
                skip_indels  => 1,    # Skip indels as their genotypes are enriched for errors
            }
        );
        print {$xargsfilehandle} $SEMICOLON . $SPACE;

        rhocall_annotate(
            {
                filehandle   => $xargsfilehandle,
                infile_path  => $infile_path{$contig},
                outfile_path => $outfile_path{$contig},
                rohfile_path => $roh_outfile_path,
                q{v14}       => 1,

            }
        );
        say {$xargsfilehandle} $NEWLINE;
    }

    close $filehandle or $log->logcroak(q{Could not close filehandle});
    close $xargsfilehandle
      or $log->logcroak(q{Could not close xargsfilehandle});

    if ( $recipe{mode} == 1 ) {

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                path             => $outfile_paths[0],
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
            }
        );
        submit_recipe(
            {
                base_command                      => $profile_base_command,
                case_id                           => $case_id,
                dependency_method                 => q{sample_to_case},
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

sub analysis_rhocall_viz {

## Function : Detect runs of homo/autozygosity and generate bed file for chromograph
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Recipe name
##          : $sample_id               => Sample id
##          : $sample_info_href        => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_id;
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

    use MIP::File::Path qw{ remove_file_path_suffix };
    use MIP::File_info qw{ get_io_files parse_io_outfiles };
    use MIP::Program::Bcftools qw{ bcftools_index bcftools_roh bcftools_view };
    use MIP::Program::Gnu::Coreutils qw{ gnu_mv };
    use MIP::Program::Picardtools qw{ picardtools_updatevcfsequencedictionary };
    use MIP::Program::Rhocall qw{ rhocall_viz };
    use MIP::Program::Ucsc qw{ ucsc_wig_to_big_wig };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Reference qw{ write_contigs_size_file };
    use MIP::Sample_info
      qw{ set_file_path_to_store set_recipe_metafile_in_sample_info set_recipe_outfile_in_sample_info };
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
            recipe_name    => q{variant_annotation},
            stream         => q{out},
        }
    );
    my $infile_name_prefix = $io{out}{file_name_prefix};
    my $infile_path_prefix = $io{out}{file_path_prefix};
    my $infile_path        = $infile_path_prefix . q{.vcf.gz};

    my %recipe = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id         => $recipe{job_id_chain},
                id               => $sample_id,
                file_info_href   => $file_info_href,
                file_name_prefix => $infile_name_prefix,
                iterators_ref    => [$sample_id],
                outdata_dir      => $active_parameter_href->{outdata_dir},
                parameter_href   => $parameter_href,
                recipe_name      => $recipe_name,
            }
        )
    );
    my $outdir_path    = $io{out}{dir_path};
    my $outfile_path   = $io{out}{file_path};
    my $outfile_suffix = $io{out}{file_suffix};

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $recipe{core_number},
            directory_id          => $sample_id,
            filehandle            => $filehandle,
            job_id_href           => $job_id_href,
            memory_allocation     => $recipe{memory},
            process_time          => $recipe{time},
            recipe_directory      => $recipe_name,
            recipe_name           => $recipe_name,
        }
    );

    ### SHELL:

    say {$filehandle} q{## } . $recipe_name;

    my $sample_outfile_path_prefix = remove_file_path_suffix(
        {
            file_path         => $outfile_path,
            file_suffixes_ref => [$outfile_suffix],
        }
    );
    my $sample_vcf = $sample_outfile_path_prefix . q{.vcf.gz};
    bcftools_view(
        {
            filehandle   => $filehandle,
            infile_path  => $infile_path,
            min_ac       => 1,
            outfile_path => $sample_vcf,
            output_type  => q{z},
            samples_ref  => [$sample_id],
        }
    );
    say {$filehandle} $NEWLINE;

    bcftools_index(
        {
            filehandle  => $filehandle,
            infile_path => $sample_vcf,
            output_type => q{tbi},
        }
    );
    say {$filehandle} $NEWLINE;

    bcftools_roh(
        {
            af_tag       => q{GNOMADAF},
            filehandle   => $filehandle,
            infile_path  => $sample_vcf,
            outfile_path => $sample_outfile_path_prefix . q{.roh},
            skip_indels  => 1,    # Skip indels as their genotypes are enriched for errors
        }
    );
    say {$filehandle} $NEWLINE;

    picardtools_updatevcfsequencedictionary(
        {
            filehandle   => $filehandle,
            infile_path  => $sample_vcf,
            java_jar     => catfile( $active_parameter_href->{picardtools_path}, q{picard.jar} ),
            outfile_path => $sample_outfile_path_prefix . q{.vcf},
            sequence_dictionary => $active_parameter_href->{human_genome_reference},
        }
    );
    say {$filehandle} $NEWLINE;

    rhocall_viz(
        {
            af_tag       => q{GNOMADAF},
            filehandle   => $filehandle,
            infile_path  => $sample_outfile_path_prefix . q{.vcf},
            outdir_path  => $outdir_path,
            rohfile_path => $sample_outfile_path_prefix . q{.roh},
            wig          => 1,
        }
    );
    say {$filehandle} $NEWLINE;

    ## Rename files
  FILE_SUFFIX:
    foreach my $file_suffix (qw{ .bed .wig }) {

        gnu_mv(
            {
                filehandle   => $filehandle,
                infile_path  => catfile( $outdir_path, q{output} . $file_suffix ),
                outfile_path => $sample_outfile_path_prefix . $file_suffix,
            }
        );
        print {$filehandle} $NEWLINE;
    }
    print {$filehandle} $NEWLINE;

    ## Create chromosome name and size file
    my $contigs_size_file_path = catfile( $outdir_path, q{contigs_size_file} . $DOT . q{tsv} );
    write_contigs_size_file(
        {
            fai_file_path => $active_parameter_href->{human_genome_reference} . $DOT . q{fai},
            outfile_path  => $contigs_size_file_path,
        }
    );

    say {$filehandle} q{## Create wig index files};
    ucsc_wig_to_big_wig(
        {
            clip                   => 1,
            contigs_size_file_path => $contigs_size_file_path,
            filehandle             => $filehandle,
            infile_path            => $sample_outfile_path_prefix . $DOT . q{wig},
            outfile_path           => $outfile_path,
        }
    );
    say {$filehandle} $NEWLINE;

    ## Close filehandle
    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe{mode} == 1 ) {

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                infile           => $infile_path,
                path             => $outfile_path,
                recipe_name      => $recipe_name,
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        set_file_path_to_store(
            {
                format           => q{bw},
                id               => $sample_id,
                path             => $outfile_path,
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
            }
        );

        submit_recipe(
            {
                base_command                      => $profile_base_command,
                case_id                           => $case_id,
                dependency_method                 => q{case_to_sample},
                job_id_chain                      => $recipe{job_id_chain},
                job_id_href                       => $job_id_href,
                job_reservation_name              => $active_parameter_href->{job_reservation_name},
                log                               => $log,
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
