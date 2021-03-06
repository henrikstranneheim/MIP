package MIP::Recipes::Download::Gnomad_pli_per_gene;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname fileparse };
use File::Spec::Functions qw{ catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $COMMA $DASH $NEWLINE $PIPE $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ download_gnomad_pli_per_gene };

}

sub download_gnomad_pli_per_gene {

## Function : Download Gnomad pLI scores per gene
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this download hash {REF}
##          : $genome_version        => Human genome version
##          : $job_id_href           => The job_id hash {REF}
##          : $profile_base_command  => Submission profile base command
##          : $recipe_name           => Recipe name
##          : $reference_href        => Reference hash {REF}
##          : $reference_version     => Reference version
##          : $quiet                 => Quiet (no output)
##          : $temp_directory        => Temporary directory for recipe
##          : $verbose               => Verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $genome_version;
    my $job_id_href;
    my $recipe_name;
    my $reference_href;
    my $reference_version;

    ## Default(s)
    my $profile_base_command;
    my $quiet;
    my $temp_directory;
    my $verbose;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        genome_version => {
            store       => \$genome_version,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
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
        reference_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$reference_href,
            strict_type => 1,
        },
        reference_version => {
            defined     => 1,
            required    => 1,
            store       => \$reference_version,
            strict_type => 1,
        },
        quiet => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$quiet,
            strict_type => 1,
        },
        temp_directory => {
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Program::Gnu::Coreutils qw{ gnu_cat gnu_cut };
    use MIP::Program::Gnu::Software::Gnu_grep qw{ gnu_grep };
    use MIP::Processmanagement::Slurm_processes qw{ slurm_submit_job_no_dependency_dead_end };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Recipes::Download::Get_reference qw{ get_reference };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Constants
    Readonly my $HGNC_SYMBOL_COL_NR => 1;
    Readonly my $PLI_COL_NR         => 21;

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger( uc q{mip_download} );

    ## Unpack parameters
    my $reference_dir = $active_parameter_href->{reference_dir};

    my %recipe = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

## Filehandle(s)
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $recipe{core_number},
            directory_id                    => q{mip_download},
            filehandle                      => $filehandle,
            info_file_id                    => $genome_version . $UNDERSCORE . $reference_version,
            job_id_href                     => $job_id_href,
            memory_allocation               => $recipe{memory},
            outdata_dir                     => $reference_dir,
            outscript_dir                   => $reference_dir,
            process_time                    => $recipe{time},
            recipe_data_directory_path      => $active_parameter_href->{reference_dir},
            recipe_directory                => $recipe_name . $UNDERSCORE . $reference_version,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => $recipe{load_env_ref},
        }
    );

    ### SHELL:

    say {$filehandle} q{## } . $recipe_name;

    get_reference(
        {
            filehandle     => $filehandle,
            recipe_name    => $recipe_name,
            reference_dir  => $reference_dir,
            reference_href => $reference_href,
            quiet          => $quiet,
            verbose        => $verbose,
        }
    );

    say {$filehandle} q{## Get HGNC symbol and pLI score};

    ## Remove ".gz" from outfile to create infile
    my $gnu_cat_infile = fileparse( $reference_href->{outfile}, qr/[.]gz/xsm );

    my $gnu_cat_infile_path = catfile( $reference_dir, $gnu_cat_infile );

    gnu_cat(
        {
            filehandle       => $filehandle,
            infile_paths_ref => [ $gnu_cat_infile_path, ],
        }
    );

    print {$filehandle} $SPACE . $PIPE . $SPACE;

    my $list_fields = join $COMMA, ( $HGNC_SYMBOL_COL_NR, $PLI_COL_NR );
    gnu_cut(
        {
            filehandle  => $filehandle,
            infile_path => $DASH,
            list        => $list_fields,
        }
    );

    print {$filehandle} $SPACE . $PIPE . $SPACE;

    my $reformated_outfile = join $UNDERSCORE,
      ( $recipe_name, $DASH, $reference_version . $DASH . q{.txt} );
    my $reformated_outfile_path = catfile( $reference_dir, $reformated_outfile );

    gnu_grep(
        {
            filehandle      => $filehandle,
            invert_match    => 1,
            pattern         => q{NA},
            stdoutfile_path => $reformated_outfile_path,
            word_regexp     => 1,
        }
    );

    ## Close filehandleS
    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe{mode} == 1 ) {

        ## No upstream or downstream dependencies
        slurm_submit_job_no_dependency_dead_end(
            {
                base_command     => $profile_base_command,
                job_id_href      => $job_id_href,
                log              => $log,
                sbatch_file_name => $recipe_file_path,
            }
        );
    }
    return 1;
}

1;
