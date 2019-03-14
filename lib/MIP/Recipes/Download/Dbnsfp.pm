package MIP::Recipes::Download::Dbnsfp;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants
  qw{ $ASTERISK $BACKWARD_SLASH $DASH $NEWLINE $PIPE $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ download_dbnsfp };

}

## Constants
Readonly my $GRCH37_CHR_POS    => 8;
Readonly my $GRCH37_REGION_POS => 9;

sub download_dbnsfp {

## Function : Download dbnsfp
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

    use Cwd;
    use MIP::Get::Parameter qw{ get_recipe_parameters };
    use MIP::Gnu::Coreutils qw{ gnu_cat gnu_head gnu_sort };
    use MIP::Gnu::Software::Gnu_grep qw{ gnu_grep };
    use MIP::Language::Awk qw{ awk };
    use MIP::Program::Utility::Htslib qw{ htslib_bgzip htslib_tabix };
    use MIP::Recipes::Download::Get_reference qw{ get_reference };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_no_dependency_dead_end };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger( uc q{mip_download} );

    ## Unpack parameters
    my $reference_dir = $active_parameter_href->{reference_dir};
    my @reference_genome_versions =
      @{ $active_parameter_href->{reference_genome_versions} };

    my ( $core_number, $time, @source_environment_cmds ) = get_recipe_parameters(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    ## Set recipe mode
    my $recipe_mode = $active_parameter_href->{$recipe_name};

    ## Filehandle(s)
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $core_number,
            directory_id          => q{mip_download},
            FILEHANDLE            => $FILEHANDLE,
            job_id_href           => $job_id_href,
            log                   => $log,
            outdata_dir           => $reference_dir,
            outscript_dir         => $reference_dir,
            process_time          => $time,
            recipe_directory      => $recipe_name . $UNDERSCORE . $reference_version,
            recipe_name           => $recipe_name,
            temp_directory        => $temp_directory,
            source_environment_commands_ref => \@source_environment_cmds,
        }
    );

    ### SHELL:

    say {$FILEHANDLE} q{## } . $recipe_name;

    get_reference(
        {
            FILEHANDLE     => $FILEHANDLE,
            recipe_name    => $recipe_name,
            reference_dir  => $reference_dir,
            reference_href => $reference_href,
            quiet          => $quiet,
            verbose        => $verbose,
        }
    );

    say {$FILEHANDLE} q{## Build dbnsfp header file after unzip};

    ## Build dbnsfp chr file name after unzip
    my $dbnsfp_chr_file_name =
      q{dbNSFP} . $reference_version . $UNDERSCORE . q{variant.chr} . $ASTERISK;
    my $dbnsfp_chr1_file_name =
      q{dbNSFP} . $reference_version . $UNDERSCORE . q{variant.chr1};

    my $dbnsfp_chr_file_path  = catfile( $reference_dir, $dbnsfp_chr_file_name );
    my $dbnsfp_chr1_file_path = catfile( $reference_dir, $dbnsfp_chr1_file_name );
    my $header_file_path      = catfile( $reference_dir, q{dbnsfp_header.txt} );
    my $reformated_outfile    = join $UNDERSCORE,
      (
        $genome_version, $recipe_name, q{reformated}, q{-} . $reference_version . q{-.txt}
      );
    my $reformated_outfile_path = catfile( $reference_dir, $reformated_outfile );
    gnu_head(
        {
            FILEHANDLE      => $FILEHANDLE,
            infile_path     => $dbnsfp_chr1_file_path,
            lines           => 1,
            stdoutfile_path => $header_file_path,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    say {$FILEHANDLE} q{## Switch columns to use with genome build prior to version 38};
    ## Read files
    gnu_cat(
        {
            FILEHANDLE       => $FILEHANDLE,
            infile_paths_ref => [$dbnsfp_chr_file_path],
        }
    );
    say {$FILEHANDLE} $PIPE . $SPACE . $BACKWARD_SLASH;

    ## Skip header line starting with "#chr"
    gnu_grep(
        {
            FILEHANDLE   => $FILEHANDLE,
            invert_match => 1,
            pattern      => q{^#chr},
        }
    );
    say {$FILEHANDLE} $PIPE . $SPACE . $BACKWARD_SLASH;

    ## If not DOT in column 8
    awk(
        {
            FILEHANDLE => $FILEHANDLE,
            statement  => q{$8 != "."},
        }
    );
    say {$FILEHANDLE} $PIPE . $SPACE . $BACKWARD_SLASH;

## Sort nummerical on column 8 and 9
    gnu_sort(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $DASH,
            keys_ref    => [ q{8,8}, q{9,9n} ],
        }
    );
    say {$FILEHANDLE} $PIPE . $SPACE . $BACKWARD_SLASH;

    ## Concatenate header and input stream (DASH)
    gnu_cat(
        {
            FILEHANDLE       => $FILEHANDLE,
            infile_paths_ref => [ $header_file_path, $DASH ],
            stdoutfile_path  => $reformated_outfile_path,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    say {$FILEHANDLE} q{## Compress and index file};
## Compress file
    htslib_bgzip(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $reformated_outfile_path,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

## Index file using tabix
    htslib_tabix(
        {
            begin       => $GRCH37_REGION_POS,
            end         => $GRCH37_REGION_POS,
            FILEHANDLE  => $FILEHANDLE,
            force       => 1,
            infile_path => $reformated_outfile_path . q{.gz},
            sequence    => $GRCH37_CHR_POS,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Close FILEHANDLES
    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});

    if ( $recipe_mode == 1 ) {

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
