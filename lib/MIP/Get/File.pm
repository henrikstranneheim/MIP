package MIP::Get::File;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use Readonly;
use List::MoreUtils qw{ any };

## MIPs lib/
use MIP::Constants qw{ $COMMA $DOT $LOG_NAME $NEWLINE $PIPE $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.15;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      get_exome_target_bed_file
      get_io_files
      get_merged_infile_prefix
      get_path_entries
    };
}

sub get_exome_target_bed_file {

## Function : Get exome_target_bed file for specfic sample_id and add file_ending from file_info hash if supplied
## Returns  : $exome_target_bed_file
## Arguments: $exome_target_bed_href => Exome target bed files lnked to sample ids
##          : $file_ending           => File ending to add to file
##          : $sample_id             => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $exome_target_bed_href;
    my $file_ending;
    my $sample_id;

    my $tmpl = {
        exome_target_bed_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$exome_target_bed_href,
        },
        file_ending => { store => \$file_ending },
        sample_id   => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sample_id
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## To track sample_ids with capture kits
    my %seen;

    ## Make local hash copy to keep in scope
    my %exome_target_bed = %{$exome_target_bed_href};

  BED_FILE:
    while ( my ( $exome_target_bed_file, $sample_id_string ) = each %exome_target_bed ) {

        my @capture_kit_samples = split $COMMA, $sample_id_string;

        ## Count number of times sample_id has been seen
      SAMPLE:
        foreach my $samples (@capture_kit_samples) {

            $seen{$samples}++;
        }

        ## If capture_kit sample_id is associated with exome_target_bed file
        if ( any { $_ eq $sample_id } @capture_kit_samples ) {

            if ( defined $file_ending ) {

                $exome_target_bed_file .= $file_ending;
            }
            return $exome_target_bed_file;
        }
    }
    ## Found exome target bed file for sample
    return if ( $seen{$sample_id} );

    $log->fatal(
        q{Could not detect }
          . $sample_id
          . q{ in '-exome_target_bed' associated files in sub routine get_exome_target_bed_file},
        $NEWLINE
    );
    exit 1;
}

sub get_io_files {

## Function : Get the io files per chain, id and stream
## Returns  : %io
## Arguments: $id             => Id (sample or case)
##          : $file_info_href => File info hash {REF}
##          : $parameter_href => Parameter hash {REF}
##          : $recipe_name    => Recipe name
##          : $stream         => Stream (in or out or temp)
##          : $temp_directory => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $id;
    my $file_info_href;
    my $parameter_href;
    my $recipe_name;
    my $stream;
    my $temp_directory;

    my $tmpl = {
        id => {
            defined     => 1,
            required    => 1,
            store       => \$id,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
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
        stream => {
            allow       => [qw{ in temp out }],
            defined     => 1,
            required    => 1,
            store       => \$stream,
            strict_type => 1,
        },
        temp_directory => {
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Avoid autovivification of variable
    use Data::Diver qw{ Dive };
    use List::MoreUtils qw{ before };
    use MIP::Set::File qw{ set_io_files };

    ## Constants
    Readonly my $CHAIN_MAIN => q{CHAIN_MAIN};

    ## Unpack
    my $chain_id = $parameter_href->{$recipe_name}{chain};

    ## Not first in chain - return file features
    if ( Dive( $file_info_href, ( q{io}, $chain_id, $id, $recipe_name, $stream ) ) ) {

        return %{ $file_info_href->{io}{$chain_id}{$id}{$recipe_name} };
    }
    else {
        ## First in chain - need to find out stream file features of
        ## correct upstream recipe

        my $upstream_direction = q{out};

        ## Unpack
        my @order_recipes =
          @{ $parameter_href->{cache}{order_recipes_ref} };

        ## Find upstream recipes starting from (and not including) recipe_name
        my @upstream_recipes =
          reverse before { $_ eq $recipe_name } @order_recipes;

      UPSTREAM_RECIPE:
        foreach my $upstream_recipe (@upstream_recipes) {

            # Get chain id
            my $upstream_chain_id = $parameter_href->{$upstream_recipe}{chain};

            ## No io file features found in chain and stream
            next UPSTREAM_RECIPE
              if (
                not Dive(
                    $file_info_href,
                    (
                        q{io},            $upstream_chain_id, $id,
                        $upstream_recipe, $upstream_direction
                    )
                )
              );

            ## PARALLEL CHAIN with multiple recipes
            # second in chain
            if ( $upstream_chain_id eq $chain_id ) {

                ## Switch upstream out to recipe in - i.e. inherit from upstream
                _inherit_upstream_io_files(
                    {
                        chain_id           => $chain_id,
                        id                 => $id,
                        file_info_href     => $file_info_href,
                        recipe_name        => $recipe_name,
                        stream             => $stream,
                        temp_directory     => $temp_directory,
                        upstream_direction => $upstream_direction,
                        upstream_chain_id  => $upstream_chain_id,
                        upstream_recipe    => $upstream_recipe,
                    }
                );

                ##  Return set file features
                return %{ $file_info_href->{io}{$chain_id}{$id}{$recipe_name} };
            }

            ## Do not inherit from other chains than self or MAIN
            next UPSTREAM_RECIPE if ( $upstream_chain_id ne q{MAIN} );

            ## Found io file features found in chain, id, recipe and stream
            if (
                Dive(
                    $file_info_href,
                    (
                        q{io},            $upstream_chain_id, $id,
                        $upstream_recipe, $upstream_direction
                    )
                )
              )
            {

                ## Switch upstream out to recipe in - i.e. inherit from upstream
                _inherit_upstream_io_files(
                    {
                        chain_id           => $chain_id,
                        id                 => $id,
                        file_info_href     => $file_info_href,
                        recipe_name        => $recipe_name,
                        stream             => $stream,
                        temp_directory     => $temp_directory,
                        upstream_direction => $upstream_direction,
                        upstream_chain_id  => $upstream_chain_id,
                        upstream_recipe    => $upstream_recipe,
                    }
                );

                ##  Return set file features
                return %{ $file_info_href->{io}{$chain_id}{$id}{$recipe_name} };
            }
        }
    }

    ## At root of initation map - add base
    # Build infiles path
    my @base_file_paths =
      map { catfile( $file_info_href->{$id}{mip_infiles_dir}, $_ ) }
      @{ $file_info_href->{$id}{mip_infiles} };

    set_io_files(
        {
            chain_id       => $CHAIN_MAIN,
            id             => $id,
            file_paths_ref => \@base_file_paths,
            file_info_href => $file_info_href,
            recipe_name    => $recipe_name,
            stream         => $stream,
            temp_directory => $temp_directory,
        }
    );
    return %{ $file_info_href->{io}{$CHAIN_MAIN}{$id}{$recipe_name} };
}

sub get_merged_infile_prefix {

## Function : Get the merged infile prefix for sample id
## Returns  : $merged_infile_prefix
## Arguments: $file_info_href => File info hash {REF}
##          : $sample_id      => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_info_href;
    my $sample_id;

    my $tmpl = {
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href
        },
        sample_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sample_id
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return $file_info_href->{$sample_id}{merged_infile};
}

sub get_path_entries {

## Function  : Collects all recipes outfile path(s) created by MIP as Path->value located in %sample_info.
## Returns   :
## Arguments : $paths_ref        => Holds the collected paths {REF}
##           : $sample_info_href => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $paths_ref;
    my $sample_info_href;

    my $tmpl = {
        paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$paths_ref,
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

    ## Copy hash to enable recursive removal of keys
    my %info = %{$sample_info_href};

    ## Temporary array for collecting outdirectories within the same recipe
    my @outdirectories;

    ## Temporary array for collecting outfile within the same recipe
    my @outfiles;

  KEY_VALUE_PAIR:
    while ( my ( $key, $value ) = each %info ) {

        if ( ref $value eq q{HASH} ) {

            get_path_entries(
                {
                    paths_ref        => $paths_ref,
                    sample_info_href => $value,
                }
            );
        }
        else {

            ## Required for first dry-run
            next KEY_VALUE_PAIR if ( not $value );

            ## Check if key is "path" and adds value to @paths_ref if true.
            _check_and_add_to_array(
                {
                    key       => $key,
                    paths_ref => $paths_ref,
                    value     => $value,
                }
            );

            ## Check if key is "outdirectory" or "outfile"  and adds joined value to @paths_ref if true.
            _collect_outfile(
                {
                    key                => $key,
                    paths_ref          => $paths_ref,
                    outdirectories_ref => \@outdirectories,
                    outfiles_ref       => \@outfiles,
                    value              => $value,
                }
            );

            delete $info{$value};
        }
    }
    return;
}

sub _check_and_add_to_array {

## Function  : Check if Key name is "path" and adds to @paths_ref if true.
## Returns   :
## Arguments : $keyName   => Hash key
##           : $paths_ref => Holds the collected paths {REF}
##           : $value     => Hash value

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $key;
    my $paths_ref;
    my $value;

    my $tmpl = {
        key       => { defined => 1, required => 1, store => \$key, strict_type => 1, },
        paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$paths_ref,
            strict_type => 1,
        },
        value => { required => 1, store => \$value, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return if ( $key ne q{path} );

    ## Do not add same path twice
    if ( not any { $_ eq $value } @{$paths_ref} ) {

        push @{$paths_ref}, $value;
    }
    return;
}

sub _collect_outfile {

## Function  : Check if Key name is "outdirectory" or "outfile"  and adds to @paths_ref if true.
## Returns   :
## Arguments : $key                => Hash key
##           : $outdirectories_ref => Holds temporary outdirectory path(s) {Optional, REF}
##           : $outfiles_ref       => Holds temporary outdirectory path(s) {Optional, REF}
##           : $value              => Hash value
##           : $paths_ref          => Holds the collected paths {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $key;
    my $outdirectories_ref;
    my $outfiles_ref;
    my $paths_ref;
    my $value;

    my $tmpl = {
        paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$paths_ref,
            strict_type => 1,
        },
        outdirectories_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$outdirectories_ref,
            strict_type => 1,
        },
        outfiles_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$outfiles_ref,
            strict_type => 1,
        },
        value => { defined => 1, required => 1, store => \$value, },
        key   => { defined => 1, store => \$key, required => 1, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    if ( $key eq q{outdirectory} ) {

        push @{$outdirectories_ref}, $value;
    }
    if ( $key eq q{outfile} ) {

        push @{$outfiles_ref}, $value;
    }

    ## Both outdirectory and outfile have been collected, time to join
    if ( @{$outdirectories_ref} && @{$outfiles_ref} ) {

        my $path = catfile( $outdirectories_ref->[0], $outfiles_ref->[0] );

        ## Do not add same path twice
        if ( not any { $_ eq $path } @{$paths_ref} ) {

            push @{$paths_ref}, catfile( $outdirectories_ref->[0], $outfiles_ref->[0] );

            ## Restart
            @{$outdirectories_ref} = ();
            @{$outfiles_ref}       = ();
        }
    }
    return;
}

sub _inherit_upstream_io_files {

## Function : Switch upstream out to recipe in - i.e. inherit from upstream
## Returns  : %io
## Arguments: $chain_id           => Chain id
##          : $id                 => Id (sample or case)
##          : $file_info_href     => File info hash {REF}
##          : $recipe_name        => Recipe name
##          : $stream             => Stream (in or out or temp)
##          : $temp_directory     => Temporary directory
##          : $upstream_direction => Upstream direction
##          : $upstream_chain_id  => Upstream chain id
##          : $upstream_recipe    => Upstream recipe

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $chain_id;
    my $id;
    my $file_info_href;
    my $recipe_name;
    my $stream;
    my $temp_directory;
    my $upstream_direction;
    my $upstream_chain_id;
    my $upstream_recipe;

    my $tmpl = {
        chain_id => {
            defined     => 1,
            required    => 1,
            store       => \$chain_id,
            strict_type => 1,
        },
        id => {
            defined     => 1,
            required    => 1,
            store       => \$id,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
        stream => {
            allow       => [qw{ in temp out }],
            defined     => 1,
            required    => 1,
            store       => \$stream,
            strict_type => 1,
        },
        temp_directory => {
            store       => \$temp_directory,
            strict_type => 1,
        },
        upstream_direction => {
            defined     => 1,
            required    => 1,
            store       => \$upstream_direction,
            strict_type => 1,
        },
        upstream_chain_id => {
            defined     => 1,
            required    => 1,
            store       => \$upstream_chain_id,
            strict_type => 1,
        },
        upstream_recipe => {
            defined     => 1,
            required    => 1,
            store       => \$upstream_recipe,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Switch upstream out to recipe in - i.e. inherit from upstream
    my @upstream_outfile_paths =
      @{ $file_info_href->{io}{$upstream_chain_id}{$id}
          {$upstream_recipe}{$upstream_direction}{file_paths} };
    set_io_files(
        {
            chain_id       => $chain_id,
            id             => $id,
            file_paths_ref => \@upstream_outfile_paths,
            file_info_href => $file_info_href,
            recipe_name    => $recipe_name,
            stream         => $stream,
            temp_directory => $temp_directory,
        }
    );
    return;
}

1;
