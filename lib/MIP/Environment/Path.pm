package MIP::Environment::Path;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile splitdir };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib
use MIP::Constants qw{ $COLON $COMMA $DOUBLE_QUOTE $EQUALS $LOG_NAME $SEMICOLON $SPACE };

## Constants
Readonly my $MINUS_TWO => -2;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      get_conda_path
      reduce_dir_paths
    };
}

sub get_conda_path {

## Function : Get path to conda directory
## Returns  : $conda_path
## Arguments: $bin_file => Bin file to test

    my ($arg_href) = @_;

    ## Default(s)
    my $bin_file;

    my $tmpl = {
        bin_file => {
            default     => q{conda},
            store       => \$bin_file,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use IPC::Cmd qw{ can_run };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Find path to conda bin
    my $conda_path = can_run($bin_file);

    if ( not $conda_path ) {

        $log->fatal( q{Failed to find path to: } . $bin_file );
        exit 1;
    }

    ## Split dirs to array
    my @conda_path_dirs = File::Spec->splitdir($conda_path);

    ## Traverse to conda dir from binary
    splice @conda_path_dirs, $MINUS_TWO;

    ## Return path to conda main directory
    return catdir(@conda_path_dirs);
}

sub reduce_dir_paths {

## Function : Parses directory paths and reduces them to a non-overlapping array. No check for existing files or directories
## Returns  : @reduced_dir_paths
## Arguments: $dir_paths_ref => Directory paths to parse {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $dir_paths_ref;

    my $tmpl = {
        dir_paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$dir_paths_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @dir_paths;

    ## Split to dir path to array
  DIR_PATH:
    foreach my $dir_path ( @{$dir_paths_ref} ) {

        next DIR_PATH if ( not defined $dir_path );

        push @dir_paths, [ splitdir($dir_path) ];
    }

    ## Sort according to size
    @dir_paths = sort { @{$a} <=> @{$b} } @dir_paths;

    ## Reformat to strings
    @dir_paths = map { catdir( @{$_} ) } @dir_paths;

    my @reduced_dir_paths;

  BIND_PATH:
    while (@dir_paths) {

        ## Shift array
        my $dir_path = shift @dir_paths;

        ## Save path
        push @reduced_dir_paths, $dir_path;

        ## Get indexes of all the paths in the array that have an identical beginning to one we are testing
        ## The \Q and \E in the regex turns of interpolation
        my @match_idxs =
          grep { $dir_paths[$_] =~ / ^\Q$dir_path\E.* /xms } 0 .. $#dir_paths;

      MATCH_IDX:
        foreach my $match_idx ( reverse @match_idxs ) {

            ## Remove those paths with matching starts
            splice @dir_paths, $match_idx, 1;
        }
    }

    return @reduced_dir_paths;
}

1;
