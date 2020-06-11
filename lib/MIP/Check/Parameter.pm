package MIP::Check::Parameter;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie;
use Email::Valid;
use Readonly;
use List::MoreUtils qw { any uniq };

## MIPs lib/
use MIP::Constants qw{ $COMMA $DOLLAR_SIGN $DOT $LOG_NAME $NEWLINE $SINGLE_QUOTE $SPACE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.47;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      check_active_installation_parameters
      check_recipe_fastq_compatibility
    };
}

sub check_active_installation_parameters {

## Function : Some active_parameter checks that are common to both installations. Returns "1" if all is OK
## Returns  : 1 or exit
## Arguments: $project_id => Project id
##          : sbatch_mode => Sbatch mode boolean

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $project_id;
    my $sbatch_mode;

    my $tmpl = {
        project_id => {
            store       => \$project_id,
            strict_type => 1,
        },
        sbatch_mode => {
            store       => \$sbatch_mode,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Check that a project id has been set if SBATCH mode
    if ( $sbatch_mode
        and not $project_id )
    {
        $log->fatal(
q{The parameter "project_id" must be set when a sbatch installation has been requested}
        );
        exit 1;
    }
    return 1;
}

sub check_recipe_fastq_compatibility {

## Function : Check that the recipe is compatible with the fastq sequence modes. Turn of downstream applications otherwise
## Returns  :
## Arguments: $active_parameter_href   => Active parameter hash {REF}
##          : $file_info_href          => File info hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $recipe_name             => Recipe name

    my ($arg_href) = @_;

    ## Flatten arguments
    my $active_parameter_href;
    my $file_info_href;
    my $parameter_href;
    my $recipe_name;

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
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        recipe_name => {
            required    => 1,
            defined     => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Active_parameter qw{ set_recipe_mode };
    use MIP::Dependency_tree
      qw{ get_recipe_dependency_tree_chain get_recipes_for_dependency_tree_chain };
    use MIP::File_info qw{ get_consensus_sequence_run_type };

    ## Check if program is going to run
    return if ( $active_parameter_href->{$recipe_name} == 0 );

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my $is_compatible = get_consensus_sequence_run_type(
        {
            file_info_href => $file_info_href,
            sample_ids_ref => $active_parameter_href->{sample_ids},
        }
    );

    if ( not $is_compatible ) {

        ## Turn of current recipe and downstream recipes
        $log->warn(q{Multiple sequence run types detected});
        $log->warn(qq{Turning off $recipe_name and downstream recipes});

        my $recipe_chain;
        get_recipe_dependency_tree_chain(
            {
                recipe               => $recipe_name,
                dependency_tree_href => $parameter_href->{dependency_tree_href},
                chain_id_ref         => \$recipe_chain,
            }
        );

        my @chain_recipes = get_recipes_for_dependency_tree_chain(
            {
                dependency_tree_href    => $parameter_href->{dependency_tree_href},
                chain_initiation_point  => $recipe_chain,
                recipe_initiation_point => $recipe_name,
            }
        );

        ## Turn of recipes
        set_recipe_mode(
            {
                active_parameter_href => $active_parameter_href,
                recipes_ref           => \@chain_recipes,
                mode                  => 0,
            }
        );
    }

    return $is_compatible;
}

1;
