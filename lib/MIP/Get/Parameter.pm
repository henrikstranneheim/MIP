package MIP::Get::Parameter;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use Readonly;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.27;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      get_gatk_intervals
      get_install_parameter_attribute
      get_package_source_env_cmds
      get_program_version
      get_programs_for_shell_installation
      get_recipe_resources
      get_recipe_attributes
    };
}

## Constants
Readonly my $TWO => 2;

sub get_gatk_intervals {

## Function : Generate and return interval hash
## Returns  : %gatk_intervals
## Arguments: $analysis_type         => Analysis type
##          : $contigs_ref           => Contigs to split in file
##          : $exome_target_bed_href => Exome target bed files lnked to sample ids
##          : $file_ending           => File ending to add {Optional}
##          : $filehandle            => Filehandle to write to
##          : $log                   => Log
##          : $max_cores_per_node    => Maximum core per node
##          : $outdirectory          => Outdirectory
##          : $reference_dir         => MIP reference directory
##          : $sample_id             => Sample_id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $analysis_type;
    my $contigs_ref;
    my $exome_target_bed_href;
    my $file_ending;
    my $filehandle;
    my $log;
    my $outdirectory;
    my $reference_dir;
    my $sample_id;

    ## Default(s)
    my $max_cores_per_node;

    my $tmpl = {
        analysis_type => {
            defined     => 1,
            required    => 1,
            store       => \$analysis_type,
            strict_type => 1,
        },
        contigs_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$contigs_ref,
            strict_type => 1,
        },
        exome_target_bed_href => {
            store => \$exome_target_bed_href,
        },
        file_ending => {
            store       => \$file_ending,
            strict_type => 1,
        },
        filehandle => { store => \$filehandle, },
        log        => {
            store => \$log,
        },
        max_cores_per_node => {
            allow       => qr/ ^\d+$ /sxm,
            default     => 1,
            store       => \$max_cores_per_node,
            strict_type => 1,
        },
        outdirectory => {
            store       => \$outdirectory,
            strict_type => 1,
        },
        reference_dir => {
            store       => \$reference_dir,
            strict_type => 1,
        },
        sample_id => {
            store       => \$sample_id,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::File qw{ get_exome_target_bed_file };
    use MIP::File::Interval qw{ generate_contig_interval_file };

    ## Store gatk interval
    my %gatk_intervals;

    if ( $analysis_type eq q{wes} ) {

        my $exome_target_bed_file = get_exome_target_bed_file(
            {
                exome_target_bed_href => $exome_target_bed_href,
                file_ending           => $file_ending,
                sample_id             => $sample_id,
            }
        );

        ## Generate contig specific interval_list and return hash with paths
        %gatk_intervals = generate_contig_interval_file(
            {
                contigs_ref           => $contigs_ref,
                exome_target_bed_file => $exome_target_bed_file,
                filehandle            => $filehandle,
                max_cores_per_node    => $max_cores_per_node,
                outdirectory          => $outdirectory,
                reference_dir         => $reference_dir,
            }
        );
    }
    else {
        ## Key-value pairs are identical for WGS/WTS
        %gatk_intervals = map { $_ => [$_] } @{$contigs_ref};
    }

    return %gatk_intervals;
}

sub get_install_parameter_attribute {

## Function : Return parameter attribute from hash
## Returns  :
## Arguments: $parameter_href => Holds all parameters {REF}
##          : $parameter_name => Name of key to return

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $parameter_name;

    my $tmpl = {
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        parameter_name => {
            store       => \$parameter_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Make sure that the supplied key exists
    croak(qq{Could not find parameter_name key: '$parameter_name' in hash})
      if ( not exists $parameter_href->{$parameter_name} );

    ## Hash attribute
    if ( ref $parameter_href->{$parameter_name} eq q{HASH} ) {

        return %{ $parameter_href->{$parameter_name} };
    }
    ## ARRAY attribute
    if ( ref $parameter_href->{$parameter_name} eq q{ARRAY} ) {

        return @{ $parameter_href->{$parameter_name} };
    }
    ## Scalar attribute
    return $parameter_href->{$parameter_name};
}

sub get_package_source_env_cmds {

## Function : Get package source environment commands
## Returns  : @source_environment_cmds
## Arguments: $active_parameter_href => The active parameters for this analysis hash {REF}
##          : $package_name          => Package name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $package_name;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        package_name => {
            defined     => 1,
            required    => 1,
            store       => \$package_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Active_parameter qw{ get_package_env_attributes };
    use MIP::Environment::Manager qw{ get_env_method_cmds };
    use MIP::Environment::Container qw{ parse_container_bind_paths };

    ## Initilize variable
    my @source_environment_cmds;

    my ( $env_name, $env_method ) = get_package_env_attributes(
        {
            load_env_href => $active_parameter_href->{load_env},
            package_name  => $package_name,
        }
    );

    ## Could not find recipe within env
    if ( not $env_name ) {

        ## Fall back to MIPs MAIN env
        ( $env_name, $env_method ) = get_package_env_attributes(
            {
                load_env_href => $active_parameter_href->{load_env},
                package_name  => q{mip},
            }
        );
    }
    ## Prior to env command special case
    ## for recipes needing addtional processing
    my $prior_to_load_cmd = $active_parameter_href->{load_env}{$env_name}{$package_name};
    if ($prior_to_load_cmd) {

        push @source_environment_cmds, $prior_to_load_cmd;
    }

    ## Append container bind variable to source_environment_cmds_ref
    parse_container_bind_paths(
        {
            active_parameter_href       => $active_parameter_href,
            package_name                => $package_name,
            source_environment_cmds_ref => \@source_environment_cmds,
        }
    );

    ## Get env load command
    my @env_method_cmds = get_env_method_cmds(
        {
            action     => q{load},
            env_name   => $env_name,
            env_method => $env_method,
        }
    );
    push @source_environment_cmds, @env_method_cmds;

    return @source_environment_cmds;
}

sub get_programs_for_shell_installation {

## Function  : Get the programs that are to be installed via SHELL
## Returns   : @shell_programs
## Arguments : $conda_programs_href        => Hash with conda progrmas {REF}
##           : $log                        => Log
##           : $prefer_shell               => Path to conda environment
##           : $shell_install_programs_ref => Array with programs selected for shell installation {REF}
##           : $shell_programs_href        => Hash with shell programs {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_programs_href;
    my $log;
    my $prefer_shell;
    my $shell_install_programs_ref;
    my $shell_programs_href;

    my $tmpl = {
        conda_programs_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$conda_programs_href,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        prefer_shell => {
            allow       => [ undef, 0, 1 ],
            required    => 1,
            store       => \$prefer_shell,
            strict_type => 1,
        },
        shell_install_programs_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$shell_install_programs_ref,
            strict_type => 1,
        },
        shell_programs_href => {
            default  => {},
            required => 1,
            store    => \$shell_programs_href,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use Array::Utils qw{ intersect array_minus unique };

    return if not keys %{$shell_programs_href};

    my @shell_programs = keys %{$shell_programs_href};
    my @conda_programs = keys %{$conda_programs_href};

    if ($prefer_shell) {

        # Only get the selected programs otherwise leave the array unaltered
        if ( @{$shell_install_programs_ref} ) {

            # Get the intersect between the two arrays
            @shell_programs =
              intersect( @shell_programs, @{$shell_install_programs_ref} );
        }
    }
    elsif ( @{$shell_install_programs_ref} ) {

        # Get elements in @shell_programs that are not part of the conda hash
        my @shell_only_programs = array_minus( @shell_programs, @conda_programs );

        # Add the selected program(s) and remove possible duplicates
        @shell_programs = unique( @shell_only_programs, @{$shell_install_programs_ref} );
    }
    else {
        # If no shell preferences only add programs lacking conda counterpart
        @shell_programs = array_minus( @shell_programs, @conda_programs );
    }

    return @shell_programs;
}

sub get_recipe_attributes {

## Function : Return recipe attributes
## Returns  : $attribute | %attribute
## Arguments: $attribute      => Attribute key
##          : $parameter_href => Holds all parameters
##          : $recipe_name    => Recipe name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $attribute;
    my $parameter_href;
    my $recipe_name;

    my $tmpl = {
        attribute => {
            store       => \$attribute,
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

    ## Get attribute value
    if ( defined $attribute && $attribute ) {

        return $parameter_href->{$recipe_name}{$attribute};
    }

    ## Get recipe attribute hash
    return %{ $parameter_href->{$recipe_name} };
}

sub get_recipe_resources {

## Function : Return recipe resources
## Returns  : $recipe_resource | %recipe_resource
## Arguments: $active_parameter_href => The active parameters for this analysis hash {REF}
##          : $recipe_name           => Recipe name
##          : $recipe_resource       => Recipe parameter key

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $recipe_name;
    my $resource;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
        resource => {
            allow       => [qw{ core_number load_env_ref memory time }],
            store       => \$resource,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Environment::Cluster qw{ check_recipe_memory_allocation };

    ## Initilize variable
    my @source_environment_cmds = get_package_source_env_cmds(
        {
            active_parameter_href => $active_parameter_href,
            package_name          => $recipe_name,
        }
    );

    my $core_number    = $active_parameter_href->{recipe_core_number}{$recipe_name};
    my $process_memory = $active_parameter_href->{recipe_memory}{$recipe_name};
    my $memory;

    ## Multiply memory with processes that are to be launched in the recipe
    if ( $process_memory and $core_number ) {
        $memory = $process_memory * $core_number;
    }
    ## Set default recipe memory allocation if it hasn't been specified
    elsif ( not $process_memory and $core_number ) {
        $memory = $core_number * $active_parameter_href->{core_ram_memory};
    }
    elsif ( not $process_memory and not $core_number ) {
        $memory = $active_parameter_href->{core_ram_memory};
    }
    else {
        $memory = $process_memory;
    }

    check_recipe_memory_allocation(
        {
            node_ram_memory          => $active_parameter_href->{node_ram_memory},
            recipe_memory_allocation => $memory,
        }
    );

    my %recipe_resource = (
        core_number  => $core_number,
        load_env_ref => \@source_environment_cmds,
        memory       => $memory,
        time         => $active_parameter_href->{recipe_time}{$recipe_name},
    );

    ## Return specified recipe resource
    if ( defined $resource && $resource ) {
        return $recipe_resource{$resource};
    }

    ## Return recipe resource hash
    return %recipe_resource;

}

1;
