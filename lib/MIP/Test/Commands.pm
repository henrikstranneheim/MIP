package MIP::Test::Commands;

#### Copyright 2017 Henrik Stranneheim

use strict;
use warnings;
use warnings qw(FATAL utf8);
use utf8;    #Allow unicode characters in this script
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );
use Carp;
use English qw(-no_match_vars);
use autodie;
use Test::More;
use List::Util qw(any);

BEGIN {

    use base qw(Exporter);
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw(test_function);
}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;    #Do not convert to lower case

sub test_function {

##test_function

##Function : Test module function by generating arguments and testing output
##Returns  : "@commands"
##Arguments: $argument_href, $required_arguments_href, $module_function_cref, $function_base_command
##         : $argument_href           => Parameters to submit to module method
##         : $required_arguments_href => Required arguments
##         : $module_function_cref    => Module method to test
##         : $function_base_command   => Function base command

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $argument_href;
    my $required_arguments_href;
    my $function_base_command;
    my $module_function_cref;

    my $tmpl = {
        argument_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$argument_href
        },
        required_arguments_href => {
            default     => {},
            strict_type => 1,
            store       => \$required_arguments_href
        },
        function_base_command => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$function_base_command
        },
        module_function_cref =>
          { required => 1, defined => 1, store => \$module_function_cref },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw(Could not parse arguments!);

  ARGUMENT:
    foreach my $argument ( keys %{$argument_href} ) {

        ### Parameter to test in this loop check if array ref or scalar
        my $input_value;
        my $input_values_ref;

        ## SCALAR
        if ( exists $argument_href->{$argument}{input} ) {

            $input_value = $argument_href->{$argument}{input};
        }
        ## ARRAY
        elsif ( exists $argument_href->{$argument}{inputs_ref} ) {

            $input_values_ref = $argument_href->{$argument}{inputs_ref};
        }

        ## Store commands from module function
        my @commands;

        ## Some functions have mandatory arguments
        if ( %{$required_arguments_href} ) {

            my @args;
            if ($input_values_ref) {

                @args = _build_call(
                    {
                        required_arguments_href => $required_arguments_href,
                        argument                => $argument,
                        input_values_ref        => $input_values_ref,
                    }
                );
            }
            else {

                @args = _build_call(
                    {
                        required_arguments_href => $required_arguments_href,
                        argument                => $argument,
                        input_value             => $input_value,
                    }
                );
            }

            ## Special case for test of FILEHANDLE. Does not return @commands
            if ( $argument eq 'FILEHANDLE' ) {

                _test_write_to_file(
                    {
                        args_ref             => \@args,
                        module_function_cref => $module_function_cref,
                        base_command         => $function_base_command,
                    }
                );
            }
            else {

                ## Submit arguments to coderef sub
                @commands = $module_function_cref->( {@args} );
            }
        }
        else {

            ## Special case for test of FILEHANDLE. Does not return @commands
            if ( $argument eq 'FILEHANDLE' ) {

                _test_write_to_file(
                    {
                        args_ref => [ $argument, $input_value ],
                        module_function_cref => $module_function_cref,
                        base_command         => $function_base_command,
                    }
                );
            }
            else {

                ## Array
                if ($input_values_ref) {

                    @commands =
                      $module_function_cref->(
                        { $argument => $input_values_ref, } );
                }
                else {

                    ## Submit arguments to coderef sub
                    @commands =
                      $module_function_cref->( { $argument => $input_value, } );
                }
            }
        }

        ## Expected return value from sub call
        my $expected_return = $argument_href->{$argument}{expected_output};

        ### Perform tests

        if (@commands) {

            ## Test function_base_command
            _test_base_command(
                {
                    base_command          => $commands[0],
                    expected_base_command => $function_base_command,
                }
            );

            ## Test argument
            ok( ( any { $_ eq $expected_return } @commands ),
                'Argument: ' . $argument );
        }
    }
    return;
}

sub _build_call {

##_build_call

##Function : Build arguments to function
##Returns  : "@arguments"
##Arguments: $required_arguments_href, argument, input_value
##         : $required_arguments_href => Required arguments
##         : $argument                => Argument key to test
##         : $input_value             => Argument value to test

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $required_arguments_href;
    my $input_values_ref;
    my $argument;
    my $input_value;

    my $tmpl = {
        required_arguments_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$required_arguments_href
        },
        input_values_ref => {
            default     => [],
            strict_type => 1,
            store       => \$input_values_ref
        },
        argument => {
            strict_type => 1,
            store       => \$argument
        },
        input_value => {
            strict_type => 1,
            store       => \$input_value
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw(Could not parse arguments!);

    ## Collect required keys and values to generate args
    my @keys = keys %{$required_arguments_href};
    my @values =
      map { ( $required_arguments_href->{$_}{input} ) } @keys;

    ### Combine the specific and required argument keys and values to test
    ## SCALAR
    if ( $argument && $input_value ) {
        push @keys,   $argument;
        push @values, $input_value;
    }
    ## ARRAY
    elsif ( $argument && @{$input_values_ref} ) {
        push @keys,   $argument;
        push @values, $input_values_ref;
    }

    ## Build arguments to submit to function
    my @args;
    while ( my ( $key_index, $key ) = each @keys ) {

        push @args, $keys[$key_index], $values[$key_index];
    }
    return @args;
}

sub _test_base_command {

##_test_base_command

##Function : Test the function base command. Executable, ".jar" etc.
##Returns  : ""
##Arguments: $base_command, $expected_base_command
##         : $base_command          => First word in command line usually name of executable
##         : $expected_base_command => Expected base command

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $base_command;
    my $expected_base_command;

    my $tmpl = {
        base_command => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$base_command
        },
        expected_base_command => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$expected_base_command
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw(Could not parse arguments!);

    if ( $base_command ne $expected_base_command ) {

        ok(
            $base_command eq $expected_base_command,
            'Argument: ' . $expected_base_command
        );
        exit 1;
    }
    return;
}

sub _test_write_to_file {

## _test_write_to_file

##Function : Test of writing to file using supplied FILEHANDLE
##Returns  : ""
##Arguments: $module_function_cref, $args_ref, $base_command
##         : $module_function_cref => Module method to test
##         : $args_ref             => Arguments to function call
##         : $base_command         => First word in command line usually name of executable

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $module_function_cref;
    my $args_ref;
    my $base_command;

    my $tmpl = {
        module_function_cref =>
          { required => 1, defined => 1, store => \$module_function_cref },
        args_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$args_ref
        },
        base_command => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$base_command
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw(Could not parse arguments!);

    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Add new FILEHANDLE to args
    push @{$args_ref}, 'FILEHANDLE', $FILEHANDLE;

    # For storing info to write
    my $file_content;

    ## Store file content in memory by using referenced variable
    open $FILEHANDLE, '>', \$file_content
      or croak 'Cannot write to ' . $file_content . ': ' . $OS_ERROR;

    $module_function_cref->( { @{$args_ref} } );

    close $FILEHANDLE;

    ## Perform test
    ok( $file_content =~ /^$base_command/, 'Write commands to file' );

    return;
}

1;
