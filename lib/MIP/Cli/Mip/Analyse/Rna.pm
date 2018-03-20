package MIP::Cli::Mip::Analyse::Rna;

use Carp;
use File::Spec::Functions qw{ catfile };
use FindBin qw{ $Bin };
use open qw( :encoding(UTF-8) :std );
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use MooseX::App::Command;
use MooseX::Types::Moose qw{ Str Int HashRef Num Bool ArrayRef };
use Moose::Util::TypeConstraints;

## MIPs lib
use MIP::Main::Analyse qw{ mip_analyse };

our $VERSION = 0.01;

extends(qw{ MIP::Cli::Mip::Analyse });

command_short_description(q{Rna analysis});

command_long_description(q{Rna analysis on wts sequence data});

command_usage(q{mip <analyse> <rna> <family_id> --config <config_file> });

## Define, check and get Cli supplied parameters
_build_usage();

sub run {
    my ($arg_href) = @_;

    ## Remove Moose::App extra variable
    delete $arg_href->{extra_argv};

    ## Input from Cli
    my %active_parameter = %{$arg_href};

    use MIP::File::Format::Parameter qw{ parse_definition_file  };
    use MIP::File::Format::Yaml qw{ order_parameter_names };
    use MIP::Get::Analysis qw{ print_program };

    ## Mip analyse rna parameters
    ## The order of files in @definition_files should follow commands inheritance
    my @definition_files = (
        catfile( $Bin, qw{ definitions mip_parameters.yaml } ),
        catfile( $Bin, qw{ definitions rna_parameters.yaml } ),
    );

    ## Non mandatory parameter definition keys to check
    my $non_mandatory_parameter_keys_path =
      catfile( $Bin, qw{ definitions non_mandatory_parameter_keys.yaml } );

    ## Mandatory parameter definition keys to check
    my $mandatory_parameter_keys_path =
      catfile( $Bin, qw{ definitions mandatory_parameter_keys.yaml } );

    ### %parameter holds all defined parameters for MIP
    ### mip analyse rna
    my %parameter;
    foreach my $definition_file (@definition_files) {

        %parameter = (
            %parameter,
            parse_definition_file(
                {
                    define_parameters_path => $definition_file,
                    non_mandatory_parameter_keys_path =>
                      $non_mandatory_parameter_keys_path,
                    mandatory_parameter_keys_path =>
                      $mandatory_parameter_keys_path,
                }
            )
        );
    }

    ## Print programs and exit
    if ( $active_parameter{print_programs} ) {

        print_program(
            {
                define_parameters_files_ref => \@definition_files,
                parameter_href              => \%parameter,
                print_program_mode => $active_parameter{print_program_mode},
            }
        );
        exit;
    }

    ### To add/write parameters in the correct order
    ## Adds the order of first level keys from yaml file to array
    my @order_parameters;
    foreach my $define_parameters_file (@definition_files) {

        push @order_parameters,
          order_parameter_names(
            {
                file_path => $define_parameters_file,
            }
          );
    }

    ## File info hash
    my %file_info = (

        # Human genome meta files
        human_genome_reference_file_endings => [qw{ .dict .fai }],
    );

    mip_analyse(
        {
            active_parameter_href => \%active_parameter,
            file_info_href        => \%file_info,
            parameter_href        => \%parameter,
            order_parameters_ref  => \@order_parameters,
        }
    );

    return;
}

sub _build_usage {

## Function : Get and/or set input parameters
## Returns  :
## Arguments:
    option(
        q{psalmon_quant} => (
            cmd_aliases   => [qw{ psqt }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{quantify transcripts using salmon},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );
    option(
        q{lib} => (
            cmd_aliases   => [qw{ psqt_bob }],
            cmd_tags      => [q{Default: ISF}],
            documentation => q{Library orientation and strandedness},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{pstar_aln} => (
            cmd_aliases   => [qw{ pstn }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Align reads using Star aln},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{align_intron_max} => (
            cmd_aliases   => [qw{ stn_aim }],
            cmd_tags      => [q{Default: 100,000}],
            documentation => q{Maximum intron size},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{align_mates_gap_max} => (
            cmd_aliases   => [qw{ stn_amg }],
            cmd_tags      => [q{Default: 100,000}],
            documentation => q{Maximum gap between two mates},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{align_sjdb_overhang_min} => (
            cmd_aliases => [qw{ stn_asom }],
            cmd_tags    => [q{Default: 10}],
            documentation =>
              q{Minimum overhang (i.e. block size) for spliced alignments},
            is  => q{rw},
            isa => Int,
        )
    );

    option(
        q{chim_junction_overhang_min} => (
            cmd_aliases   => [qw{ stn_cjom }],
            cmd_tags      => [q{Default: 12}],
            documentation => q{Minimum overhang for a chimeric junction},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{chim_segment_min} => (
            cmd_aliases   => [qw{ stn_csm }],
            cmd_tags      => [q{Default: 12}],
            documentation => q{Minimum length of chimaeric segment},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{two_pass_mode} => (
            cmd_aliases   => [qw{ stn_tpm }],
            cmd_tags      => [q{Default: Basic}],
            documentation => q{Two pass mode setting},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{pstar_fusion} => (
            cmd_aliases   => [qw{ pstf }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Detect fusion transcripts with star fusion},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    return;
}

1;
