package MIP::Contigs;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename };
use File::Spec::Functions qw{ catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use List::Util qw{ any };

## MIPs lib/
use MIP::Constants qw{ $AMPERSAND $LOG_NAME $NEWLINE %PRIMARY_CONTIG $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      check_select_file_contigs
      delete_contig_elements
      delete_non_wes_contig
      delete_male_contig
      generate_contig_interval_file
      get_contig_set
      set_contigs
      sort_contigs_to_contig_set
      update_contigs_for_run
    };

}

sub check_select_file_contigs {

## Function : Check that select file contigs is a subset of primary contigs
## Returns  :
## Arguments: $contigs_ref             => Primary contigs of the human genome reference
##          : $select_file_contigs_ref => Select file contigs

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $contigs_ref;
    my $select_file_contigs_ref;

    my $tmpl = {
        contigs_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$contigs_ref,
            strict_type => 1,
        },
        select_file_contigs_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$select_file_contigs_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use Array::Utils qw{ array_minus };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Check that select file contigs are a subset of primary contigs
    my @unique_select_contigs =
      array_minus( @{$select_file_contigs_ref}, @{$contigs_ref} );

    if (@unique_select_contigs) {

        $log->fatal( q{Option 'vcfparser_select_file' contig(s): } . join $SPACE,
            @unique_select_contigs );
        $log->fatal( q{Is not a subset of the human genome reference contigs: } . join $SPACE,
            @{$contigs_ref} );
        exit 1;
    }
    return 1;
}

sub delete_contig_elements {

## Function : Return new cleansed contig array without remove contigs elements while leaving original contigs_ref
##            untouched
## Returns  : @cleansed_contigs
## Arguments: contigs_ref         => Contigs array to remove element(s) from {REF}
##          : $remove_contigs_ref => Remove these contigs {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $contigs_ref;
    my $remove_contigs_ref;

    my $tmpl = {
        contigs_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$contigs_ref,
            strict_type => 1,
        },
        remove_contigs_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$remove_contigs_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my %contig_to_remove = map { $_ => 1 } @{$remove_contigs_ref};

    ## Special case for Mitochondria
    if ( exists $contig_to_remove{MT} ) {
        $contig_to_remove{M} = 1;
    }

    ### Make sure that contig is removed independent of genome source i.e prefix or not
    ## If contigs has prefix
    if ( defined $contigs_ref->[0]
        && $contigs_ref->[0] =~ / \A chr /xsm )
    {

        ## And remove contigs has not
        if ( defined $remove_contigs_ref->[0]
            && $remove_contigs_ref->[0] !~ / \A chr /xsm )
        {

            ## Add prefix to keys to match contigs prefix
            %contig_to_remove = map { q{chr} . $_ => 1 } keys %contig_to_remove;
        }
    }
    my @cleansed_contigs = grep { not $contig_to_remove{$_} } @{$contigs_ref};
    return @cleansed_contigs;
}

sub delete_non_wes_contig {

## Function : Delete contig chrM | MT from contigs array if consensus analysis type is wes
## Returns  : @contigs
## Arguments: $consensus_analysis_type => Consensus analysis_type
##          : $contigs_ref             => Contigs array to update {REF}
##          : $contig_names_ref        => Contig names to remove {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $consensus_analysis_type;
    my $contigs_ref;
    my $contig_names_ref;

    my $tmpl = {
        consensus_analysis_type => {
            defined     => 1,
            required    => 1,
            store       => \$consensus_analysis_type,
            strict_type => 1,
        },
        contigs_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$contigs_ref,
            strict_type => 1,
        },
        contig_names_ref => {
            allow => [
                sub {
                    check_allowed_array_values(
                        {
                            allowed_values_ref => [qw{ chrM M MT }],
                            values_ref         => $arg_href->{contig_names_ref},
                        }
                    );
                }
            ],
            default     => [qw{ chrM M MT }],
            store       => \$contig_names_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::List qw{ check_allowed_array_values };
    use MIP::Contigs qw{ delete_contig_elements };

    return @{$contigs_ref}
      if ( $consensus_analysis_type eq q{wgs} || $consensus_analysis_type eq q{wts} );

    ## Removes contig M | chrMT from contigs
    my @contigs = delete_contig_elements(
        {
            contigs_ref        => \@{$contigs_ref},
            remove_contigs_ref => $contig_names_ref,
        }
    );
    return @contigs;
}

sub delete_male_contig {

## Function : Delete contig chrY | Y from contigs array if no male or other found
## Returns  : @contigs
## Arguments: $contigs_ref => Contigs array to update {REF}
##          : $include_y   => Male(s) was included in the analysis

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $contigs_ref;
    my $include_y;

    my $tmpl = {
        contigs_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$contigs_ref,
            strict_type => 1,
        },
        include_y => {
            allow       => [ 0, 1 ],
            defined     => 1,
            required    => 1,
            store       => \$include_y,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return @{$contigs_ref} if ($include_y);

    ## Removes contig Y | chrY from contigs if no males or 'other' found in analysis
    my @contigs = delete_contig_elements(
        {
            contigs_ref        => $contigs_ref,
            remove_contigs_ref => [qw{ Y chrY }],
        }
    );
    return @contigs;
}

sub generate_contig_interval_file {

## Function : Generate contig specific interval list files
## Returns  : %bed_file_path
## Arguments: $contigs_ref           => Contigs to split in files
##          : $exome_target_bed_file => Interval file to split
##          : $file_ending           => File ending to add {Optional}
##          : $filehandle            => Filehandle to write to
##          : $max_process_number    => Maximum processes to use
##          : $outdirectory          => Outdirectory
##          : $reference_dir         => MIP reference directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $contigs_ref;
    my $exome_target_bed_file;
    my $file_ending;
    my $filehandle;
    my $outdirectory;
    my $reference_dir;

    ## Default(s)
    my $max_process_number;

    my $tmpl = {
        contigs_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$contigs_ref,
            strict_type => 1,
        },
        exome_target_bed_file => {
            defined     => 1,
            required    => 1,
            store       => \$exome_target_bed_file,
            strict_type => 1,
        },
        file_ending        => { strict_type => 1, store => \$file_ending },
        filehandle         => { store       => \$filehandle, },
        max_process_number => {
            allow       => qr/ \A \d+ \z /sxm,
            default     => 1,
            store       => \$max_process_number,
            strict_type => 1,
        },
        outdirectory => {
            defined     => 1,
            required    => 1,
            store       => \$outdirectory,
            strict_type => 1,
        },
        reference_dir => {
            defined     => 1,
            required    => 1,
            store       => \$reference_dir,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Processmanagement::Processes qw{ print_wait };

    my %bed_file_path;
    my $process_batches_count = 1;

    say {$filehandle} q{## Generate contig specific interval_list}, $NEWLINE;

  CONTIG:
    while ( my ( $contig_index, $contig ) = each @{$contigs_ref} ) {

        $process_batches_count = print_wait(
            {
                filehandle            => $filehandle,
                max_process_number    => $max_process_number,
                process_batches_count => $process_batches_count,
                process_counter       => $contig_index,
            }
        );

        ## Splits a target file into new contig specific target file
        my $contig_bed_file_path = _split_interval_file_contigs(
            {
                contig       => $contig,
                file_ending  => $file_ending,
                filehandle   => $filehandle,
                indirectory  => $reference_dir,
                infile       => basename($exome_target_bed_file),
                outdirectory => $outdirectory,
            }
        );
        $bed_file_path{$contig} = [$contig_bed_file_path];
    }
    say {$filehandle} q{wait}, $NEWLINE;

    return %bed_file_path;
}

sub get_contig_set {

## Function : Get contig set per genome build version
## Returns  : @{$PRIMARY_CONTIG{$version}{$contig_set}} | %{$PRIMARY_CONTIG{$version}{$contig_set}} | %{$PRIMARY_CONTIG{$version}}
## Arguments: $contig_set => Name of contig set to get
##          : $version    => Version of the human genome reference

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $contig_set;
    my $version;

    my $tmpl = {
        contig_set => {
            allow       => [qw{ contigs contigs_size_ordered synonyms_map }],
            defined     => 1,
            store       => \$contig_set,
            strict_type => 1,
        },
        version => {
            defined     => 1,
            required    => 1,
            store       => \$version,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return if ( not exists $PRIMARY_CONTIG{$version} );

    if ($contig_set) {

        return @{ $PRIMARY_CONTIG{$version}{$contig_set} }
          if ( ref $PRIMARY_CONTIG{$version}{$contig_set} eq q{ARRAY} );

        return %{ $PRIMARY_CONTIG{$version}{$contig_set} }
          if ( ref $PRIMARY_CONTIG{$version}{$contig_set} eq q{HASH} );

    }
    return %{ $PRIMARY_CONTIG{$version} };
}

sub set_contigs {

## Function : Set contig prefix and contig names depending on reference used.
## Returns  :
## Arguments: $file_info_href => File info hash {REF}
##          : $version        => Version of the human genome reference

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_info_href;
    my $version;

    my $tmpl = {
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        version => {
            defined     => 1,
            required    => 1,
            store       => \$version,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File_info qw{ set_alt_loci_contigs set_primary_contigs };

    ## Make a modifiable copy for downstream use of global constant
    my %primary_contig_clone = Readonly::Clone %PRIMARY_CONTIG;

    ## Get hash of genome build version primary assembly contigs
    my %primary_contig;
    @primary_contig{ @{ $primary_contig_clone{$version}{contigs} } } = ();

    ## Set alternative loci contig set
    set_alt_loci_contigs(
        {
            alt_contig_set_name => q{alt_loci},
            file_info_href      => $file_info_href,
            primary_contig_href => \%primary_contig,
        }
    );

    ## Set contigs sets for bam level processing downstream
    my %contig_set = (
        bam_contigs              => q{contigs},
        bam_contigs_size_ordered => q{contigs_size_ordered},
        contigs                  => q{contigs},
        contigs_size_ordered     => q{contigs_size_ordered},
        primary_contigs          => q{contigs},
    );

  CONTIG_SET:
    while ( my ( $contig_set_name, $contig_set ) = each %contig_set ) {

        set_primary_contigs(
            {
                file_info_href      => $file_info_href,
                primary_contigs_ref => \@{ $primary_contig_clone{$version}{$contig_set} },
                contig_set_name     => $contig_set_name,
            }
        );
    }
    return;
}

sub sort_contigs_to_contig_set {

## Function : Sorts array depending on reference array. NOTE: Only entries present in reference array will survive in sorted array.
## Returns  : @sorted_contigs
## Arguments: $consensus_analysis_type    => Consensus analysis_type {REF}
##          : $sort_contigs_ref           => Contigs to sort according to reference contig set
##          : $sort_reference_contigs_ref => Contigs to use as reference map when sorting

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $consensus_analysis_type;
    my $sort_contigs_ref;
    my $sort_reference_contigs_ref;

    my $tmpl = {
        consensus_analysis_type => {
            defined     => 1,
            required    => 1,
            store       => \$consensus_analysis_type,
            strict_type => 1,
        },
        sort_contigs_ref => {
            default  => [],
            required => 1,
            store    => \$sort_contigs_ref,
        },
        sort_reference_contigs_ref => {
            default  => [],
            defined  => 1,
            required => 1,
            store    => \$sort_reference_contigs_ref,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my @sorted_contigs;

    ## Sanity check
    if ( not @{$sort_contigs_ref} ) {

        $log->fatal(q{Nothing to sort in supplied sort contigs set });
        exit 1;
    }
    if ( not @{$sort_reference_contigs_ref} ) {

        $log->fatal(q{No contigs to use in supplied sort reference contig set });
        exit 1;
    }

    ## Sort the contigs depending on reference array
  REF_ELEMENT:
    foreach my $element ( @{$sort_reference_contigs_ref} ) {

        ## If present in hash of array to sort push to @sorted_contigs
        if ( any { $_ eq $element } @{$sort_contigs_ref} ) {

            push @sorted_contigs, $element;
        }
    }

    ## Test if all contigs collected from select file was sorted by reference contig array
    if ( @sorted_contigs
        and scalar @{$sort_contigs_ref} != scalar @sorted_contigs )
    {

      SORT_ELEMENT:
        foreach my $element ( @{$sort_contigs_ref} ) {

            ## If element is not part of array
            if ( not any { $_ eq $element } @sorted_contigs ) {

                ## Special case when analysing wes since Mitochondrial contigs have no baits in exome capture kits
                next SORT_ELEMENT
                  if (  $consensus_analysis_type eq q{wes}
                    and $element =~ / MT$ | M$ /sxm );

                $log->fatal( q{Could not detect 'contig'= }
                      . $element
                      . q{ from column 1 in '-vcfparser_select_file' in reference contigs collected from '-human_genome_reference'}
                );
                exit 1;
            }
        }
    }
    return @sorted_contigs;
}

sub update_contigs_for_run {

## Function : Update contigs depending on settings in run
## Returns  :
## Arguments: $consensus_analysis_type => Consensus analysis_type
##          : $exclude_contigs_ref     => Exclude contigs from analysis {REF}
##          : $file_info_href          => File info hash {REF}
##          : $include_y               => Male was included in the analysis

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $consensus_analysis_type;
    my $exclude_contigs_ref;
    my $file_info_href;
    my $include_y;

    my $tmpl = {
        consensus_analysis_type => {
            defined     => 1,
            required    => 1,
            store       => \$consensus_analysis_type,
            strict_type => 1,
        },
        exclude_contigs_ref => {
            default     => [],
            required    => 1,
            store       => \$exclude_contigs_ref,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        include_y => {
            allow       => [ 0, 1 ],
            defined     => 1,
            required    => 1,
            store       => \$include_y,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @contig_sets = (
        \@{ $file_info_href->{bam_contigs} },
        \@{ $file_info_href->{bam_contigs_size_ordered} },
        \@{ $file_info_href->{contigs} },
        \@{ $file_info_href->{contigs_size_ordered} },
        \@{ $file_info_href->{select_file_contigs} },
    );

  CONTIG_REF:
    foreach my $contigs_ref (@contig_sets) {

        ## Delete user specified contigs from contigs array
        @{$contigs_ref} = delete_contig_elements(
            {
                contigs_ref        => $contigs_ref,
                remove_contigs_ref => $exclude_contigs_ref,
            }
        );

        ## Delete contig chrM|MT from contigs array if consensus analysis type is wes
        @{$contigs_ref} = delete_non_wes_contig(
            {
                consensus_analysis_type => $consensus_analysis_type,
                contigs_ref             => $contigs_ref,
            }
        );
    }

    my @male_contig_arrays = (
        \@{ $file_info_href->{contigs} },
        \@{ $file_info_href->{contigs_size_ordered} },
        \@{ $file_info_href->{select_file_contigs} },
    );

  ARRAY_REF:
    foreach my $male_contigs_ref (@male_contig_arrays) {

        ## Removes contig_names from contigs array if no male or 'other' found
        @{$male_contigs_ref} = delete_male_contig(
            {
                contigs_ref => $male_contigs_ref,
                include_y   => $include_y,
            }
        );
    }
    return;
}

sub _split_interval_file_contigs {

## Function : Splits a target file into new contig specific target file
## Returns  : $outfile_path
## Arguments: $file_ending  => File ending to add
##          : $contig       => Contig to extract
##          : $filehandle   => filehandle to write to
##          : $indirectory  => Indirectory
##          : $infile       => Target file
##          : $outdirectory => Outdirectory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $contig;
    my $file_ending;
    my $filehandle;
    my $indirectory;
    my $infile;
    my $outdirectory;

    my $tmpl = {
        contig      => { store   => \$contig, },
        file_ending => { store   => \$file_ending, strict_type => 1, },
        filehandle  => { defined => 1, required => 1, store => \$filehandle, },
        indirectory => {
            defined     => 1,
            required    => 1,
            store       => \$indirectory,
            strict_type => 1,
        },
        infile => {
            defined     => 1,
            required    => 1,
            store       => \$infile,
            strict_type => 1,
        },
        outdirectory => {
            defined     => 1,
            required    => 1,
            store       => \$outdirectory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Language::Perl qw{ perl_nae_oneliners };

    my $outfile_path =
      defined $contig ? catfile( $outdirectory, $contig . $UNDERSCORE . $infile ) : undef;

    if ( defined $contig ) {

        $outfile_path = defined $file_ending ? $outfile_path . $file_ending : $outfile_path;

        perl_nae_oneliners(
            {
                filehandle         => $filehandle,
                oneliner_name      => q{write_header_for_contig},
                oneliner_parameter => $contig,
                stdinfile_path     => catfile( $indirectory, $infile ),
                stdoutfile_path    => $outfile_path . $SPACE . $AMPERSAND,
                use_container      => 1,
            }
        );
        say {$filehandle} $NEWLINE;
    }
    return $outfile_path;
}

1;
