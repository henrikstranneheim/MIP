package MIP::Program::Salmon;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $SPACE };
use MIP::Environment::Executable qw{ get_executable_base_command };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ salmon_index salmon_quant };
}

Readonly my $BASE_COMMAND => q{salmon};

sub salmon_index {

## Function  : Perl wrapper for Salmon index, version 1.4.0.
## Returns   : @commands
## Arguments : $decoy_path             => Decoy sequence ids
##           : $fasta_path             => Input reference fasta path, note salmon does not use the genome reference fasta, it uses a fasta file of transcripts
##           : $filehandle             => Filehandle to write to
##           : $gencode                => Transcripts are in gencode format
##           : $outfile_path           => Outfile path
##           : $stderrfile_path        => Stderrfile path
##           : $stderrfile_path_append => Append stderr info to file path
##           : $stdoutfile_path        => Stdoutfile path
##           : $temp_directory         => Temporary directory
##           : $threads                => Threads used for indexing

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $decoy_path;
    my $fasta_path;
    my $filehandle;
    my $gencode;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $temp_directory;
    my $threads;

    my $tmpl = {
        decoy_path => {
            store       => \$decoy_path,
            strict_type => 1,
        },
        fasta_path => {
            defined     => 1,
            required    => 1,
            store       => \$fasta_path,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        gencode => {
            allow       => [ undef, 0, 1 ],
            store       => \$gencode,
            strict_type => 1,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        stderrfile_path => {
            store       => \$stderrfile_path,
            strict_type => 1,
        },
        stderrfile_path_append => {
            store       => \$stderrfile_path_append,
            strict_type => 1,
        },
        stdoutfile_path => {
            store       => \$stdoutfile_path,
            strict_type => 1,
        },
        temp_directory => {
            store       => \$temp_directory,
            strict_type => 1,
        },
        threads => {
            allow       => [ undef, qr/\A \d+ \z/xms ],
            store       => \$threads,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands =
      ( get_executable_base_command( { base_command => $BASE_COMMAND, } ), qw{ index } );

    push @commands, q{--transcripts} . $SPACE . $fasta_path;

    push @commands, q{--index} . $SPACE . $outfile_path;

    if ($decoy_path) {

        push @commands, q{--decoy} . $SPACE . $decoy_path;
    }

    if ($gencode) {

        push @commands, q{--gencode};
    }

    if ($temp_directory) {

        push @commands, q{--tmpdir} . $SPACE . $temp_directory;
    }

    if ($threads) {
        push @commands, q{--threads} . $SPACE . $threads,;
    }

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
            stdoutfile_path        => $stdoutfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            filehandle   => $filehandle,
            separator    => $SPACE,

        }
    );
    return @commands;
}

sub salmon_quant {

## Function  : Perl wrapper for Salmon quant, version 0.9.1.
## Returns   : @commands
## Arguments : $filehandle             => Filehandle to write to
##           : $gc_bias                => Correct for GC-bias
##           : $index_path             => Path to the index folder
##           : $libi_type              => Library visit the salmon website for more  info
##           : $outdir_path            => Path of the output directory
##           : $read_1_fastq_paths_ref => Read 1 Fastq paths
##           : $read_2_fastq_paths_ref => Read 2 Fastq paths
##           : $read_files_command     => command applied to the input FASTQ files
##           : $stderrfile_path        => Stderrfile path
##           : $stderrfile_path_append => Append stderr info to file path
##           : $stdoutfile_path        => Stdoutfile path
##           : $validate_mappings      => Validate mappings

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $gc_bias;
    my $index_path;
    my $outdir_path;
    my $read_1_fastq_paths_ref;
    my $read_2_fastq_paths_ref;
    my $read_files_command;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $lib_type;
    my $validate_mappings;

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
        },
        gc_bias => {
            store       => \$gc_bias,
            strict_type => 1,
        },
        index_path => {
            defined     => 1,
            required    => 1,
            store       => \$index_path,
            strict_type => 1,
        },
        lib_type => {
            allow       => [qw{ A ISF ISR MSF MSR OSR OSF }],
            default     => q{A},
            store       => \$lib_type,
            strict_type => 1,
        },
        outdir_path => {
            defined     => 1,
            required    => 1,
            store       => \$outdir_path,
            strict_type => 1,
        },
        read_1_fastq_paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$read_1_fastq_paths_ref,
            strict_type => 1,
        },
        read_2_fastq_paths_ref => {
            default     => [],
            defined     => 1,
            store       => \$read_2_fastq_paths_ref,
            strict_type => 1,
        },
        read_files_command => {
            required    => 1,
            defined     => 1,
            store       => \$read_files_command,
            strict_type => 1,
        },
        stderrfile_path => {
            store       => \$stderrfile_path,
            strict_type => 1,
        },
        stderrfile_path_append => {
            store       => \$stderrfile_path_append,
            strict_type => 1,
        },
        stdoutfile_path => {
            store       => \$stdoutfile_path,
            strict_type => 1,
        },
        validate_mappings => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$validate_mappings,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands =
      ( get_executable_base_command( { base_command => $BASE_COMMAND, } ), qw{ quant } );

    if ($gc_bias) {
        push @commands, q{--gcBias};
    }

    push @commands, q{--index} . $SPACE . $index_path;

# Library type, defines if the library is stranded or not, and the orientation of the reads, according to the documentation http://salmon.readthedocs.io/en/latest/library_type.html
    push @commands, q{--libType} . $SPACE . $lib_type;

    push @commands, q{--output} . $SPACE . $outdir_path;

    if ($validate_mappings) {
        push @commands, q{--validateMappings};
    }

# The input Fastq files, either single reads or paired. Salmon uses a bash command to stream the reads. Here, the default is <( pigz -dc file.fastq.gz )
    push @commands,
        q{-1}
      . $SPACE . q{<(}
      . $read_files_command
      . $SPACE
      . join( $SPACE, @{$read_1_fastq_paths_ref} )
      . $SPACE . q{)};

    if ( @{$read_2_fastq_paths_ref} ) {
        push @commands,
            q{-2}
          . $SPACE . q{<(}
          . $read_files_command
          . $SPACE
          . join( $SPACE, @{$read_2_fastq_paths_ref} )
          . $SPACE . q{)};
    }

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
            stdoutfile_path        => $stdoutfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            filehandle   => $filehandle,
            separator    => $SPACE,

        }
    );
    return @commands;
}

1;
