package MIP::Language::Perl;

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
use MIP::Constants qw{ $BACKWARD_SLASH $DASH $DOT $COLON $NEWLINE $SPACE $SINGLE_QUOTE };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    our @EXPORT_OK =
      qw{ check_modules_existance get_cpan_file_modules perl_base perl_nae_oneliners };
}

Readonly my $MINUS_ONE => -1;

sub check_modules_existance {

## Function : Evaluate that all perl modules required by MIP are installed
## Returns  :
## Arguments: $modules_ref  => Array of module names
##          : $program_name => Program name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $modules_ref;
    my $program_name;

    my $tmpl = {
        modules_ref => {
            default     => [],
            required    => 1,
            store       => \$modules_ref,
            strict_type => 1,
        },
        program_name => {
            defined     => 1,
            required    => 1,
            store       => \$program_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    require Try::Tiny;
    use Try::Tiny;

  MODULE:
    foreach my $module ( @{$modules_ref} ) {

        ## Special case for Readonly::XS since it is not a standalone module
        $module =~ s{Readonly::XS}{Readonly}sxmg;

        ## Replace "::" with "/" since the automatic replacement magic only occurs for bare words.
        $module =~ s{::}{/}sxmg;

        ## Add perl module ending for the same reason
        $module .= $DOT . q{pm};

        try {
            require $module;
        }
        catch {
            say {*STDERR} qq{FATAL: $module not installed - Please install to run $program_name};
            croak(q{FATAL: Aborting!});
        };
    }
    return 1;
}

sub get_cpan_file_modules {

## Function : Get perl modules from cpan file
## Returns  : @cpanm_modules
## Arguments: $cpanfile_path => Path to cpanfile

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $cpanfile_path;

    my $tmpl = {
        cpanfile_path => {
            default     => 1,
            required    => 1,
            store       => \$cpanfile_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use Module::CPANfile;
    use CPAN::Meta::Prereqs;

    ## Load cpanfile
    my $file = Module::CPANfile->load($cpanfile_path);
    ## Get hash_ref without objects
    my $file_href = $file->prereqs->as_string_hash;

    ## Get cpanm modules
    my @cpanm_modules = sort keys %{ $file_href->{runtime}{requires} };

    return @cpanm_modules;
}

sub perl_base {

## Function : Perl base and switches
## Returns  :
## Arguments: $autosplit     => Turns on autosplit mode when used with a -n or -p
##          : $command_line  => Enter one line of program
##          : $inplace       => In place edit
##          : $n             => Iterate over filename arguments
##          : $print         => Print line
##          : $print_newline => Print newline at the end of line
##          : $use_container => Use container perl instead of main (this) process perl

    my ($arg_href) = @_;

    ## Flatten argument(s)

    ## Default(s)
    my $autosplit;
    my $command_line;
    my $inplace;
    my $n;
    my $print;
    my $print_newline;
    my $use_container;

    my $tmpl = {
        autosplit => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$autosplit,
            strict_type => 1,
        },
        command_line => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$command_line,
            strict_type => 1,
        },
        inplace => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$inplace,
            strict_type => 1,
        },
        n => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$n,
            strict_type => 1,
        },
        print => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$print,
            strict_type => 1,
        },
        print_newline => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$print_newline,
            strict_type => 1,
        },
        use_container => => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$use_container,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Environment::Executable qw{ get_executable_base_command };

    my @commands =
      $use_container
      ? ( get_executable_base_command( { base_command => q{perl}, } ), )
      : qw{ perl };

    if ($n) {

        push @commands, q{-n};
    }
    if ($autosplit) {

        push @commands, q{-a};
    }
    if ($inplace) {

        push @commands, q{-i};
    }
    if ($print) {

        push @commands, q{-p};
    }
    if ($print_newline) {

        push @commands, q{-l};
    }
    if ($command_line) {

        push @commands, q{-e};
    }

    return @commands;
}

sub perl_nae_oneliners {

## Function : Return predifined one liners
## Returns  : @commands
## Arguments: $autosplit              => Turns on autosplit mode when used with a -n or -p
##          : $command_line           => Enter one line of program
##          : $escape_oneliner        => Escape perl oneliner program
##          : $filehandle             => Filehandle to write to
##          : $n                      => Iterate over filename arguments
##          : $oneliner_cmd           => Command to execute
##          : $oneliner_name          => Perl oneliner name
##          : $oneliner_parameter     => Feed a parameter to the oneliner program
##          : $print_newline          => Print newline at end of line
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdinfile_path         => Stdinfile path
##          : $stdoutfile_path        => Stdoutfile path
##          : $stdoutfile_path_append => Append stdout info to file path
##          : $use_container          => Use container perl instead of main (this) process perl

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $escape_oneliner;
    my $filehandle;
    my $oneliner_cmd;
    my $oneliner_name;
    my $oneliner_parameter;
    my $print_newline;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdinfile_path;
    my $stdoutfile_path;
    my $stdoutfile_path_append;

    ## Default(s)
    my $autosplit;
    my $command_line;
    my $n;
    my $use_container;

    my $tmpl = {
        autosplit => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$autosplit,
            strict_type => 1,
        },
        command_line => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$command_line,
            strict_type => 1,
        },
        escape_oneliner => {
            allow       => [ undef, 0, 1 ],
            store       => \$escape_oneliner,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        n => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$n,
            strict_type => 1,
        },
        oneliner_cmd => {
            store       => \$oneliner_cmd,
            strict_type => 1,
        },
        oneliner_name => {
            store       => \$oneliner_name,
            strict_type => 1,
        },
        oneliner_parameter => {
            store       => \$oneliner_parameter,
            strict_type => 1,
        },
        print_newline => {
            allow       => [ undef, 0, 1 ],
            store       => \$print_newline,
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
        stdinfile_path  => { store => \$stdinfile_path, strict_type => 1, },
        stdoutfile_path => {
            store       => \$stdoutfile_path,
            strict_type => 1,
        },
        stdoutfile_path_append => {
            store       => \$stdoutfile_path_append,
            strict_type => 1,
        },
        use_container => => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$use_container,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Oneliner dispatch table
    my %oneliner = (
        add_binary                           => \&_add_binary,
        add_pct_mapped_reads_from_samtools   => \&_add_pct_mapped_reads_from_samtools,
        bcftools_norm_check                  => \&_bcftools_norm_check,
        build_md5sum_check                   => \&_build_md5sum_check,
        genepred_to_refflat                  => \&_genepred_to_refflat,
        get_dict_contigs                     => \&_get_dict_contigs,
        q{get_fastq_header_v1.4}             => \&_get_fastq_header_v1_4,
        q{get_fastq_header_v1.4_interleaved} => \&_get_fastq_header_v1_4_interleaved,
        q{get_fastq_header_v1.8}             => \&_get_fastq_header_v1_8,
        q{get_fastq_header_v1.8_interleaved} => \&_get_fastq_header_v1_8_interleaved,
        get_fastq_read_length                => \&_get_fastq_read_length,
        get_gene_panel_header                => \&_get_gene_panel_info,
        get_gene_panel_hgnc_symbols          => \&_get_gene_panel_hgnc_symbols,
        get_rrna_transcripts                 => \&_get_rrna_transcripts,
        get_select_contigs_by_col            => \&_get_select_contigs_by_col,
        get_vcf_header_id_line               => \&_get_vcf_header_id_line,
        get_vcf_loqusdb_headers              => \&_get_vcf_loqusdb_headers,
        get_vcf_sample_ids                   => \&_get_vcf_sample_ids,
        reformat_sacct_headers               => \&_reformat_sacct_headers,
        remove_decomposed_asterisk_records   => \&_remove_decomposed_asterisk_records,
        star_sj_tab_to_bed                   => \&_star_sj_tab_to_bed,
        synonyms_grch37_to_grch38            => \&_synonyms_grch37_to_grch38,
        synonyms_grch38_to_grch37            => \&_synonyms_grch38_to_grch37,
        write_header_for_contig              => \&_write_header_for_contig,
        write_contigs_size_file              => \&_write_contigs_size_file,
    );

    my %oneliner_option = (
        add_binary => {
            binary => $oneliner_parameter,
        },
        build_md5sum_check => {
            file_path => $oneliner_parameter,
        },
        get_vcf_header_id_line => {
            id => $oneliner_parameter,
        },
        reformat_sacct_headers  => { sacct_header => $oneliner_parameter, },
        write_header_for_contig => => {
            contig => $oneliner_parameter,
        },
    );

    ## Stores commands depending on input parameters
    my @commands = perl_base(
        {
            autosplit     => $autosplit,
            command_line  => $command_line,
            n             => $n,
            print_newline => $print_newline,
            use_container => $use_container,
        }
    );

    ## Fetch oneliner from dispatch table
    if (    defined $oneliner_name
        and exists $oneliner{$oneliner_name}
        and not $oneliner_cmd )
    {

        $oneliner_cmd = $oneliner{$oneliner_name}->(
            defined $oneliner_option{$oneliner_name}
            ? $oneliner_option{$oneliner_name}
            : undef
        );
    }

    ## Quote oneliner for use with xargs
    if ( $oneliner_cmd and $escape_oneliner ) {

        $oneliner_cmd = $BACKWARD_SLASH . $oneliner_cmd;
        substr $oneliner_cmd, $MINUS_ONE, 0, $BACKWARD_SLASH;
    }

    if ($oneliner_cmd) {

        push @commands, $oneliner_cmd;
    }

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
            stdinfile_path         => $stdinfile_path,
            stdoutfile_path        => $stdoutfile_path,
            stdoutfile_path_append => $stdoutfile_path_append,
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

sub _add_binary {

## Function : Add "binary: " to stdin
## Returns  : $add_binary
## Arguments: $binary => Binary to add

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $binary;

    my $tmpl = {
        binary => {
            defined     => 1,
            required    => 1,
            store       => \$binary,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $add_binary = q?'say STDOUT qq{? . $binary . $COLON . $SPACE . q?$_}'?;

    return $add_binary;
}

sub _add_pct_mapped_reads_from_samtools {

## Function : Calculate the percentage of raw total sequences and reads mapped from samtools stats
## Returns  : $add_pct_mapped_reads_from_samtools
## Arguments:

    my ($arg_href) = @_;

    # Initiate variables
    my $add_pct_mapped_reads_from_samtools = q?'$raw; $map; ?;

    # Remove newline
    $add_pct_mapped_reads_from_samtools .= q?chomp $_; ?;

    # Always relay incoming line to stdout
    $add_pct_mapped_reads_from_samtools .= q?print $_, qq{\n}; ?;

    # Find raw total sequences
    $add_pct_mapped_reads_from_samtools .= q?if ($_=~/raw total sequences:\s+(\d+)/) { ?;

    # Assign raw total sequence
    $add_pct_mapped_reads_from_samtools .= q?$raw = $1; ?;

    # End if
    $add_pct_mapped_reads_from_samtools .= q?} ?;

    # Find reads mapped
    $add_pct_mapped_reads_from_samtools .= q?elsif ($_=~/reads mapped:\s+(\d+)/) { ?;

    # Assign reads mapped
    $add_pct_mapped_reads_from_samtools .= q?$map = $1; ?;

    # Calculate percentage
    $add_pct_mapped_reads_from_samtools .= q?my $percentage = ($map / $raw ) * 100; ?;

    # Write calculation to stdout
    $add_pct_mapped_reads_from_samtools .=
      q?print qq{percentage mapped reads:\t} . $percentage . qq{\n}?;

    # End elsif
    $add_pct_mapped_reads_from_samtools .= q?}?;

    # End oneliner
    $add_pct_mapped_reads_from_samtools .= q?'?;

    return $add_pct_mapped_reads_from_samtools;
}

sub _bcftools_norm_check {

## Function : Return vcf header line matching bcftools norm command
## Returns  : $bcftools_norm_line
## Arguments:

    my ($arg_href) = @_;

    ## Find vcf_key
    my $bcftools_norm_line = q?'if($_=~/\A[#]{2}? . q{bcftools_norm} . q?/) { ?;

    ## Write to stdout
    $bcftools_norm_line .= q?print $_} ?;

    ## If header is finished quit
    $bcftools_norm_line .= q?if($_=~ /\A#CHROM/) {last}'?;

    return $bcftools_norm_line;
}

sub _build_md5sum_check {

## Function : Build a md5sum check file on format: "md5sum  file"
## Returns  : $build_md5sum_check
## Arguments: file_path => Path to file that md5sum was calculated for

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_path;

    my $tmpl = {
        file_path => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Print first md5sum from $md5_file_path
    my $build_md5sum_check = $SINGLE_QUOTE . q?print $F[0].q{? . $SPACE x 2;

    # Print file name in the same line that correspond to md5sum hash
    $build_md5sum_check .= $file_path . q?}? . $SINGLE_QUOTE;

    return $build_md5sum_check;
}

sub _genepred_to_refflat {

## Function : Convert extended genePred format to refFlat format
## Returns  : $genepred_to_refflat
## Arguments:

    my ($arg_href) = @_;

    ## Put gene name first followed by the original first 10 fields
    my $genepred_to_refflat = q?'say STDOUT join qq{\t}, ($F[11], @F[0..9])'?;

    return $genepred_to_refflat;
}

sub _get_dict_contigs {

## Function : Return predifined one liners to get contig names from dict file
## Returns  : $get_dict_contigs
## Arguments:

    my ($arg_href) = @_;

    # Find contig line
    my $get_dict_contigs = q?'if($F[0]=~/^\@SQ/) { ?;

    # Collect contig name
    $get_dict_contigs .= q? if($F[1]=~/SN\:(\S+)/) { ?;

    # Alias capture
    $get_dict_contigs .= q?my $contig_name = $1; ?;

    # Write to STDOUT
    $get_dict_contigs .= q?print $contig_name, q{,};} }'?;

    return $get_dict_contigs;
}

sub _get_fastq_header_v1_4 {

## Function : Return regexp for header elements for fastq file format version 1.4
## Returns  : $get_fastq_header_regexp
## Arguments:

    my ($arg_href) = @_;

    ## Remove newline
    my $get_fastq_header_regexp = q?'chomp; ?;

    ## Get header features
    $get_fastq_header_regexp .=
q?my ($instrument_id, $run_number, $flowcell, $lane, $tile, $x_pos, $y_pos, $direction) = /^(@[^:]*):(\d+):(\w+):(\d+):(\d+):(\d+):(\d+)[\/](\d+)/; ?;

    ## If  we found correct header version
    $get_fastq_header_regexp .= q?if($instrument_id) { ?;

    # Print header features
    $get_fastq_header_regexp .=
q?print join " ", ($instrument_id, $run_number, $flowcell, $lane, $tile, $x_pos, $y_pos, $direction);} ?;

    # Process line one and then exit
    $get_fastq_header_regexp .= q?if($.=1) {last;}' ?;

    return $get_fastq_header_regexp;
}

sub _get_fastq_header_v1_4_interleaved {

## Function : Return reg exp to get read direction for interleaved fastq file format version 1.4
## Returns  : $get_fastq_header_regexp
## Arguments:

    my ($arg_href) = @_;

    ## Remove newline
    my $get_fastq_header_regexp = q?'chomp; ?;

    ## If line one or five
    $get_fastq_header_regexp .= q?if($.==1 or $.==5) { ?;

    ## Get header features
    $get_fastq_header_regexp .=
q?my ($instrument_id, $run_number, $flowcell, $lane, $tile, $x_pos, $y_pos, $direction) = /^(@[^:]*):(\d+):(\w+):(\d+):(\d+):(\d+):(\d+)[\/](\d+)/; ?;

    ## If  we found correct header version
    $get_fastq_header_regexp .= q?if($instrument_id) { ?;

    # Print direction feature
    $get_fastq_header_regexp .= q?print $direction;} } ?;

    # Process to line six and then exit
    $get_fastq_header_regexp .= q?elsif ($.==6) {last;}' ?;

    return $get_fastq_header_regexp;
}

sub _get_fastq_header_v1_8 {

## Function : Return regexp for header elements for fastq interleaved file format version 1.8
## Returns  : $get_fastq_header_regexp
## Arguments:

    my ($arg_href) = @_;

    ## Remove newline
    my $get_fastq_header_regexp = q?'chomp; ?;

    ## Get header features
    $get_fastq_header_regexp .=
q?my ($instrument_id, $run_number, $flowcell, $lane, $tile, $x_pos, $y_pos, $direction, $filtered, $control_bit, $index,) = /^(@[^:]*):(\d+):(\w+):(\d+):(\d+):(\d+):(\d+)\s(\d+):(\w+):(\d+):(\w+)/; ?;

    ## If  we found correct header version
    $get_fastq_header_regexp .= q?if($instrument_id) { ?;

    # Print header features
    $get_fastq_header_regexp .=
q?print join " ", ($instrument_id, $run_number, $flowcell, $lane, $tile, $x_pos, $y_pos, $direction, $filtered, $control_bit, $index); } ?;

    # Process line one and then exit
    $get_fastq_header_regexp .= q?if($.=1) {last;}' ?;

    return $get_fastq_header_regexp;
}

sub _get_fastq_header_v1_8_interleaved {

## Function : Return read direction for interleaved fastq file format version 1.8
## Returns  : $get_fastq_header_regexp
## Arguments:

    my ($arg_href) = @_;

    ## Remove newline
    my $get_fastq_header_regexp = q?'chomp; ?;

    ## If line one or five
    $get_fastq_header_regexp .= q?if($.==1 or $.==5) { ?;

    ## Get header features
    $get_fastq_header_regexp .=
q?my ($instrument_id, $run_number, $flowcell, $lane, $tile, $x_pos, $y_pos, $direction, $filtered, $control_bit, $index,) = /^(@[^:]*):(\d+):(\w+):(\d+):(\d+):(\d+):(\d+)\s(\d+):(\w+):(\d+):(\w+)/; ?;

    ## If  we found correct header version
    $get_fastq_header_regexp .= q?if($instrument_id) { ?;

    # Print direction feature
    $get_fastq_header_regexp .= q?print $direction;} } ?;

    # Process to line six and then exit
    $get_fastq_header_regexp .= q?elsif ($.==6) {last;}' ?;

    return $get_fastq_header_regexp;
}

sub _get_fastq_read_length {

## Function : Return read length from a fastq infile
## Returns  : $read_length_regexp
## Arguments:

    my ($arg_href) = @_;

    ## Prints sequence length and exits

    # Skip header line
    my $read_length_regexp = q?'if ($_!~/@/) {?;

    # Remove newline
    $read_length_regexp .= q?chomp;?;

    # Count chars
    $read_length_regexp .= q?my $seq_length = length;?;

    # Print and exit
    $read_length_regexp .= q?print $seq_length;last;}' ?;

    return $read_length_regexp;
}

sub _get_gene_panel_info {

## Function : Return gene panel data from header
## Returns  : $gene_panel_info_regexp
## Arguments:

    my ($arg_href) = @_;

    ## Create array for storing
    my $gene_panel_info_regexp = q?'BEGIN{ @panels } ?;

    ## Store gene panel to array
    $gene_panel_info_regexp .= q?if (/ \A [#]{2} (gene_panel=[^,]*) /xms ){ push @panels, $1;} ?;

    ## Skip rest if comment
    $gene_panel_info_regexp .= q?elsif (/ \A [^#] /xms) {last;} ?;

    ## Join array to string separated by ':' and print
    $gene_panel_info_regexp .= q?END {print join q{:}, @panels }'?;

    return $gene_panel_info_regexp;
}

sub _get_gene_panel_hgnc_symbols {

## Function : Return HGNC symbols from gene panel
## Returns  : $gene_panel_hgnc_symbols
## Arguments:

    my ($arg_href) = @_;

    # If line starts with gene panel comment
    my $gene_panel_hgnc_symbols = q?'next if (/ \A [#] /xms ); ?;

    # Append ":". Skip rest if it's a comment
    $gene_panel_hgnc_symbols .= q?print $F[4]'?;

    return $gene_panel_hgnc_symbols;
}

sub _reformat_sacct_headers {

## Function : Write individual job line - skip line containing (.batch or .bat+) in the first column
## Returns  : $reformat_sacct_headers
## Arguments: $sacct_header => Sacct header to reformat

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sacct_header;

    my $tmpl = {
        sacct_header => {
            defined     => 1,
            required    => 1,
            store       => \$sacct_header,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Set headers
    my $reformat_sacct_headers = q?'my @headers=(? . $sacct_header . q?); ?;

    # Write header line
    $reformat_sacct_headers .= q?if($. == 1) { print q{#} . join(qq{\t}, @headers), qq{\n} } ?;

    # Write individual job line - skip line containing (.batch or .bat+) in the first column
    $reformat_sacct_headers .=
      q?if ($. >= 3 && $F[0] !~ /( .batch | .bat+ )\b/xms) { print join(qq{\t}, @F), qq{\n} }' ?;

    return $reformat_sacct_headers;
}

sub _get_rrna_transcripts {

## Function : Return rRNA transcripts from gtf file
## Returns  : $rrna_transcripts
## Arguments:

    my ($arg_href) = @_;

    # Print header line
    my $rrna_transcripts = q?'if (/^#/) {print $_} ?;

    # For rRNA, rRNA_pseudogenes or Mt_rRNA
    $rrna_transcripts .= q?elsif ($_ =~ / gene_type \s \"(rRNA|rRNA_pseudogene|Mt_rRNA)\" /nxms)?;

    # Print
    $rrna_transcripts .= q?{print $_}'?;

    return $rrna_transcripts;
}

sub _get_select_contigs_by_col {

## Function : Return contig names from column one of bed file
## Returns  : $get_select_contigs
## Arguments:

    my ($arg_href) = @_;

    # Initilize hash
    my $get_select_contigs = q?'my %contig; ?;

    # Loop per line
    $get_select_contigs .= q?while (<>) { ?;

    # Get contig name
    $get_select_contigs .= q?my ($contig_name) = $_=~/(\w+)/sxm; ?;

    # Skip if on black list
    $get_select_contigs .=
      q?next if($contig_name =~/browser|contig|chromosome|gene_panel|track/); ?;

    # Set contig name in hash
    $get_select_contigs .= q?if($contig_name) { $contig{$contig_name}=undef } } ?;

    # Print contig names to string
    $get_select_contigs .= q?print join ",", sort keys %contig; ?;

    # Exit perl program
    $get_select_contigs .= q?last;' ?;

    return $get_select_contigs;
}

sub _get_vcf_loqusdb_headers {

## Function : Return vcf header line for "Nrcases" and "software line"
## Returns  : $vcf_loqusdb_headers
## Arguments: $id => Header info id

    my ($arg_href) = @_;

    ## If NrCases in line
    my $vcf_loqusdb_headers = q?'if($_=~ /\A [#]{2}NrCases/xsm ?;

    ## or
    $vcf_loqusdb_headers .= q?|| ?;

    ## software in line
    $vcf_loqusdb_headers .= q?$_=~ /\A [#]{2}Software=<ID=loqusdb/xsm) {?;

    ## Print line
    $vcf_loqusdb_headers .= q?print $_}'?;

    return $vcf_loqusdb_headers;
}

sub _get_vcf_header_id_line {

## Function : Return vcf header line matching given id
## Returns  : $vcf_header_line
## Arguments: $id => Header info id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $id;

    my $tmpl = {
        id => {
            defined     => 1,
            required    => 1,
            store       => \$id,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Find vcf_key
    my $vcf_header_line = q?'if($_=~/\A[#]{2}INFO=<ID=? . $id . q?,/) { ?;

    ## Write to stdout
    $vcf_header_line .= q?print $_} ?;

    ## If header is finished quit
    $vcf_header_line .= q?if($_=~ /\A#CHROM/) {last}'?;

    return $vcf_header_line;
}

sub _get_vcf_sample_ids {

## Function : Return sample ids from a vcf file
## Returns  : $get_vcf_sample_ids
## Arguments:

    my ($arg_href) = @_;

    # Find VCF column header line
    my $get_vcf_sample_ids = q?'if ($_ =~ /^#CHROM/ and $F[8] eq q{FORMAT}) {?;

    # Print all sample ids
    $get_vcf_sample_ids .= q?print "@F[9..$#F]"}'?;

    return $get_vcf_sample_ids;
}

sub _remove_decomposed_asterisk_records {

## Function : Remove decomposed '*' records
## Returns  : $remove_star_regexp
## Arguments:

    my ($arg_href) = @_;

    ## As long as the fourth column isn't an asterisk
    my $remove_star_regexp = q?'unless( $F[4] eq q{*} ) { ?;

    ## Print record
    $remove_star_regexp .= q?print $_ }'?;

    return $remove_star_regexp;
}

sub _star_sj_tab_to_bed {

## Function : Convert star splice junction file to bed format
## Returns  : $star_sj_bed
## Arguments:

    my ($arg_href) = @_;

    ## Define motif array
    my $star_sj_bed =
q?'BEGIN { @motifs = qw{ non-canonical GT/AG CT/AC GC/AG CT/GC AT/AC GT/AT }; @strands = qw{ . + - } }; ?;

    ## Only consider junctions with support from uniquely mapped reads
    $star_sj_bed .= q?next if ( $F[6] == 0 ); ?;

    ## Reorder line
    $star_sj_bed .= q?my @elements; ?;

    ## Add contig name
    $star_sj_bed .= q?push @elements, $F[0]; ?;

    ## Add start position
    $star_sj_bed .= q?push @elements, $F[1] - 1; ?;

    ## Add end position
    $star_sj_bed .= q?push @elements, $F[2]; ?;

    ## Build name column
    $star_sj_bed .=
q?my @names = ( qq{motif=$motifs[$F[4]]}, qq{uniquely_mapped=$F[6]}, qq{multi_mapped=$F[7]}, qq{maximum_spliced_alignment_overhang=$F[8]} ); ?;

    ## Join and add name
    $star_sj_bed .= q?push @elements, join q{;}, @names; ?;

    ## Add score column
    $star_sj_bed .= q?push @elements, $F[6]; ?;

    ## Add strand
    $star_sj_bed .= q?push @elements, $strands[$F[3]]; ?;

    ## Print line
    $star_sj_bed .= q?print join qq{\t}, @elements;'?;

    return $star_sj_bed;
}

sub _synonyms_grch37_to_grch38 {

## Function : Return predifined one liners to modify chr prefix from genome version 37 to 38
## Returns  : $modify_chr_prefix
## Arguments:

    my ($arg_href) = @_;

    ## Add "chr" prefix to chromosome name and rename "MT" to "chrM"
    my $modify_chr_prefix = q?'if($_=~s/^MT/chrM/g) {} ?;

    ## Add "chr" prefix to chromosome name
    $modify_chr_prefix .= q?elsif ($_=~s/^([^#])/chr$1/g) {} ?;

    ## Print line
    $modify_chr_prefix .= q?print $_'?;

    return $modify_chr_prefix;
}

sub _synonyms_grch38_to_grch37 {

## Function : Return predifined one liners to modify chr prefix from genome version 37 to 38
## Returns  : $modify_chr_prefix
## Arguments:

    my ($arg_href) = @_;

    ## Remove "chr" prefix from chromosome name and rename "chrM" to "MT"
    my $modify_chr_prefix = q?'if($_=~s/^chrM/MT/g) {} ?;

    ## Remove "chr" prefix from chromosome name
    $modify_chr_prefix .= q?elsif ($_=~s/^chr(.+)/$1/g) {} ?;

    ## Print line
    $modify_chr_prefix .= q?print $_'?;

    return $modify_chr_prefix;
}

sub _write_header_for_contig {

## Function : Write only header line for specific contig and rest of header
## Returns  : $write_header_for_contig
## Arguments: $contig => Contig to create header for

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $contig;

    my $tmpl = {
        contig => {
            defined     => 1,
            required    => 1,
            store       => \$contig,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## If header line
    my $write_header_for_contig = q?'if($_=~/^\@/) { ?;

    ## Write to STDOUT
    $write_header_for_contig .= q?print $_;} ?;

    ## If line begin with specific $contig
    $write_header_for_contig .= q?elsif($_=~/^? . $contig . q?\s+/) { ?;

    ## Write to STDOUT
    $write_header_for_contig .= q?print $_;}'?;

    return $write_header_for_contig;
}

sub _write_contigs_size_file {

## Function : Return predifined one liners to write contig names and length from fai file
## Returns  : $write_contigs_size
## Arguments:

    my ($arg_href) = @_;

    ## Contig name ($F[0]), contig length ($F[1])
    my $write_contigs_size = q?'say STDOUT $F[0] . "\t" . $F[1] '?;

    return $write_contigs_size;
}

1;
