package MIP::Program::Bcftools;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
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
use MIP::Constants qw{ $COMMA $DOUBLE_QUOTE $NEWLINE $PIPE $SPACE };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.22;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      bcftools_annotate
      bcftools_base
      bcftools_call
      bcftools_concat
      bcftools_create_reheader_samples_file
      bcftools_filter
      bcftools_index
      bcftools_merge
      bcftools_mpileup
      bcftools_norm
      bcftools_query
      bcftools_reheader
      bcftools_rename_vcf_samples
      bcftools_roh
      bcftools_stats
      bcftools_view
      bcftools_view_and_index_vcf
    };
}

sub bcftools_annotate {

## Function : Perl wrapper for writing bcftools annotate recipe to $filehandle or return commands array. Based on bcftools 1.9.
## Returns  : @commands
## Arguments: $annotations_file_path  => VCF file or tabix-indexed file path with annotations: CHR\tPOS[\tVALUE]+
##          : $columns_name           => List of columns in the annotation file, e.g. CHROM,POS,REF,ALT,-,INFO/TAG
##          : $filehandle             => Filehandle to write to
##          : $headerfile_path        => File with lines which should be appended to the VCF header
##          : $include                => Include only sites for which the expression is true
##          : $infile_path            => Infile path to read from
##          : $outfile_path           => Outfile path to write to
##          : $output_type            => 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
##          : $regions_ref            => Regions to process {REF}
##          : $remove_ids_ref         => List of annotations to remove
##          : $samples_file_path      => File of samples to annotate
##          : $samples_ref            => Samples to include or exclude if prefixed with "^"
##          : $set_id                 => Set ID column
##          : $stderrfile_path        => Stderr file path to write to {OPTIONAL}
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $annotations_file_path;
    my $columns_name;
    my $filehandle;
    my $include;
    my $infile_path;
    my $headerfile_path;
    my $outfile_path;
    my $regions_ref;
    my $remove_ids_ref;
    my $samples_file_path;
    my $samples_ref;
    my $set_id;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $output_type;

    my $tmpl = {
        annotations_file_path => { store => \$annotations_file_path, strict_type => 1, },
        columns_name          => { store => \$columns_name,          strict_type => 1, },
        filehandle            => { store => \$filehandle, },
        headerfile_path       => { store => \$headerfile_path,       strict_type => 1, },
        include               => { store => \$include,               strict_type => 1, },
        infile_path           => { store => \$infile_path,           strict_type => 1, },
        outfile_path          => { store => \$outfile_path,          strict_type => 1, },
        output_type           => {
            allow       => [qw{ b u z v}],
            default     => q{v},
            store       => \$output_type,
            strict_type => 1,
        },
        regions_ref    => { default => [], store => \$regions_ref, strict_type => 1, },
        remove_ids_ref => {
            default     => [],
            defined     => 1,
            store       => \$remove_ids_ref,
            strict_type => 1,
        },
        samples_file_path => { store => \$samples_file_path, strict_type => 1, },
        samples_ref       => {
            default     => [],
            store       => \$samples_ref,
            strict_type => 1,
        },
        set_id          => { store => \$set_id,          strict_type => 1, },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ bcftools annotate };

    ## Bcftools base args
    @commands = bcftools_base(
        {
            commands_ref      => \@commands,
            outfile_path      => $outfile_path,
            output_type       => $output_type,
            regions_ref       => $regions_ref,
            samples_file_path => $samples_file_path,
            samples_ref       => $samples_ref,
        }
    );

    if ($annotations_file_path) {

        push @commands, q{--annotations} . $SPACE . $annotations_file_path;
    }

    if ($columns_name) {

        push @commands, q{--columns} . $SPACE . $columns_name;
    }

    if ($include) {

        push @commands, q{--include} . $SPACE . $include;
    }

    if ( @{$remove_ids_ref} ) {

        push @commands, q{--remove} . $SPACE . join $COMMA, @{$remove_ids_ref};
    }

    if ($set_id) {

        push @commands, q{--set-id} . $SPACE . $set_id;
    }

    if ($headerfile_path) {

        push @commands, q{--header-lines} . $SPACE . $headerfile_path;
    }

    ## Infile
    if ($infile_path) {

        push @commands, $infile_path;
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

sub bcftools_base {

## Function : Perl wrapper for bcftools base. Based on Bcftools 1.9
## Returns  : @commands
## Arguments: $commands_ref      => List of commands added earlier
##          : $filehandle        => Filehandle to write to
##          : $outfile_path      => Outfile path
##          : $output_type       => 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
##          : $regions_file_path => Bed or vcf file with reions
##          : $regions_ref       => Regions to process {REF}
##          : $samples_file_path => File of samples to annotate
##          : $samples_ref       => Samples to include or exclude if prefixed with "^"
##          : $threads           => Extra compression threds in addition to main thread

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $commands_ref;
    my $filehandle;
    my $outfile_path;
    my $output_type;
    my $regions_file_path;
    my $regions_ref;
    my $samples_file_path;
    my $samples_ref;
    my $threads;

    my $tmpl = {
        commands_ref => {
            default     => [],
            store       => \$commands_ref,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        outfile_path => {
            store       => \$outfile_path,
            strict_type => 1,
        },
        output_type => {
            allow       => [ undef, qw{ b u z v} ],
            store       => \$output_type,
            strict_type => 1,
        },
        regions_file_path => {
            store       => \$regions_file_path,
            strict_type => 1,
        },
        regions_ref => {
            default     => [],
            store       => \$regions_ref,
            strict_type => 1,
        },
        samples_file_path => {
            store       => \$samples_file_path,
            strict_type => 1,
        },
        samples_ref => {
            default     => [],
            store       => \$samples_ref,
            strict_type => 1,
        },
        threads => {
            allow       => qr/\A \d+ \z | undef /xms,
            store       => \$threads,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = @{$commands_ref};

    if ($samples_file_path) {

        push @commands, q{--samples-file} . $SPACE . $samples_file_path;
    }
    if ( @{$samples_ref} ) {

        push @commands, q{--samples} . $SPACE . join $COMMA, @{$samples_ref};
    }
    if ($regions_file_path) {

        push @commands, q{--regions-file} . $SPACE . $regions_file_path;
    }
    if ( @{$regions_ref} ) {

        push @commands, q{--regions} . $SPACE . join $COMMA, @{$regions_ref};
    }
    if ($outfile_path) {

        push @commands, q{--output} . $SPACE . $outfile_path;
    }

    if ($output_type) {

        push @commands, q{--output-type} . $SPACE . $output_type;
    }

    if ($threads) {

        push @commands, q{--threads} . $SPACE . $threads;
    }

    unix_write_to_file(
        {
            commands_ref => \@commands,
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub bcftools_call {

## Function : Perl wrapper for writing bcftools call recipe to $filehandle or return commands array. Based on bcftools 1.6.
## Returns  : @commands
## Arguments: $constrain              => One of: alleles, trio
##          : $filehandle             => Filehandle to write to
##          : $form_fields_ref        => Output format fields {REF}
##          : $infile_path            => Infile path to read from
##          : $multiallelic_caller    => Alternative model for multiallelic and rare-variant calling
##          : $outfile_path           => Outfile path to write to
##          : $output_type            => 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
##          : $regions_ref            => Regions to process {REF}
##          : $samples_file_path      => PED file or a file with an optional column with sex
##          : $samples_ref            => Samples to include or exclude if prefixed with "^"
##          : $stderrfile_path        => Stderr file path to write to {OPTIONAL}
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path
##          : $variants_only          => Output variant sites only

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $constrain;
    my $filehandle;
    my $form_fields_ref;
    my $infile_path;
    my $outfile_path;
    my $regions_ref;
    my $samples_file_path;
    my $samples_ref;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $multiallelic_caller;
    my $output_type;
    my $variants_only;

    my $tmpl = {
        constrain => {
            allow       => [ undef, qw{ alleles trio } ],
            store       => \$constrain,
            strict_type => 1,
        },
        filehandle      => { store => \$filehandle, },
        form_fields_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$form_fields_ref,
            strict_type => 1,
        },
        infile_path         => { store => \$infile_path, strict_type => 1, },
        multiallelic_caller => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$multiallelic_caller,
            strict_type => 1,
        },
        outfile_path => { store => \$outfile_path, strict_type => 1, },
        output_type  => {
            allow       => [qw{ b u z v }],
            default     => q{v},
            store       => \$output_type,
            strict_type => 1,
        },
        regions_ref => { default => [], store => \$regions_ref, strict_type => 1, },
        samples_file_path => { store => \$samples_file_path, strict_type => 1, },
        samples_ref       => {
            default     => [],
            store       => \$samples_ref,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
        variants_only   => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$variants_only,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ bcftools call };

    ## Bcftools base args
    @commands = bcftools_base(
        {
            commands_ref      => \@commands,
            outfile_path      => $outfile_path,
            output_type       => $output_type,
            regions_ref       => $regions_ref,
            samples_file_path => $samples_file_path,
            samples_ref       => $samples_ref,
        }
    );

    if ($multiallelic_caller) {

        push @commands, q{--multiallelic-caller};
    }

    if ( @{$form_fields_ref} ) {

        push @commands, q{--format-fields} . $SPACE . join $COMMA, @{$form_fields_ref};
    }

    if ($variants_only) {

        push @commands, q{--variants-only};
    }

    if ($constrain) {

        push @commands, q{--constrain} . $SPACE . $constrain;
    }

    ## Infile
    if ($infile_path) {

        push @commands, $infile_path;
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

sub bcftools_concat {

## Function : Perl wrapper for writing bcftools concat recipe to $filehandle or return commands array. Based on bcftools 1.6.
## Returns  : @commands
## Arguments: $allow_overlaps         => First coordinate of the next file can precede last record of the current file
##          : $filehandle             => Filehandle to write to
##          : $infile_paths_ref       => Infile path to read from
##          : $outfile_path           => Outfile path to write to
##          : $output_type            => 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
##          : $regions_ref            => Regions to process {REF}
##          : $rm_dups                => Output duplicate records present in multiple files only once: <snps|indels|both|all|none>
##          : $stderrfile_path        => Stderr file path to write to
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile file path to write to
##          : $threads                => Extra compression threds in addition to main thread

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_paths_ref;
    my $outfile_path;
    my $regions_ref;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $threads;

    ## Default(s)
    my $allow_overlaps;
    my $output_type;
    my $rm_dups;

    my $tmpl = {
        allow_overlaps => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$allow_overlaps,
            strict_type => 1,
        },
        filehandle => { store => \$filehandle, },
        infile_paths_ref =>
          { default => [], store => \$infile_paths_ref, strict_type => 1, },
        outfile_path => { store => \$outfile_path, strict_type => 1, },
        output_type  => {
            allow       => [qw{ b u z v }],
            default     => q{v},
            store       => \$output_type,
            strict_type => 1,
        },
        regions_ref => { default => [], store => \$regions_ref, strict_type => 1, },
        rm_dups     => {
            allow       => [qw{ 0 snps indels both all none }],
            default     => q{all},
            store       => \$rm_dups,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
        threads         => {
            allow       => [ undef, qr/ \A \d+ \z /xms ],
            default     => 0,
            store       => \$threads,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ bcftools concat };

    ## Bcftools base args
    @commands = bcftools_base(
        {
            commands_ref => \@commands,
            outfile_path => $outfile_path,
            output_type  => $output_type,
            regions_ref  => $regions_ref,
            threads      => $threads,
        }
    );

    if ($allow_overlaps) {

        push @commands, q{--allow-overlaps};
    }

    if ($rm_dups) {

        push @commands, q{--rm-dups} . $SPACE . $rm_dups;
    }

    if ( @{$infile_paths_ref} ) {

        push @commands, join $SPACE, @{$infile_paths_ref};
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

sub bcftools_filter {

## Function : Perl wrapper for writing bcftools filter recipe to $filehandle or return commands array. Based on bcftools 1.6.
## Returns  : @commands
## Arguments: $filehandle             => Filehandle to write to
##          : $exclude                => Exclude sites for which the expression is true
##          : $filter_mode            => "+": do not replace but add to existing FILTER; "x": reset filters at sites which pass
##          : $include                => Include only sites for which the expression is true
##          : $indel_gap              => Filter clusters of indels separated by <int> or fewer base pairs allowing only one to pass
##          : $infile_path            => Infile paths
##          : $outfile_path           => Outfile path to write to
##          : $output_type            => 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
##          : $regions_ref            => Regions to process {REF}
##          : $samples_file_path      => File of samples to annotate
##          : $samples_ref            => Samples to include or exclude if prefixed with "^"
##          : $soft_filter            => Annotate FILTER column with <string> or unique filter name
##          : $snp_gap                => Filter SNPs within <int> base pairs of an indel
##          : $stdoutfile_path        => Stdoutfile path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $exclude;
    my $filter_mode;
    my $filehandle;
    my $include;
    my $indel_gap;
    my $infile_path;
    my $outfile_path;
    my $regions_ref;
    my $samples_file_path;
    my $samples_ref;
    my $soft_filter;
    my $snp_gap;
    my $stdoutfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;

    ## Default(s)
    my $output_type;

    my $tmpl = {
        filehandle => { store => \$filehandle, },
        exclude    => {
            store       => \$exclude,
            strict_type => 1,
        },
        filter_mode => {
            allow       => [qw{ + x }],
            default     => q{+},
            store       => \$filter_mode,
            strict_type => 1,
        },
        infile_path => {
            store       => \$infile_path,
            strict_type => 1,
        },
        include => {
            store       => \$include,
            strict_type => 1,
        },
        indel_gap => {
            allow       => qr/ \A \d+ \z /sxm,
            store       => \$indel_gap,
            strict_type => 1,
        },
        outfile_path => {
            store       => \$outfile_path,
            strict_type => 1,
        },
        output_type => {
            allow       => [qw{ b u z v}],
            default     => q{v},
            store       => \$output_type,
            strict_type => 1,
        },
        regions_ref => {
            default     => [],
            store       => \$regions_ref,
            strict_type => 1,
        },
        samples_file_path => {
            store       => \$samples_file_path,
            strict_type => 1,
        },
        samples_ref => {
            default     => [],
            store       => \$samples_ref,
            strict_type => 1,
        },
        soft_filter => {
            store       => \$soft_filter,
            strict_type => 1,
        },
        snp_gap => {
            allow       => qr/ \A \d+ \z /sxm,
            store       => \$snp_gap,
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ bcftools filter };

    ## Bcftools base args
    @commands = bcftools_base(
        {
            commands_ref      => \@commands,
            outfile_path      => $outfile_path,
            output_type       => $output_type,
            regions_ref       => $regions_ref,
            samples_file_path => $samples_file_path,
            samples_ref       => $samples_ref,
        }
    );

    if ($exclude) {

        push @commands, q{--exclude} . $SPACE . $exclude;
    }

    if ($filter_mode) {

        push @commands, q{--mode} . $SPACE . $filter_mode;
    }

    if ($include) {

        push @commands, q{--include} . $SPACE . $include;
    }

    if ($indel_gap) {

        push @commands, q{--IndelGap} . $SPACE . $indel_gap;
    }

    ## Infile
    if ($infile_path) {

        push @commands, $infile_path;
    }

    if ($soft_filter) {

        push @commands, q{--soft-filter} . $SPACE . $soft_filter;
    }

    if ($snp_gap) {

        push @commands, q{--SnpGap} . $SPACE . $snp_gap;
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

sub bcftools_index {

## Function : Perl wrapper for writing bcftools index recipe to $filehandle or return commands array. Based on bcftools 1.6.
## Returns  : @commands
## Arguments: $filehandle             => Filehandle to write to
##          : $infile_path            => Infile path to read from
##          : $outfile_path           => Outfile path to write to
##          : $output_type            => 'csi' generate CSI-format index, 'tbi' generate TBI-format index
##          : $regions_ref            => Regions to process {REF}
##          : $samples_file_path      => File of samples to annotate
##          : $samples_ref            => Samples to include or exclude if prefixed with "^"
##          : $stderrfile_path        => Stderr file path to write to {OPTIONAL}
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $outfile_path;
    my $regions_ref;
    my $samples_file_path;
    my $samples_ref;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $output_type;

    my $tmpl = {
        filehandle  => { store => \$filehandle, },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        outfile_path => { store => \$outfile_path, strict_type => 1, },
        output_type  => {
            allow       => [qw{ csi tbi }],
            default     => q{csi},
            store       => \$output_type,
            strict_type => 1,
        },
        regions_ref => { default => [], store => \$regions_ref, strict_type => 1, },
        samples_file_path => { store => \$samples_file_path, strict_type => 1, },
        samples_ref       => {
            default     => [],
            store       => \$samples_ref,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ bcftools index };

    ## Bcftools base args
    @commands = bcftools_base(
        {
            commands_ref      => \@commands,
            outfile_path      => $outfile_path,
            regions_ref       => $regions_ref,
            samples_file_path => $samples_file_path,
            samples_ref       => $samples_ref,
        }
    );

    # Special case: 'csi' or 'tbi'
    if ($output_type) {

        #Specify output type
        push @commands, q{--} . $output_type;
    }

    push @commands, $infile_path;

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

sub bcftools_merge {

## Function : Perl wrapper for writing bcftools merge recipe to $filehandle or return commands array. Based on bcftools 1.6.
## Returns  : @commands
## Arguments: $filehandle             => Filehandle to write to
##          : $infile_paths_ref       => Infile path to read from
##          : $outfile_path           => Outfile path to write to
##          : $output_type            => 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
##          : $regions_ref            => Regions to process {REF}
##          : $samples_file_path      => File of samples to annotate
##          : $samples_ref            => Samples to include or exclude if prefixed with "^"
##          : $stderrfile_path        => Stderr file path to write to
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile file path to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_paths_ref;
    my $outfile_path;
    my $regions_ref;
    my $samples_file_path;
    my $samples_ref;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $output_type;

    my $tmpl = {
        filehandle => { store => \$filehandle, },
        infile_paths_ref =>
          { default => [], store => \$infile_paths_ref, strict_type => 1, },
        outfile_path => { store => \$outfile_path, strict_type => 1, },
        output_type  => {
            allow       => [qw{ b u z v}],
            default     => q{v},
            store       => \$output_type,
            strict_type => 1,
        },
        regions_ref => { default => [], store => \$regions_ref, strict_type => 1, },
        samples_file_path => { store => \$samples_file_path, strict_type => 1, },
        samples_ref       => {
            default     => [],
            store       => \$samples_ref,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ bcftools merge };

    ## Bcftools base args
    @commands = bcftools_base(
        {
            commands_ref      => \@commands,
            outfile_path      => $outfile_path,
            output_type       => $output_type,
            regions_ref       => $regions_ref,
            samples_file_path => $samples_file_path,
            samples_ref       => $samples_ref,
        }
    );

    if ( @{$infile_paths_ref} ) {

        push @commands, join $SPACE, @{$infile_paths_ref};
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

sub bcftools_mpileup {

## Function : Perl wrapper for writing bcftools mpileup recipe to $filehandle. Based on bcftools 1.6 (using htslib 1.6).
## Returns  : @commands
##          : $adjust_mq                        => Adjust mapping quality
##          : $filehandle                       => Sbatch filehandle to write to
##          : $infile_paths_ref                 => Infile paths {REF}
##          : $outfile_path                     => Outfile path
##          : $output_tags_ref                  => Optional tags to output {REF}
##          : $output_type                      => 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
##          : $per_sample_increased_sensitivity => Apply -m and -F per-sample for increased sensitivity
##          : $referencefile_path               => Reference sequence file
##          : $regions_ref                      => Regions to process {REF}
##          : $samples_file_path                => File of samples to annotate
##          : $samples_ref                      => Samples to include or exclude if prefixed with "^"
##          : $stderrfile_path                  => Stderrfile path
##          : $stderrfile_path_append           => Stderrfile path append
##          : $stdoutfile_path                  => Stdoutfile file path to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_paths_ref;
    my $outfile_path;
    my $output_tags_ref;
    my $referencefile_path;
    my $regions_ref;
    my $samples_file_path;
    my $samples_ref;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $adjust_mq;
    my $per_sample_increased_sensitivity;
    my $output_type;

    ## Constants
    Readonly my $ADJUST_MAPPING_QUALITY => 50;

    my $tmpl = {
        adjust_mq => {
            allow       => qr/ \A \d+ \z /sxm,
            default     => $ADJUST_MAPPING_QUALITY,
            store       => \$adjust_mq,
            strict_type => 1,
        },
        filehandle       => { store => \$filehandle, },
        infile_paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$infile_paths_ref,
            strict_type => 1,
        },
        outfile_path    => { store => \$outfile_path, strict_type => 1, },
        output_tags_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$output_tags_ref,
            strict_type => 1,
        },
        output_type => {
            allow       => [qw{ b u z v}],
            default     => q{b},
            store       => \$output_type,
            strict_type => 1,
        },
        per_sample_increased_sensitivity => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$per_sample_increased_sensitivity,
            strict_type => 1,
        },
        referencefile_path => {
            defined     => 1,
            required    => 1,
            store       => \$referencefile_path,
            strict_type => 1,
        },
        regions_ref => { default => [], store => \$regions_ref, strict_type => 1, },
        samples_file_path => { store => \$samples_file_path, strict_type => 1, },
        samples_ref       => {
            default     => [],
            store       => \$samples_ref,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ bcftools mpileup };

    ## Bcftools base args
    @commands = bcftools_base(
        {
            commands_ref      => \@commands,
            regions_ref       => $regions_ref,
            outfile_path      => $outfile_path,
            output_type       => $output_type,
            samples_file_path => $samples_file_path,
            samples_ref       => $samples_ref,
        }
    );

    push @commands, q{--adjust-MQ} . $SPACE . $adjust_mq;

    if ($per_sample_increased_sensitivity) {

        push @commands, q{--per-sample-mF};
    }

    if ( @{$output_tags_ref} ) {

        push @commands, q{--annotate} . $SPACE . join $COMMA, @{$output_tags_ref};
    }

    push @commands, q{--fasta-ref} . $SPACE . $referencefile_path;

    push @commands, join $SPACE, @{$infile_paths_ref};

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

sub bcftools_norm {

## Function : Perl wrapper for writing bcftools norm recipe to $filehandle or return commands array. Based on bcftools 1.6.
## Returns  : @commands
## Arguments: $filehandle             => Filehandle to write to
##          : $infile_path            => Infile path to read from
##          : $multiallelic           => To split/join multiallelic calls or not
##          : $multiallelic_type      => Type of multiallelic to split/join {OPTIONAL}
##          : $outfile_path           => Outfile path to write to
##          : $output_type            => 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
##          : $reference_path         => Human genome reference path
##          : $regions_ref            => Regions to process {REF}
##          : $samples_file_path      => File of samples to annotate
##          : $samples_ref            => Samples to include or exclude if prefixed with "^"
##          : $stderrfile_path        => Stderr file path to write to {OPTIONAL}
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $multiallelic;
    my $outfile_path;
    my $reference_path;
    my $regions_ref;
    my $samples_file_path;
    my $samples_ref;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $multiallelic_type;
    my $output_type;

    my $tmpl = {
        filehandle   => { store => \$filehandle, },
        infile_path  => { store => \$infile_path, strict_type => 1, },
        multiallelic => {
            allow       => [qw{ + - }],
            store       => \$multiallelic,
            strict_type => 1,
        },
        multiallelic_type => {
            allow       => [qw{ snps indels both any }],
            default     => q{both},
            store       => \$multiallelic_type,
            strict_type => 1,
        },
        outfile_path => {
            store       => \$outfile_path,
            strict_type => 1,
        },
        output_type => {
            allow       => [qw{ b u z v }],
            default     => q{v},
            store       => \$output_type,
            strict_type => 1,
        },
        reference_path => {
            defined     => 1,
            required    => 1,
            store       => \$reference_path,
            strict_type => 1,
        },
        regions_ref => { default => [], store => \$regions_ref, strict_type => 1, },
        samples_file_path => { store => \$samples_file_path, strict_type => 1, },
        samples_ref       => {
            default     => [],
            store       => \$samples_ref,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ bcftools norm };

    ## Bcftools base args
    @commands = bcftools_base(
        {
            commands_ref      => \@commands,
            regions_ref       => $regions_ref,
            outfile_path      => $outfile_path,
            output_type       => $output_type,
            samples_file_path => $samples_file_path,
            samples_ref       => $samples_ref,
        }
    );

    if ($multiallelic) {

        push @commands, q{--multiallelics} . $SPACE . $multiallelic . $multiallelic_type;
    }

    if ($reference_path) {

        push @commands, q{--fasta-ref} . $SPACE . $reference_path;
    }

    if ($infile_path) {

        push @commands, $infile_path;
    }

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
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

sub bcftools_query {

## Function : Perl wrapper for writing bcftools query recipe to $filehandle or return commands array. Based on bcftools 1.10.
## Returns  : @commands
## Arguments: $exclude                => Exclude sites for which the expression is true
##          : $filehandle             => Filehandle to write to
##          : $format                 => Format
##          : $include                => Include only sites for which the expression is true
##          : $infile_paths_ref       => Infile paths to read from {REF}
##          : $outfile_path           => Outfile path to write to
##          : $print_header           => Print header
##          : $regions_ref            => Regions to process {REF}
##          : $samples_file_path      => File of samples to annotate
##          : $samples_ref            => Samples to include or exclude if prefixed with "^"
##          : $stderrfile_path        => Stderr file path to write to
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile file path to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $exclude;
    my $filehandle;
    my $format;
    my $include;
    my $infile_paths_ref;
    my $outfile_path;
    my $print_header;
    my $regions_ref;
    my $samples_file_path;
    my $samples_ref;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $output_type;

    my $tmpl = {
        exclude    => { store => \$exclude, strict_type => 1, },
        filehandle => { store => \$filehandle, },
        format => {
            store       => \$format,
            strict_type => 1,
        },
        include          => { store => \$include, strict_type => 1, },
        infile_paths_ref => {
            default     => [],
            store       => \$infile_paths_ref,
            strict_type => 1,
        },
        outfile_path => { store => \$outfile_path, strict_type => 1, },
        print_header => {
            allow       => [ undef, 0, 1 ],
            store       => \$print_header,
            strict_type => 1,
        },
        regions_ref => { default => [], store => \$regions_ref, strict_type => 1, },
        samples_file_path => { store => \$samples_file_path, strict_type => 1, },
        samples_ref       => {
            default     => [],
            store       => \$samples_ref,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ bcftools query };

    ## Bcftools base args
    @commands = bcftools_base(
        {
            commands_ref      => \@commands,
            regions_ref       => $regions_ref,
            samples_file_path => $samples_file_path,
            samples_ref       => $samples_ref,
        }
    );

    if ($exclude) {

        push @commands, q{--exclude} . $SPACE . $exclude;
    }
    if ($format) {

        push @commands, q{--format} . $SPACE . $format;
    }
    if ($include) {

        push @commands, q{--include} . $SPACE . $include;
    }
    if ($outfile_path) {

        push @commands, q{--output} . $SPACE . $outfile_path;
    }
    if ($print_header) {

        push @commands, q{--print-header};
    }
    if ($infile_paths_ref) {

        push @commands, join $SPACE, @{$infile_paths_ref};
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

sub bcftools_reheader {

## Function : Perl wrapper for writing bcftools reheader recipe to already open $filehandle or return commands array. Based on bcftools 1.6.
## Returns  : @commands
## Arguments: $filehandle             => Filehandle to write to
##          : $infile_path            => Infile path
##          : $outfile_path           => Outfile path
##          : $regions_ref            => Regions to process {REF}
##          : $samples_file_path      => File of samples to annotate
##          : $samples_ref            => Samples to include or exclude if prefixed with "^"
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $outfile_path;
    my $regions_ref;
    my $samples_file_path;
    my $samples_ref;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        filehandle        => { store   => \$filehandle, },
        infile_path       => { store   => \$infile_path, strict_type => 1, },
        outfile_path      => { store   => \$outfile_path, strict_type => 1, },
        regions_ref       => { default => [], store => \$regions_ref, strict_type => 1, },
        samples_file_path => { store   => \$samples_file_path, strict_type => 1, },
        samples_ref => {
            default     => [],
            store       => \$samples_ref,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ bcftools reheader };

    ## Bcftools base args
    @commands = bcftools_base(
        {
            commands_ref => \@commands,
            regions_ref  => $regions_ref,
            outfile_path => $outfile_path,
        }
    );

    if ($samples_file_path) {

        push @commands, q{--samples} . $SPACE . $samples_file_path;
    }

    if ($infile_path) {

        push @commands, $infile_path;
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

sub bcftools_rename_vcf_samples {

## Function : Rename vcf samples. The samples array will replace the sample names in the same order as supplied.
## Returns  :
## Arguments: $create_sample_file  => Create sample file for bcftools reheader
##          : $filehandle          => Filehandle to write to
##          : $index               => Generate index of reformated file
##          : $index_type          => Type of index
##          : $infile              => Vcf infile to rename samples for
##          : $outfile_path_prefix => Out file path no file_ending {Optional}
##          : $output_type         => Output type
##          : $sample_ids_ref      => Samples to rename in the same order as in the vcf {REF}
##          : $temp_directory      => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $create_sample_file;
    my $filehandle;
    my $infile;
    my $outfile_path_prefix;
    my $sample_ids_ref;
    my $temp_directory;

    ## Default(s)
    my $index;
    my $index_type;
    my $output_type;

    my $tmpl = {
        create_sample_file => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$create_sample_file,
            strict_type => 1,
        },
        filehandle => { defined => 1, required => 1, store => \$filehandle, },
        index      => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$index,
            strict_type => 1,
        },
        index_type => {
            allow       => [ undef, qw{ csi tbi } ],
            default     => q{csi},
            store       => \$index_type,
            strict_type => 1,
        },
        infile => { defined => 1, required => 1, store => \$infile, strict_type => 1, },
        outfile_path_prefix => { store => \$outfile_path_prefix, strict_type => 1, },
        output_type         => {
            allow       => [qw{ b u z v }],
            default     => q{v},
            store       => \$output_type,
            strict_type => 1,
        },
        sample_ids_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$sample_ids_ref,
            strict_type => 1,
        },
        temp_directory => {
            defined     => 1,
            required    => 1,
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    if ($create_sample_file) {

        bcftools_create_reheader_samples_file(
            {
                filehandle     => $filehandle,
                sample_ids_ref => $sample_ids_ref,
                temp_directory => $temp_directory,
            }
        );
    }

    ## Rename samples in VCF
    say {$filehandle} q{## Rename sample(s) names in VCF file};
    bcftools_reheader(
        {
            filehandle        => $filehandle,
            infile_path       => $infile,
            samples_file_path => catfile( $temp_directory, q{sample_name.txt} ),
        }
    );
    ## Pipe
    print {$filehandle} $PIPE . $SPACE;

    bcftools_view_and_index_vcf(
        {
            filehandle          => $filehandle,
            index               => $index,
            index_type          => $index_type,
            infile_path         => q{-},
            outfile_path_prefix => $outfile_path_prefix,
            output_type         => $output_type,
        }
    );
    return;
}

sub bcftools_create_reheader_samples_file {

## Function : Create reheader samples file
## Returns  :
## Arguments: $filehandle     => Filehandle to write to
##          : $sample_ids_ref => Samples to rename in the same order as in the vcf {REF}
##          : $temp_directory => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $sample_ids_ref;
    my $temp_directory;

    my $tmpl = {
        filehandle     => { defined => 1, required => 1, store => \$filehandle, },
        sample_ids_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$sample_ids_ref,
            strict_type => 1,
        },
        temp_directory => {
            defined     => 1,
            required    => 1,
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Program::Gnu::Coreutils qw{ gnu_printf };

    ## Create new sample names file
    say {$filehandle} q{## Create new sample(s) names file};

    ## Get parameters
    my $format_string = $DOUBLE_QUOTE;
  SAMPLE_ID:
    foreach my $sample_id ( @{$sample_ids_ref} ) {

        $format_string .= $sample_id . q{\n};
    }
    $format_string .= $DOUBLE_QUOTE;
    gnu_printf(
        {
            filehandle      => $filehandle,
            format_string   => $format_string,
            stdoutfile_path => catfile( $temp_directory, q{sample_name.txt} ),
        }
    );
    say {$filehandle} $NEWLINE;
    return;
}

sub bcftools_roh {

## Function : Perl wrapper for writing bcftools roh recipe to $filehandle or return commands array. Based on bcftools 1.6.
## Returns  : @commands
## Arguments: $af_file_path           => Read allele frequencies from file (CHR\tPOS\tREF,ALT\tAF)
##          : $af_tag                 => Use tag in info field as allelle frequency
##          : $filehandle             => Filehandle to write to
##          : $infile_path            => Infile path to read from
##          : $outfile_path           => Outfile path to write to
##          : $output_type            => 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
##          : $regions_ref            => Regions to process {REF}
##          : $samples_file_path      => File of samples to annotate
##          : $samples_ref            => Samples to include or exclude if prefixed with "^"
##          : $skip_indels            => Skip indels as their genotypes are enriched for errors
##          : $stderrfile_path        => Stderr file path to write to
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile file path to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $af_file_path;
    my $af_tag;
    my $filehandle;
    my $infile_path;
    my $outfile_path;
    my $regions_ref;
    my $samples_file_path;
    my $samples_ref;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $skip_indels;
    my $output_type;

    my $tmpl = {
        af_file_path => { store => \$af_file_path, strict_type => 1, },
        af_tag       => {
            store       => \$af_tag,
            strict_type => 1,
        },
        filehandle  => { store => \$filehandle, },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        outfile_path => { store => \$outfile_path, strict_type => 1, },
        output_type  => {
            allow       => [qw{ b u z v}],
            default     => q{v},
            store       => \$output_type,
            strict_type => 1,
        },
        regions_ref => { default => [], store => \$regions_ref, strict_type => 1, },
        samples_file_path => { store => \$samples_file_path, strict_type => 1, },
        samples_ref       => {
            default     => [],
            store       => \$samples_ref,
            strict_type => 1,
        },
        samples_ref => {
            default     => [],
            store       => \$samples_ref,
            strict_type => 1,
        },
        skip_indels => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$skip_indels,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ bcftools roh };

    ## Bcftools base args
    @commands = bcftools_base(
        {
            commands_ref      => \@commands,
            regions_ref       => $regions_ref,
            output_type       => $output_type,
            outfile_path      => $outfile_path,
            samples_file_path => $samples_file_path,
            samples_ref       => $samples_ref,
        }
    );

    if ($af_file_path) {

        push @commands, q{--AF-file} . $SPACE . $af_file_path;
    }

    if ($af_tag) {

        push @commands, q{--AF-tag} . $SPACE . $af_tag;
    }

    if ($skip_indels) {

        push @commands, q{--skip-indels};
    }

    if ($infile_path) {

        push @commands, $infile_path;
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

sub bcftools_stats {

## Function : Perl wrapper for writing bcftools stats recipe to already open $filehandle or return commands array. Based on bcftools 1.6.
## Returns  : @commands
## Arguments: $infile_path            => Infile path
##          : $outfile_path           => Outfile path
##          : $output_type            => 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
##          : $regions_ref            => Regions to process {REF}
##          : $samples_file_path      => File of samples to annotate
##          : $samples_ref            => Samples to include or exclude if prefixed with "^"
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path
##          : $filehandle             => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_path;
    my $outfile_path;
    my $regions_ref;
    my $samples_file_path;
    my $samples_ref;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $filehandle;

    ## Default(s)
    my $output_type;

    my $tmpl = {
        filehandle   => { store => \$filehandle, },
        infile_path  => { store => \$infile_path, strict_type => 1, },
        outfile_path => { store => \$outfile_path, strict_type => 1, },
        output_type => {
            allow       => [qw{ b u z v}],
            default     => q{v},
            store       => \$output_type,
            strict_type => 1,
        },
        regions_ref => { default => [], store => \$regions_ref, strict_type => 1, },
        samples_file_path => { store => \$samples_file_path, strict_type => 1, },
        samples_ref       => {
            default     => [],
            store       => \$samples_ref,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ bcftools stats };

    ## Bcftools base args
    @commands = bcftools_base(
        {
            commands_ref      => \@commands,
            regions_ref       => $regions_ref,
            samples_file_path => $samples_file_path,
            samples_ref       => $samples_ref,
        }
    );

    if ($infile_path) {

        push @commands, $infile_path;
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

sub bcftools_view {

## Function : Perl wrapper for writing bcftools view recipe to $filehandle or return commands array. Based on bcftools 1.6.
## Returns  : @commands
## Arguments: $apply_filters_ref      => Require at least one of the listed FILTER strings
##          : $exclude_types_ref      => Exclude comma-separated list of variant types: snps,indels,mnps,other
##          : $exclude                => Exclude sites for which the expression is true
##          : $filehandle             => Filehandle to write to
##          : $genotype               => Genotype to include (hom|het|miss). Prefix with "^" for exclude
##          : $header_only            => Header only
##          : $include                => Include only sites for which the expression is true
##          : $infile_path            => Infile path to read from
##          : $max_alleles            => Max alleles listed in REF and ALT columns
##          : $min_ac                 => Min allele count for sites to be printed
##          : $min_alleles            => Min alleles listed in REF and ALT columns
##          : $no_header              => No header
##          : $outfile_path           => Outfile path to write to
##          : $output_type            => 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
##          : $regions_file_path      => Bed or vcf file with reions
##          : $regions_ref            => Regions to process {REF}
##          : $samples_file_path      => File of samples to annotate
##          : $samples_ref            => Samples to include or exclude if prefixed with "^"
##          : $stderrfile_path        => Stderr file path to write to
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile file path to write to
##          : $threads                => Number of threads to use
##          : $types                  => Comma separated variant types to include (snps|indels|mnps|other), based on based on REF,ALT

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $apply_filters_ref;
    my $exclude_types_ref;
    my $exclude;
    my $filehandle;
    my $genotype;
    my $header_only;
    my $include;
    my $infile_path;
    my $max_alleles;
    my $min_ac;
    my $min_alleles;
    my $no_header;
    my $outfile_path;
    my $regions_file_path;
    my $regions_ref;
    my $samples_file_path;
    my $samples_ref;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $threads;
    my $types;

    ## Default(s)
    my $output_type;

    my $tmpl = {
        filehandle => { store => \$filehandle, },
        apply_filters_ref =>
          { default => [], store => \$apply_filters_ref, strict_type => 1, },
        exclude_types_ref =>
          { default => [], store => \$exclude_types_ref, strict_type => 1, },
        exclude  => { store => \$exclude, strict_type => 1, },
        genotype => {
            store       => \$genotype,
            strict_type => 1,
        },
        header_only => {
            allow       => [ undef, 0, 1 ],
            store       => \$header_only,
            strict_type => 1,
        },
        include     => { store => \$include,     strict_type => 1, },
        infile_path => { store => \$infile_path, strict_type => 1, },
        max_alleles => {
            allow       => [ undef, qr/ \A \d+ \z /xms ],
            store       => \$max_alleles,
            strict_type => 1,
        },
        min_ac => {
            store       => \$min_ac,
            strict_type => 1,
        },
        min_alleles => {
            allow       => [ undef, qr/ \A \d+ \z /xms ],
            store       => \$min_alleles,
            strict_type => 1,
        },
        no_header => {
            allow       => [ undef, 0, 1 ],
            store       => \$no_header,
            strict_type => 1,
        },
        outfile_path => { store => \$outfile_path, strict_type => 1, },
        output_type  => {
            allow       => [qw{ b u z v}],
            default     => q{v},
            store       => \$output_type,
            strict_type => 1,
        },
        regions_file_path => { store => \$regions_file_path, strict_type => 1, },
        regions_ref       => { default => [], store => \$regions_ref, strict_type => 1, },
        samples_file_path => { store => \$samples_file_path, strict_type => 1, },
        samples_ref       => {
            default     => [],
            store       => \$samples_ref,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
        threads         => {
            allow       => qr/ \A \d+ \z /xms,
            store       => \$threads,
            strict_type => 1,
        },
        types => {
            store       => \$types,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ bcftools view };

    ## Bcftools base args
    @commands = bcftools_base(
        {
            commands_ref      => \@commands,
            output_type       => $output_type,
            regions_file_path => $regions_file_path,
            regions_ref       => $regions_ref,
            samples_file_path => $samples_file_path,
            samples_ref       => $samples_ref,
        }
    );

    if ( @{$apply_filters_ref} ) {

        push @commands, q{--apply-filters} . $SPACE . join $COMMA, @{$apply_filters_ref};
    }

    if ( @{$exclude_types_ref} ) {

        push @commands, q{--exclude-types} . $SPACE . join $COMMA, @{$exclude_types_ref};
    }

    if ($exclude) {

        push @commands, q{--exclude} . $SPACE . $exclude;
    }

    if ($genotype) {

        push @commands, q{--genotype} . $SPACE . $genotype;
    }

    if ($header_only) {

        push @commands, q{--header-only};
    }

    if ($include) {

        push @commands, q{--include} . $SPACE . $include;
    }

    if ($max_alleles) {

        push @commands, q{--max-alleles} . $SPACE . $max_alleles;
    }

    if ($min_ac) {

        push @commands, q{--min-ac} . $SPACE . $min_ac;
    }

    if ($min_alleles) {

        push @commands, q{--min-alleles} . $SPACE . $min_alleles;
    }

    if ($no_header) {

        push @commands, q{--no-header};
    }

    if ($outfile_path) {

        push @commands, q{--output-file} . $SPACE . $outfile_path;
    }

    if ($types) {

        push @commands, q{--types} . $SPACE . $types;
    }

    if ($threads) {

        push @commands, q{--threads} . $SPACE . $threads;
    }

    ## Infile
    if ($infile_path) {

        push @commands, $infile_path;
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

sub bcftools_view_and_index_vcf {

## Function : View variant calling file and index.
## Returns  :
## Arguments: $filehandle          => SBATCH script filehandle to print to
##          : $index               => Generate index of reformated file
##          : $index_type          => Type of index
##          : $infile_path         => Path to infile to compress and index
##          : $outfile_path_prefix => Out file path no file_ending {Optional}
##          : $output_type         => 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
##          : $threads             => Number of threads to use

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $outfile_path_prefix;

    ## Default(s)
    my $index;
    my $index_type;
    my $output_type;
    my $threads;

    my $tmpl = {
        filehandle  => { defined => 1, required => 1, store => \$filehandle, },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        index => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$index,
            strict_type => 1,
        },
        index_type => {
            allow       => [ undef, qw{ csi tbi } ],
            default     => q{csi},
            store       => \$index_type,
            strict_type => 1,
        },
        outfile_path_prefix => { store => \$outfile_path_prefix, strict_type => 1, },
        output_type         => {
            allow       => [qw{ b u z v }],
            default     => q{v},
            store       => \$output_type,
            strict_type => 1,
        },
        threads => {
            allow       => qr/ \A \d+ \z /xms,
            default     => 1,
            store       => \$threads,
            strict_type => 1,
        },

    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $outfile_path;
    my %output_type_ending = (
        b => q{.bcf},
        u => q{.bcf},
        v => q{.vcf},
        z => q{.vcf.gz},
    );

    if ( defined $outfile_path_prefix ) {

        $outfile_path = $outfile_path_prefix . $output_type_ending{$output_type};
    }

    bcftools_view(
        {
            filehandle   => $filehandle,
            infile_path  => $infile_path,
            outfile_path => $outfile_path,
            output_type  => $output_type,
            threads      => $threads,
        }
    );
    say {$filehandle} $NEWLINE;

    if ($index) {

        say {$filehandle} q{## Index};

        bcftools_index(
            {
                filehandle  => $filehandle,
                infile_path => $outfile_path,
                output_type => $index_type,
            }
        );
        say {$filehandle} $NEWLINE;
    }
    return;
}

1;
