package MIP::Environment::Executable;

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
use MIP::Constants qw{ $EMPTY_STR $LOG_NAME $PIPE $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.12;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      build_binary_version_cmd
      get_binaries_versions
      get_executable
      get_binary_version
      set_executable_container_cmd
    };
}

sub build_binary_version_cmd {

## Function : Build binary version commands
## Returns  : @version_cmds
## Arguments: $binary_path    => Executables (binary) file path
##          : $version_cmd    => Version command line option
##          : $version_regexp => Version reg exp to get version from system call

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $binary_path;
    my $version_cmd;
    my $version_regexp;

    my $tmpl = {
        binary_path => {
            defined     => 1,
            required    => 1,
            store       => \$binary_path,
            strict_type => 1,
        },
        version_cmd => {
            store       => \$version_cmd,
            strict_type => 1,
        },
        version_regexp => {
            defined     => 1,
            required    => 1,
            store       => \$version_regexp,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Language::Perl qw{ perl_nae_oneliners };

    ## Get perl wrapper around regexp
    my @perl_commands = perl_nae_oneliners(
        {
            oneliner_cmd => $version_regexp,
        }
    );

    my @version_cmds = ($binary_path);

    ## Allow for binary without any option to get version
    if ($version_cmd) {

        push @version_cmds, $version_cmd;
    }
    push @version_cmds, ( $PIPE, @perl_commands );
    return @version_cmds;
}

sub get_binaries_versions {

## Function : Get executables/binaries versions
## Returns  : %binary_version
## Arguments: $binary_info_href => Binary_Info_Href object

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $binary_info_href;

    my $tmpl = {
        binary_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$binary_info_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my %binary;
    my %executable = get_executable( {} );

  BINARY:
    while ( my ( $binary, $binary_path ) = each %{$binary_info_href} ) {

        ## No information on how to get version for this binary - skip
        next BINARY if ( not exists $executable{$binary} );

        my $binary_version = get_binary_version(
            {
                binary      => $binary,
                binary_path => $binary_path,
            }
        );

        ## Set binary version and path
        $binary{$binary} = (
            {
                path    => $binary_path,
                version => $binary_version,
            }
        );
    }
    return %binary;
}

sub get_binary_version {

## Function : Get version for executable/binary
## Returns  : $binary_version
## Arguments: $binary      => Binary to get version of
##          : $binary_path => Path to binary

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $binary;
    my $binary_path;

    my $tmpl = {
        binary => {
            defined     => 1,
            required    => 1,
            store       => \$binary,
            strict_type => 1,
        },
        binary_path => {
            defined     => 1,
            required    => 1,
            store       => \$binary_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Environment::Child_process qw{ child_process };

## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my %executable = get_executable( { executable_name => $binary, } );

    ## Get version command
    my @version_cmds = build_binary_version_cmd(
        {
            binary_path    => $binary_path,
            version_cmd    => $executable{version_cmd},
            version_regexp => $executable{version_regexp},
        }
    );

    ## Call binary and parse output to generate version
    my %process_return = child_process(
        {
            commands_ref => \@version_cmds,
            process_type => q{open3},
        }
    );

    my $binary_version = $process_return{stdouts_ref}[0];

    if ( not $binary_version ) {

        $log->warn(qq{Could not find version for binary: $binary});
    }
    return $binary_version;
}

sub get_executable {

## Function : Define the executable features and return them
## Returns  : %{ $executable{$executable_name} } or %executable
## Arguments: $executable_name => Executable name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $executable_name;

    my $tmpl = {
        executable_name => {
            store       => \$executable_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my %executable = (
        arriba => {
            version_cmd => q{-h},
            version_regexp =>
q?'my ($version) = /\AVersion:\s(\S+)/xms; if($version) {print $version;last;}'?,
        },
        q{bam2wig.py} => {
            version_cmd => q{--version},
            version_regexp =>
q?'my ($version) = /bam2wig.py\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        q{bam_stat.py} => {
            version_cmd => q{--version},
            version_regexp =>
q?'my ($version) = /bam_stat.py\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        bcftools => {
            version_cmd => q{2>&1 >/dev/null},
            version_regexp =>
q?'my ($version) = /Version:\s+(.*)/xms; if($version) {chomp $version;print $version;last;}'?,
        },
        bedtools => {
            version_cmd => q{--version},
            version_regexp =>
q?'my ($version) = /bedtools\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        bgzip => {
            version_cmd => q{-h 2>&1 >/dev/null},
            version_regexp =>
q?'my ($version) = /Version:\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        bwa => {
            version_cmd => q{2>&1 >/dev/null},
            version_regexp =>
q?'my ($version) = /Version:\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        q{bwa-mem2} => {
            version_cmd => q{ version 2>&1 >/dev/null},
            version_regexp =>
q?'my ($version) = /Version:\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        chanjo => {
            version_cmd => q{--version},
            version_regexp =>
q?'my ($version) = /version\s(\S+)/xms; if($version) {print $version;last;}'?,
        },
        q{configManta.py} => {
            version_cmd    => q{--version},
            version_regexp => q?'chomp;print $_;last;'?,
        },
        fastqc => {
            version_cmd => q{--version},
            version_regexp =>
q?'my ($version) = /FastQC\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        ExpansionHunter => {
            version_cmd => q{--version 2>&1 >/dev/null},
            version_regexp =>
q?'my ($version) = /Hunter\s+(v\d+.\d+.\d+)/xms; if($version) {print $version;last;}'?,
        },
        gatk => {
            version_cmd => q{--java-options "-Xmx1G" --version},
            version_regexp =>
q?'my ($version) = /\(GATK\)\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        q{geneBody_coverage2.py} => {
            version_cmd => q{--version},
            version_regexp =>
q?'my ($version) = /geneBody_coverage2.py\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        genmod => {
            version_cmd => q{--version},
            version_regexp =>
q?'my ($version) = /genmod\s+version:\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        gffcompare => {
            version_cmd => q{--version 2>&1 >/dev/null},
            version_regexp =>
q?'my ($version) = /gffcompare\sv(\S+)/xms; if($version) {print $version;last;}'?,
        },
        q{grep} => {
            version_cmd => q{--version},
            version_regexp =>
q?'my ($version) = /\(GNU\s+grep\)\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        gzip => {
            version_cmd => q{--version},
            version_regexp =>
              q?'my ($version) = /gzip\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        q{infer_exeperiment.py} => {
            version_cmd => q{--version},
            version_regexp =>
q?'my ($version) = /infer_exeperiment.py\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        q{inner_distance.py} => {
            version_cmd => q{--version},
            version_regexp =>
q?'my ($version) = /inner_distance.py\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        q{junction_annotation.py} => {
            version_cmd => q{--version},
            version_regexp =>
q?'my ($version) = /junction_annotation.py\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        mip => {
            version_cmd => q{version},
            version_regexp =>
q?'my ($version) = /mip\s+version\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        multiqc => {
            version_cmd => q{--version},
            version_regexp =>
q?'my ($version) = /version\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        peddy => {
            version_cmd => q{--version},
            version_regexp =>
q?'my ($version) = /version\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        preseq => {
            version_cmd => q{2>&1 >/dev/null},
            version_regexp =>
q?'my ($version) = /Version:\s+(\S+)/xms; if ($version) {print $version; last;}'?,
        },
        picard => {
            version_cmd => q{BamIndexStats 2>&1 >/dev/null},
            version_regexp =>
q?'my ($version) = /Version:\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        pigz => {
            version_cmd => q{--version 2>&1 >/dev/null},
            version_regexp =>
              q?'my ($version) = /pigz\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        plink2 => {
            version_cmd => q{--version},
            version_regexp =>
q?'my ($version) = /PLINK\s+(.*)/xms; if($version) {chomp $version;print $version;last;}'?,
        },
        q{read_distribution.py} => {
            version_cmd => q{--version},
            version_regexp =>
q?'my ($version) = /read_distribution.py\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        q{read_duplication.py} => {
            version_cmd => q{--version},
            version_regexp =>
q?'my ($version) = /read_duplication.py\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        rhocall => {
            version_cmd => q{--version},
            version_regexp =>
q?'my ($version) = /version\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        salmon => {
            version_cmd => q{--version},
            version_regexp =>
q?'my ($version) = /salmon\s+(\S+)/xms; if ($version) {print $version; last;}'?,
        },
        sambamba => {
            version_cmd => q{--version 2>&1 >/dev/null},
            version_regexp =>
q?'my ($version) = /sambamba\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        samtools => {
            version_cmd => q{2>&1 >/dev/null},
            version_regexp =>
q?'my ($version) = /Version:\s+(.*)/xms; if($version) {chomp $version;print $version;last;}'?,
        },
        sed => {
            version_cmd => q{--version},
            version_regexp =>
q?'my ($version) = /\(GNU\s+sed\)\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        STAR => {
            version_cmd => q{--version},
            version_regexp =>
              q?'my ($version) = /(\S+)/xms; if($version) {print $version; last;}'?,
        },
        q{STAR-Fusion} => {
            version_cmd => q{--version},
            version_regexp =>
q?'my ($version) = /version:\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        stranger => {
            version_cmd => q{--version},
            version_regexp =>
              q?'my ($version) = /(\S+)/xms; if($version) {print $version;last;}'?,
        },
        stringtie => {
            version_cmd => q{--version},
            version_regexp =>
              q?'my ($version) = /(\S+)/xms; if($version) {print $version;last;}'?
        },
        svdb => {
            version_cmd => q{--version},
            version_regexp =>
              q?'my ($version) = /SVDB-(\S+)/xms; if($version) {print $version;last;}'?,
        },
        tabix => {
            version_cmd => q{--version},
            version_regexp =>
q?'my ($version) = /\(htslib\)\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        telomerecat => {
            version_regexp =>
q?'BEGIN {my $match = 0}; if ($match) {my ($version) = /(\d+.\d+.\d+)/; print $version; last;} $match=1 if (/\A Version: /xms);'?,
        },
        q{TIDDIT.py} => {
            version_cmd => q{--version},
            version_regexp =>
              q?'my ($version) = /TIDDIT-(\S+)/xms; if($version) {print $version;last;}'?,
        },
        trim_galore => {
            version_cmd => q{--version},
            version_regexp =>
q?'my ($version) = /version\s(\S+)/xms; if($version) {print $version;last;}'?,
        },
        upd => {
            version_cmd => q{--version},
            version_regexp =>
              q?'my ($version) = /(\S+)/xms; if($version) {print $version;last;}'?,
        },
        vcfanno => {
            version_cmd => q{2>&1 >/dev/null},
            version_regexp =>
q?'my ($version) = /version\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        vcf2cytosure => {
            version_cmd => q{-V 2>&1 >/dev/null},
            version_regexp =>
q?'my ($version) = /to\s+cytosure\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        varg => {
            version_cmd => q{--version},
            version_regexp =>
q?'my ($version) = /version\s(\S+)/xms; if($version) {print $version;last;}'?
        },
        vep => {
            version_regexp =>
q?'my ($version) = /ensembl-vep\s+:\s(\d+)/xms; if($version) {print $version;last;}'?,
        },
        vt => {
            version_cmd => q{normalize 2>&1 >/dev/null},
            version_regexp =>
q?'my ($version) = /normalize\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        q{wigToBigWig} => {
            version_cmd => q{2>&1 >/dev/null},
            version_regexp =>
q?'my ($version) = /wigToBigWig\sv\s(\S+)/xms; if($version) {print $version;last;}'?,
        },
        just_to_enable_testing => {
            version_cmd    => q{bwa 2>&1 >/dev/null},
            version_regexp => q?'print $version;last;'?,
        },
    );

    if ( defined $executable_name and exists $executable{$executable_name} ) {

        return %{ $executable{$executable_name} };
    }
    return %executable;
}

sub set_executable_container_cmd {

## Function : Set executable command depending on container manager
## Returns  :
## Arguments: $container_href    => Containers hash {REF}
##          : $container_manager => Container manager
##          : $bind_paths_ref    => Array with paths to bind {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $container_href;
    my $container_manager;
    my $bind_paths_ref;

    my $tmpl = {
        container_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$container_href,
            strict_type => 1,
        },
        container_manager => {
            allow       => [qw{ docker singularity }],
            required    => 1,
            store       => \$container_manager,
            strict_type => 1,
        },
        bind_paths_ref => {
            default     => [],
            store       => \$bind_paths_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Program::Singularity qw{ singularity_exec };
    use MIP::Program::Docker qw{ docker_run };

    my %container_api = (
        docker => {
            arg_href => {
                bind_paths_ref => [],
                image          => undef,
            },
            method => \&docker_run,
        },
        singularity => {
            arg_href => {
                bind_paths_ref        => [],
                singularity_container => undef,
            },
            method => \&singularity_exec,
        },
    );

    my %container_cmd;
  CONTAINER_NAME:
    foreach my $container_name ( keys %{$container_href} ) {

      EXECUTABLE:
        while ( my ( $executable_name, $executable_path ) =
            each %{ $container_href->{$container_name}{executable} } )
        {

            my $container_arg =
              $container_manager eq q{singularity} ? q{singularity_container} : q{image};

            ## Set container option depending on singularity or docker
            $container_api{$container_manager}{arg_href}{$container_arg} =
              $container_href->{$container_name}{uri};

            my @cmds = $container_api{$container_manager}{method}
              ->( { %{ $container_api{$container_manager}{arg_href} } } );


            if ( $executable_path and $executable_path ne q{no_executable_in_image}) {

                push @cmds, $executable_path;
            }
            else {

                push @cmds, $executable_name;
            }
            $container_cmd{$executable_name} = join $SPACE, @cmds;
        }
    }
    return %container_cmd;
}

1;
