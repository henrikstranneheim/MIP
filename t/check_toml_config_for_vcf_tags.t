#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw { :all };
use Modern::Perl qw{ 2018 };
use Readonly;
use Test::Trap qw{ :stderr:output(systemsafe) };

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COLON $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log test_standard_cli };
use MIP::Test::Writefile qw{ write_toml_config };

my $VERBOSE = 1;
our $VERSION = 1.00;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Reference}      => [qw{ check_toml_config_for_vcf_tags }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Reference qw{ check_toml_config_for_vcf_tags };

diag(   q{Test check_toml_config_for_vcf_tags from Reference.pm v}
      . $MIP::Reference::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

# Create log object
my $log = test_log( { log_name => q{MIP}, no_screen => 0, } );

my $cluster_reference_path = catdir( dirname($Bin), qw{ t data references } );
my $toml_template_path     = catfile( $cluster_reference_path,
    q{grch37_frequency_vcfanno_filter_config_template-v1.0-.toml} );
my $toml_config_path =
  catfile( $cluster_reference_path, q{grch37_vcfanno_config-v1.0-.toml} );
my $spliceai_annotation_path =
  catfile( $cluster_reference_path, q{grch37_spliceai_scores_raw_snv_-v1.3-.vcf.gz} );

## Update path in toml config
write_toml_config(
    {
        test_reference_path => $cluster_reference_path,
        toml_config_path    => $toml_config_path,
        toml_template_path  => $toml_template_path,
    }
);

## Given a toml config with preops
my %active_parameter_test = (
    binary_path        => { bcftools => q{bcftools}, },
    variant_annotation => 1,
    vcfanno_config     => $toml_config_path,
);

my %preop_annotations = check_toml_config_for_vcf_tags(
    {
        active_parameter_href => \%active_parameter_test,
    }
);

## Then return the preop annotation hash
my %expected = (
    catfile($spliceai_annotation_path) => {
        annotation => [
            {
                file   => $spliceai_annotation_path,
                fields => [q{SpliceAI}],
                ops    => [q{lua:spliceai_max_score(vals)}],
                names  => [q{SpliceAI_DS_max}],
            },
        ],
    },
);
is_deeply( \%preop_annotations, \%expected, q{Set preops} );

## Given a vcfanno annotation request that tries to use a non-existing vcf tag
$toml_template_path =
  catfile( $cluster_reference_path, q{grch37_vcfanno_config_bad_template-v1.0-.toml} );
$toml_config_path =
  catfile( $cluster_reference_path, q{grch37_vcfanno_config-v1.0-.toml} );

write_toml_config(
    {
        test_reference_path => $cluster_reference_path,
        toml_config_path    => $toml_config_path,
        toml_template_path  => $toml_template_path,
    }
);

trap {
    check_toml_config_for_vcf_tags(
        {
            active_parameter_href => \%active_parameter_test,
        }
    )
};

## Then print fatal log message and exit
ok( $trap->exit, q{Exit when vcf annotation tags are missing } );
like( $trap->stderr, qr/FATAL/xms, q{Throw fatal log message} );

unlink $toml_config_path;
done_testing();
