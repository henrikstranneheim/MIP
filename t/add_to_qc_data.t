#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir };
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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = ( q{MIP::Qc_data} => [qw{ add_to_qc_data }], );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Qc_data qw{ add_to_qc_data };

diag(   q{Test add_to_qc_data from Qc_data.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given recipe metric when key value pair
my $attribute = q{greeting};
my $recipe    = q{japan_greetings};
my %qc_data;
my %qc_header;
my $value          = q{konnichi wa};
my %qc_recipe_data = ( $recipe => { $attribute => [$value], }, );
my %sample_info;

add_to_qc_data(
    {
        attribute           => $attribute,
        qc_data_href        => \%qc_data,
        qc_header_href      => \%qc_header,
        qc_recipe_data_href => \%qc_recipe_data,
        recipe              => $recipe,
        sample_info_href    => \%sample_info,
    }
);

## Then qc data metric should be present on case level
is( $qc_data{recipe}{$recipe}{$attribute}, $value, q{Set qc data metric key value pair} );

## Given recipe metric when key array values
my $attribute_goodbye = q{goodbye};
my $recipe_goodbye    = q{japanese goodbye};
my @goodbyes          = ( q{sayonara}, q{otsukaresama deshita}, );
push @{ $qc_recipe_data{$recipe_goodbye}{$attribute_goodbye} }, @goodbyes;

add_to_qc_data(
    {
        attribute           => $attribute_goodbye,
        qc_data_href        => \%qc_data,
        qc_header_href      => \%qc_header,
        qc_recipe_data_href => \%qc_recipe_data,
        recipe              => $recipe_goodbye,
        sample_info_href    => \%sample_info,
    }
);

## Then qc data metric should be present on case level
is_deeply( \@{ $qc_data{recipe}{$recipe_goodbye}{$attribute_goodbye} },
    \@goodbyes, q{Added qc data metrics key array pair} );

## Clean-up from previous tests
delete $qc_data{metrics};

## Given an infile and sample_id
my $infile    = q{an_infile};
my $sample_id = q{sample_1};

add_to_qc_data(
    {
        attribute           => $attribute,
        infile              => $infile,
        qc_data_href        => \%qc_data,
        qc_header_href      => \%qc_header,
        qc_recipe_data_href => \%qc_recipe_data,
        recipe              => $recipe,
        sample_info_href    => \%sample_info,
        sample_id           => $sample_id,
    }
);
my %expected_metric = (
    metrics => [
        {
            header => undef,
            id     => $sample_id,
            input  => $infile,
            name   => $attribute,
            value  => $value,
            step   => $recipe,
        },
    ],
);

is_deeply( $qc_data{metrics}, $expected_metric{metrics}, q{Set metrics in qc_data} );

done_testing();
