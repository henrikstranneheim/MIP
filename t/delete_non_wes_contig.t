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
use MIP::Test::Fixtures qw{ test_log };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Contigs}        => [qw{ delete_non_wes_contig }],
        q{MIP::Test::Fixtures} => [qw{ test_log }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Contigs qw{ delete_non_wes_contig };

diag(   q{Test delete_non_wes_contig from Contigs.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( {} );

## Given contigs, when no prefix
my @contigs = qw{ 1 2 3 4 MT};

## When consensus analysis type is wes
my @no_wes_contigs = delete_non_wes_contig(
    {
        consensus_analysis_type => q{wes},
        contigs_ref             => \@contigs,
        contig_names_ref        => [qw{ M MT }],
    }
);

## Define expected outcome
my @expected_no_wes_contigs = qw{ 1 2 3 4 };

## Then remove the non wes contigs
is_deeply( \@no_wes_contigs, \@expected_no_wes_contigs, q{Removed non wes contigs} );

## Given contigs, when a non consensus type of run i.e. mixed
@no_wes_contigs = delete_non_wes_contig(
    {
        consensus_analysis_type => q{mixed},
        contigs_ref             => \@contigs,
    }
);

## Then remove the non wes contigs
is_deeply( \@no_wes_contigs, \@expected_no_wes_contigs,
    q{Remove non wes contigs for mixed analysis run} );

## Given contigs, when a wgs consensus type of run
my @has_wes_contigs = delete_non_wes_contig(
    {
        consensus_analysis_type => q{wgs},
        contigs_ref             => \@contigs,
    }
);

## Define expected outcome
my @expected_wes_contigs = qw{ 1 2 3 4 MT };

## Then keep the non wes contigs
is_deeply( \@has_wes_contigs, \@expected_wes_contigs, q{Keept non wes contigs for wgs} );

## Given contigs, when a wts consensus type of run
@has_wes_contigs = delete_non_wes_contig(
    {
        consensus_analysis_type => q{wts},
        contigs_ref             => \@contigs,
    }
);

## Then keep the non wes contigs
is_deeply( \@has_wes_contigs, \@expected_wes_contigs, q{Keept non wes contigs for wts} );
done_testing();
