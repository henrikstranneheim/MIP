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


## Constants
Readonly my $COMMA => q{,};
Readonly my $COLON => q{:};
Readonly my $SPACE => q{ };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Language::Shell} => [qw{ clear_trap }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Language::Shell qw{ clear_trap };

diag(   q{Test clear_trap from SHELL.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

# Create anonymous filehandle
my $filehandle = IO::Handle->new();

# For storing info to write
my $file_content;

## Store file content in memory by using referenced variable
open $filehandle, q{>}, \$file_content
  or croak q{Cannot write to} . $SPACE . $file_content . $COLON . $SPACE . $OS_ERROR;

## Given a filehandle
clear_trap( { filehandle => $filehandle } );

# Close the filehandle
close $filehandle;

## Then trap comment and trap should be written to file
my ($clear_trap_command) = $file_content =~ /^(## Clear trap)/ms;

ok( $clear_trap_command, q{Wrote clear trap title} );

my ($trap_command) = $file_content =~ /^(trap\s+['-'])/mxs;

ok( $trap_command, q{Wrote trap command} );

done_testing();
