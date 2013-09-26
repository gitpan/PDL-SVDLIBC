
#
# GENERATED WITH PDL::PP! Don't modify!
#
package PDL::SVDLIBC;

@EXPORT_OK  = qw( PDL::PP _svdccsencode  svdlas2a PDL::PP svdlas2  svdlas2ad PDL::PP svdlas2d );
%EXPORT_TAGS = (Func=>[@EXPORT_OK]);

use PDL::Core;
use PDL::Exporter;
use DynaLoader;



   $PDL::SVDLIBC::VERSION = 0.08;
   @ISA    = ( 'PDL::Exporter','DynaLoader' );
   push @PDL::Core::PP, __PACKAGE__;
   bootstrap PDL::SVDLIBC $VERSION;




=pod

=head1 NAME

PDL::SVDLIBC - PDL interface to Doug Rohde's SVD C Library

=head1 SYNOPSIS

 use PDL;
 use PDL::SVDLIBC;

 ##---------------------------------------------------------------------
 ## Input matrix (dense)
 ##---------------------------------------------------------------------
 $n = 100;                  ##-- number of columns
 $m = 50;                   ##-- number of rows
 $a = random(double,$n,$m); ##-- random matrix

 ##---------------------------------------------------------------------
 ## Output pdls
 ##---------------------------------------------------------------------
 $d  = $n;                   ##-- max number of output dimensions
 $ut = zeroes(double,$m,$d); ##-- left singular components
 $s  = zeroes(double,$d);    ##-- singular values (diagnonal vector)
 $vt = zeroes(double,$n,$n); ##-- right singular components

 ##---------------------------------------------------------------------
 ## Singular Value Decomposition (dense)
 ##---------------------------------------------------------------------
 svdlas2d($a, $maxiters, $end, $kappa, $ut, $s, $vt);

 ##---------------------------------------------------------------------
 ## Singular Value Decomposition (sparse)
 ##---------------------------------------------------------------------
 use PDL::CCS;
 ($ptr,$rowids,$nzvals) = ccsencode($a);
 $ptr->reshape($ptr->nelem+1);
 $ptr->set(-1, $rowids->nelem);

 svdlas2($ptr, $rowids, $nzvals, $m, $maxiters, $end, $kappa, $ut, $s, $vt);


=head1 DESCRIPTION

PDL::SVDLIBC provides a PDL interface to the SVDLIBC routines
for singular value decomposition of large sparse matrices.
SVDLIBC is available from http://tedlab.mit.edu/~dr/SVDLIBC/

=cut







=head1 FUNCTIONS



=cut





=pod

=head1 SVDLIBC Globals

There are several global data structures still lurking in the
SVDLIBC code, so expect problems if you are trying to run more
than one 'las2' procedure at once (even in different processes).

PDL::SVDLIBC provides access to (some of) the SVDLIBC globals
through the following functions, which are not exported.

=cut



=pod

=head2 PDL::SVDLIBC::verbosity()

=head2 PDL::SVDLIBC::verbosity($level)

Get/set the current SVDLIBC verbosity level.
Valid values for $level are between 0 (no messages) and
2 (many messages).

=cut




=pod

=head2 PDL::SVDLIBC::svdVersion()

Returns a string representing the SVDLIBC version
this module was compiled with.

=cut




=pod

=head1 SVD Utilities

=cut





=head2 _svdccsencode

=for sig

  Signature: (double a(n,m); int      [o]ptr(n1); int      [o]rowids(nnz); double [o]nzvals(nnz))


=for ref

info not available


=for bad

_svdccsencode does not process bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*_svdccsencode = \&PDL::_svdccsencode;




=pod

=head2 svdlas2a

=for sig

    indx    ptr(nplus1);
    indx    rowids(nnz);
    double  nzvals(nnz);
    indx    nrows();      ##-- default: max($rowids)+1
    int     d();          ##-- default: nplus1-1
    int     iterations(); ##-- default: 2*$d
    double  end(2);       ##-- default: [-1e-30,1e-30]
    double  kappa();      ##-- default: 1e-6
    double  [o]ut(m,d);   ##-- default: new
    double  [o] s(d);     ##-- default: new
    double  [o]vt(n,d);   ##-- default: new

Uses a variant of the single-vector Lanczos method (Lanczos, 1950)
to compute the singular value decomposition of a sparse matrix with
$nrows() rows and data encoded
in Harwell-Boeing sparse format in the input parameters $ptr(), $rowids(),
and $nzvals().  See L<"PDL::CCS"> for a way to acquire these parameters
from a dense input matrix, but note that for svdlas2(), the
column pointer $ptr() is of size ($n+1) for a dense matrix $a with
$n columns, where $ptr($n)==$nnz is the total number of nonzero
values in $a.

$iterations() is the maximum number of Lanczos iterations to perform.

$end() specifies two endpoints of an interval within which all unwanted
eigenvalues lie.

$kappa() is a double containing the relative accuracy of Ritz
values acceptable as eigenvalues.

The left singular components are returned in the matrix $ut(),
the singular values themselved in the vector $s(), and the
right singular components in the matrix $vt().  Note that
$ut() and $vt() are transposed, and must be specified explicitly
in the call, so that the degree of reduction (the size parameter $d)
can be determined.  If $d==$n, then a full decomposition
will be computed, and on return, $ut() and $vt() should be transposed
instances of the matrices $u() and $v() as returned by PDL::MatrixOps::svd().

The Lanczos method as used here seems to be consistently the
fastest. This algorithm has the drawback that the low order singular
values may be relatively imprecise, but that is not a problem for most
users who only want the higher-order values or who can tolerate some
imprecision.

See also: svdlas2d()

=cut

## ($iters,$end,$kappa,$ut,$s,$vt) = svddefaults($nrows,$ncols,$d, $iters,...)
## + returns default values
sub svddefaults {
    my ($nrows,$ncols,$d, $iters,$end,$kappa,$ut,$s,$vt) = @_;
    $nrows = $nrows->at(0) if (UNIVERSAL::isa($nrows,'PDL'));
    $ncols = $ncols->at(0) if (UNIVERSAL::isa($ncols,'PDL'));
    $d     = $ncols if (!defined($d));
    $iters = 2*$d if (!defined($iters));
    $end   = pdl(double,[-1e-30,1e-30]) if (!defined($end));
    $kappa = pdl(double,1e-6) if (!defined($kappa));
    $ut    = PDL->zeroes(double,$nrows,$d) if (!defined($ut));
    $s     = PDL->zeroes(double,$d) if (!defined($s));
    $vt    = PDL->zeroes(double,$ncols,$d) if (!defined($vt));
    return ($iters,$end,$kappa,$ut,$s,$vt);
}

sub svdlas2a {
    my ($ptr,$rowids,$nzvals, $nrows,$d, @args) = @_;
    $nrows = $rowids->flat->max+1 if (!defined($nrows));
    @args = svddefaults($ptr->dim(0)-1,$nrows,$d,@args);
    svdlas2($ptr,$ropwids,$nzvals,$nrows,@args);
    return @args[3..5];
}





=head2 svdlas2

=for sig

  Signature: (
    int   ptr(nplus1);
    int   rowids(nnz);
    double  nzvals(nnz);
    int   nrows();
    int     iterations();
    double  end(2);
    double  kappa();
    double  [o]ut(m,d);
    double  [o] s(d);
    double  [o]vt(n,d);
    )


Guts for svdlas2a().
No default instantiation, and slightly different calling conventions.


=for bad

svdlas2 does not process bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*svdlas2 = \&PDL::svdlas2;




=pod

=head2 svdlas2ad

=for sig

    double  a(n,m);
    int     d();          ##-- default: $n
    int     iterations(); ##-- default: 2*$d
    double  end(2);       ##-- default: [-1e-30,1e-30]
    double  kappa();      ##-- default: 1e-6
    double  [o]ut(m,d);   ##-- default: new
    double  [o] s(d);     ##-- default: new
    double  [o]vt(n,d);   ##-- default: new

As for svdlas2(), but implicitly converts the dense input matrix
$a() to sparse format before computing the decomposition.

=cut

sub svdlas2ad {
    my ($a,$d, @args) = @_;
    @args = svddefaults($a->dim(1),$a->dim(0),$d,@args);
    svdlas2d($a,@args);
    return @args[3..5];
}





=head2 svdlas2d

=for sig

  Signature: (
    double  a(n,m);
    int     iterations();
    double  end(2);
    double  kappa();
    double  [o]ut(m,d);
    double  [o] s(d);
    double  [o]vt(n,d);
    )


Guts for _svdlas2d().


=for bad

svdlas2d does not process bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*svdlas2d = \&PDL::svdlas2d;




##---------------------------------------------------------------------
=pod

=head1 ACKNOWLEDGEMENTS

Perl by Larry Wall.

PDL by Karl Glazebrook, Tuomas J. Lukka, Christian Soeller, and others.

SVDLIBC by Dough Rohde.

SVDPACKC by Michael Berry, Theresa Do, Gavin O'Brien, Vijay Krishna and Sowmini Varadhan.

=cut

##----------------------------------------------------------------------
=pod

=head1 KNOWN BUGS

Globals still lurk in the depths of SVDLIBC.

=cut


##---------------------------------------------------------------------
=pod

=head1 AUTHOR

Bryan Jurish E<lt>moocow@cpan.orgE<gt>

=head2 Copyright Policy

Copyright (C) 2005-2013, Bryan Jurish. All rights reserved.

This package is free software, and entirely without warranty.
You may redistribute it and/or modify it under the same terms
as Perl itself.

=head1 SEE ALSO

perl(1), PDL(3perl), PDL::CCS(3perl), SVDLIBC documentation.

=cut



;



# Exit with OK status

1;

		   