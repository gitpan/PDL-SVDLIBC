# -*- Mode: CPerl -*-
# t/01_svd.t: test las2 svd

$TEST_DIR = './t';
#use lib qw(../blib/lib ../blib/arch); $TEST_DIR = '.'; # for debugging

# load common subs
use Test;
do "$TEST_DIR/common.plt";
use PDL;
use PDL::SVDLIBC;

BEGIN { plan tests=>6, todo=>[]; }

##-- setup
$a = pdl(double,
	 [[10,0,0,0,-2,0,0],
	  [3,9,0,0,0,3,1],
	  [0,7,8,7,0,0,0],
	  [3,0,8,7,5,0,1],
	  [0,8,0,9,9,13,0],
	  [0,4,0,0,2,-1,1]]);

$ptr=pdl(long,[0,3,7,9,12,16,19, 22]);
$rowids=pdl(long,[0,1,3,1,2,4,5,2,3,2,3,4,0,3,4,5,1,4,5,1,3,5]);
$nzvals=pdl(double,[10,3,3,9,7,8,4,8,8,7,7,9,-2,5,9,2,3,13,-1,1,1,1]);

($n,$m) = $a->dims;



##-- common pars
$iters = pdl(long,14);
$end   = pdl(double,[-1e-30,1e-30]);
$kappa = pdl(double,1e-6);

##-- $d==$n: expect
$d = $n;

$ut_want =
  pdl(double,
      [[0.000536870080275633,0.272178136497484,0.402223741637742,0.342152709202008,0.79777433508708,0.103065364251536],
       [0.28527967615187,-0.101703564414206,0.539263413191387,0.608398065357372,-0.497242386011596,-0.00828650566923618],
       [0.769356159612389,0.57278591222696,-0.0758981711367609,-0.248163283906345,-0.0632588945976401,0.0930599957769441],
       [0.472065824089417,-0.458560287003903,-0.509978454684947,0.383227968381176,0.285120934392026,-0.280429445117102],
       [-0.14350173411983,0.181133269960328,-0.451985822294512,0.418727072563739,-0.10964420848225,0.744951403456183],
       [0.28856096035857,-0.586860245243611,0.277995130671446,-0.359168853225782,0.137799246742958,0.589114109915952],
       [0,0,0,0,0,0]]);

$vt_want =
  pdl(double,
      [[0.0792510967013729,0.517073465756313,0.255329527055914,0.531264640948551,0.389994814017901,0.475265099143429,0.030759374976173],
       [0.337930883175723,-0.0889921877945774,0.709519096957111,0.274992613247936,-0.156128150575161,-0.522479963841427,0.0385163653800307],
       [0.78833914428785,0.408380265790316,-0.235798012061466,-0.258106221868997,-0.287664388100387,0.073030050581229,0.0379899903637539],
       [0.494549455431307,-0.719339483237642,-0.1115712708961,0.184722915705758,0.327586707049301,0.287324708955066,-0.0391446152142804],
       [0.0948079103988546,0.147960483656045,-0.0691937827562251,-0.317169803014294,0.749942274518314,-0.42309612560424,0.349729817537494],
       [0.0411788312667813,0.106645671467875,-0.562706510404694,0.582278834218919,0.0393876958323224,-0.483779794610716,-0.30927247974065],
       [0,0,0,0,0,0,0]]);

$s_want = pdl(double,
	      [23.32284744104,12.9401616781924,10.9945440916999,9.08839598479768,3.84528764361343,1.1540470359863,0]);


##-- test 1..3 : svdlas2
svdlas2($ptr,$rowids,$nzvals, $m,
	$iters, $end, $kappa,
	($ut=zeroes(double,$m,$d)),
	($s=zeroes(double,$d)),
	($vt=zeroes(double,$d,$n)),
       );
isok("svdlas2,d=n:ut", all($ut->approx($ut_want)));
isok("svdlas2,d=n:s", all($s->approx($s_want)));
isok("svdlas2,d=n:vt", all($s->approx($vt_want)));


##-- test 4..6 : svdlas2d
svdlas2d($a,
	 $iters, $end, $kappa,
	 ($ut=zeroes(double,$m,$d)),
	 ($s=zeroes(double,$d)),
	 ($vt=zeroes(double,$d,$n)),
	);
isok("svdlas2d,d=n:ut", all($ut->approx($ut_want)));
isok("svdlas2d,d=n:s",  all($s->approx($s_want)));
isok("svdlas2d,d=n:vt", all($s->approx($vt_want)));


print "\n";
# end of t/01_svd.t
