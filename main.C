#include <conf.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <alg/alg_plaq.h>
#include <alg/alg_pbp.h>
#include <alg/alg_eig.h>
#include <alg/alg_lanczos.h>
#include <alg/alg_meas.h>
#include <alg/alg_s_vacpol.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/command_line.h>
#include <util/ReadLatticePar.h>
#include <util/qio_readLattice.h>
#include <util/qcdio.h>
#include <util/site.h>
#include <util/enum_func.h>
#include <alg/s_spect_arg.h>
#include <alg/alg_rnd_gauge.h>

#include <util/eigen_container.h>


USING_NAMESPACE_CPS
using namespace std;

DoArg do_arg;
NoArg no_arg;
PbpArg pbp_arg;
StagQuarkArg sq_arg;
StagQuarkArg sq_arg2;
CommonArg common_arg;
CommonArg hvp_common_arg;
LanczosArg eig_arg;
EigArg ritz_arg;
MatrixPolynomialArg cheby_arg;

void movefloattoFloat(Float* out, float* in, int f_size){

  float flt;
  for(int i=0;i<f_size;i++){
    flt = in[i];
    out[i] = (Float)flt;
  }
};

//------------------------------------------
// should go to eigen_container.C later, someday
//------------------------------------------
// needed to declare globally
std::vector<EigenCache*> cps::EigenCacheList(0);

//Search contents that match to arguments, return 0 if not found
EigenCache* cps::EigenCacheListSearch( char* fname_root_bc, int neig )
{
  EigenCache* ecache=0;

  for(int i=0; i< EigenCacheList.size(); ++i  )
    if( EigenCacheList[i]-> is_cached( fname_root_bc, neig ) )
      ecache = EigenCacheList[i];

  return ecache;
}

static Float link_trace(Lattice& lattice)
{
  Float linktrace(0);
  int is;
  Matrix *m =  lattice.GaugeField(); 
  for(is=0;is< GJP.VolSites()*4; is++){
    linktrace += m->ReTr();
    m++;
  }
  linktrace /= (GJP.VolSites()*12.0);
  printf(" link trace %e\n",linktrace );


  return linktrace ; 
}

int main(int argc,char *argv[])
{

  char *cname = argv[0] ;
  char *fname = "main()" ;
  char* filename;
  filename = (char*) smalloc( 128*sizeof(char) );
  
  Start(&argc, &argv);
  
  //------------------------------
  //set log directory
  //------------------------------
  char log_dir[255];
  sprintf(&log_dir[0],"IOLOG");
  
  if ( argc!=10) { 
    if(!UniqueID())printf("(exe) do_arg lanczos_arg cheby_arg sq_arg start inc tstart tinc work-directory\n");
    exit(-1);
  }
   
  // load defaults
  
  // move to working sub-directory (first arg on command line after exec.)
  // should be relative to /host.../username
  
  chdir(argv[9]);

  if ( !do_arg.Decode(argv[1],"do_arg") ) 
    { 
      ERR.General(fname,fname,"Decoding of do_arg failed\n");
    }
  if ( !eig_arg.Decode(argv[2],"eig_arg") ) 
    { 
      ERR.General(fname,fname,"Decoding of eig_arg failed\n");
    }
  if ( !cheby_arg.Decode(argv[3],"cheby_arg") ) 
    { 
      ERR.General(fname,fname,"Decoding of cheby_arg failed\n");
    }
  if ( !sq_arg.Decode(argv[4],"sq_arg") ) 
    { 
      sq_arg.Encode("sq_arg.dat","sq_arg");
      ERR.General(fname,fname,"Decoding of sq_arg failed\n");
    }
  eig_arg.matpoly_arg = (Float*)&cheby_arg; // ugly thing to workaround vml 

  sq_arg2 = sq_arg;

  int start = atoi(argv[5]);
  int inc = atoi(argv[6]);
  int tstart = atoi(argv[7]);
  int tinc = atoi(argv[8]);

  
  VRB.Level(do_arg.verbose_level);
  GJP.Initialize(do_arg);

  // has to be explicitly Fhisq so we can point to smeared links 
  GwilsonFhisq lattice;
  //GwilsonFp4 lattice;
  //GwilsonFstag lattice;
  qio_readLattice rl(do_arg.start_conf_filename, lattice, argc, argv);
  common_arg.set_filename("plaq.dat");
  AlgPlaq plaq(lattice,&common_arg,&no_arg);
  plaq.run();

  do_arg.Encode("do_arg.dat","do_arg");
  eig_arg.Encode("eig_arg.dat","eig_arg");
  sq_arg.Encode("sq_arg.dat","sq_arg");

  /*
   *  Compute eigenvectors and values
   */

  int traj = GJP.Traj();

  char eig_fname[1024];
  char eig_fname2[1024];
  char fname_bc[1024];
  char fname_bc2[1024];
  
  Float etime = time_elapse();
  char cache_name[1024];

#if 1
  // no twist
  EigenCache *ecache;
  snprintf(cache_name,1024,"cache");
  ecache = new EigenCache(cache_name);
  snprintf(eig_fname,1024, "%s.twist%s%s%s%s", 
	   sq_arg.cg.fname_eigen,
	   "0","0","0","0");
  eig_arg.file = eig_fname;
  sq_arg.cg.fname_eigen = eig_fname;// need this later for CG in HVP
  snprintf(fname_bc,1024, "%s.bc%d%d%d%d", 
	   eig_arg.file,
	   GJP.Bc(0),GJP.Bc(1),GJP.Bc(2),GJP.Bc(3));
  EigenCacheList. push_back( ecache );

  AlgLanczos eig(lattice, &common_arg, &eig_arg, ecache);
  int Ncb = eig.NumChkb(eig_arg.RitzMat_lanczos);
  int fsize = GJP.VolNodeSites() * lattice.FsiteSize() * Ncb / 2 /2 ;
  int neig;
#else

  ritz_arg.N_eig = 20;
  ritz_arg.pattern_kind = LIN;
  ritz_arg.Mass_init = 0.1;
  ritz_arg.Mass_final = 0.2;
  ritz_arg.Mass_step = 0.2;
  ritz_arg.RsdR_a = 1.0E-8;
  ritz_arg.RsdR_r = 1.0E-8;
  ritz_arg.Rsdlam = 1.0E-8;
  ritz_arg.Kalk_Sim = 1;
  ritz_arg.Cv_fact = 0.1;
  ritz_arg.N_min = 10;
  ritz_arg.N_max = 200;
  ritz_arg.N_KS_max = 50;
  ritz_arg.n_renorm = 15;
  ritz_arg.MaxCG = 500;
  ritz_arg.ProjApsiP = 1;
  ritz_arg.RitzMatOper = MATPCDAG_MATPC;
  ritz_arg.print_hsum = 1;
  ritz_arg.hsum_dir = 3;

  AlgEig eig(lattice, &common_arg, &ritz_arg);
  //eig.run();

#endif
#if 1
  if( eig_arg.nk_lanczos_vectors > 0) {

    neig  = eig_arg.nk_lanczos_vectors+eig_arg.np_lanczos_vectors;
    ecache->alloc( fname_bc, neig, fsize );

    eig.run( 0, 0, "" );

    neig = sq_arg.cg.neig;
    ecache->set_neig(neig);

    etime = time_elapse();
    if(!UniqueID())printf("Time for Lanczos %g\n",etime);
    
  } else {
    
    const int n_fields =  GJP.SnodeSites(); 
    const int f_size_per_site = lattice.FsiteSize() / n_fields / 2  ;
    // the two is for checkerboard, only 1/2 sites
    
    neig = sq_arg.cg.neig;
    ecache->alloc( fname_bc, neig, fsize );
    EigenContainer eigcon( lattice, fname_bc, neig, f_size_per_site, n_fields, ecache );
    
    Float* eval;
    if(neig) eval = eigcon.load_eval();
    Float* evecFloat = (Float*)smalloc(2*fsize * sizeof(Float));
    
    for(int iev=0; iev< neig; iev++){

      Vector* evec= eigcon.nev_load( iev );

      movefloattoFloat(evecFloat,(float*)evec,2*fsize);
      
      if(iev % 1 == 0){  // Let's  check
	
	Float res_chk, eval_chk;
	Float mass = eig_arg.mass;  	
	
	eigcon. nev_check( (Vector*)evecFloat, mass, &res_chk, &eval_chk );

	Float ev = *(eval+iev);
	if( fabs(res_chk) > 1e-5 ) 
	  VRB.Warn(cname,fname,"nev_check index %d eigval %g mass %g res %e > 1e-5\n",
			   iev, ev, mass, res_chk);

	if( fabs(eval_chk- ev) > 1e-5 ) 
	  VRB.Warn(cname,fname,"nev_check index %d mass %g eval_chk %e eval %e, abs_err %e > 1e-5\n", 
		   iev, mass, eval_chk, ev, fabs(eval_chk-ev) );
	
      }
    }

#if 0
Vector* evec= eigcon.nev_load( 0 );
movefloattoFloat(evecFloat,(float*)evec,2*fsize);
for(int i=0;i<6*GJP.VolSites()/2;i+=2){
   printf("EVEC 0 at site 0: %d %10g %10g\n",i,evecFloat[i],evecFloat[i+1]);
}
#endif
  }
 
#endif

  // now calc HVP

  char dfname[1024];

#if 0
  // Exact
  {
    //sq_arg. cg. Inverter = LOWMODEAPPROX;
    snprintf(dfname, 1024,"hvpExact.%d.dat", sq_arg.cg.neig);
    hvp_common_arg.set_filename(dfname);
    AlgVacPolStag hvp(lattice, &hvp_common_arg, &sq_arg, &sq_arg2);
    //hvp.VacPolStagConsConsAMA(start,GJP.Sites(0),tstart,GJP.Sites(3)/8);
    // do all the even sites
    //hvp.VacPolStagConsConsAMAnoTwist(start,GJP.Sites(0),tstart,GJP.Sites(3));
    hvp.VacPolStagConsConsAMAnoTwist(0,1,0,4);
  }
#endif
  sq_arg. cg. stop_rsd = sq_arg. cg. ama_stop_rsd;
  sq_arg2. cg. stop_rsd = sq_arg. cg. ama_stop_rsd;
#if 0
  // ama on new sources for sub
  {
    snprintf(dfname, 1024,"hvpSUB.%d.dat", sq_arg.cg.neig);
    hvp_common_arg.set_filename(dfname);
    AlgVacPolStag hvp(lattice, &hvp_common_arg, &sq_arg,&sq_arg2);
    //hvp.VacPolStagConsConsAMA(start,GJP.Sites(0),tstart,GJP.Sites(3)/8);
    hvp.VacPolStagConsConsAMAnoTwist(start,GJP.Sites(0),tstart,GJP.Sites(3)/8);
  } 
#endif
#if 0
  // ama on new sources for sub
  {
    snprintf(dfname, 1024,"hvpAMA.%d.dat", sq_arg.cg.neig);
    hvp_common_arg.set_filename(dfname);
    AlgVacPolStag hvp(lattice, &hvp_common_arg, &sq_arg,&sq_arg2);
    //hvp.VacPolStagConsConsAMA(start,inc,tstart,tinc);
    hvp.VacPolStagConsConsAMAnoTwist(start,inc,tstart,tinc);
  } 
#endif
#if 1
  // lma volume average
  {
    //snprintf(dfname, 1024,"hvpAMA.%d.dat", sq_arg.cg.neig);
    //hvp_common_arg.set_filename(dfname);
    AlgVacPolStag hvp(lattice, &hvp_common_arg, &sq_arg,&sq_arg2);
    hvp.VacPolStagConsConsLMA(5);
  } 
#endif

  sfree(filename);

  EigenCacheListCleanup( );

  End();
  
  return 0;
} 


// Cleanup list EigenCache, it also destroies contents pointed by the elements.
void cps::EigenCacheListCleanup( )
{
  for(size_t i=0;i< EigenCacheList. size(); ++i){
    EigenCacheList[i]-> dealloc();
  }
  EigenCacheList.clear() ;
}
