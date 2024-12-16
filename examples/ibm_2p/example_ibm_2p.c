/** 2D normal impact of a drop on a cylinder 
 */

#include "axi.h"
#include "navier-stokes/centered.h"
#include "ibm/ibminclude.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "view.h"

#define LARGE 1e36
#define D 0.5

double max_level = 10 ;//max_level for droplet
double max_level_ibm = 0 ;//max_level for IBM
double L = 4. ;//domain size
double t_out = 0.1 ;       
double t_end = 10.01 ;    

/** dimensionless properties, normalized by scaling variables rhol, D, sigma
 */
double rhog=0.0012 ; 
double mul =8.45e-3;
double mug =1.52e-4;
double u0  =3.207;          //falling velocity
double xd0 =0.; 
double yd0 =0.; 

//adapt vars
double femax = 0.01;
double uemax = 0.01;

//accelerations
double gravity = 9.8;     //gravity acceleration (m/s2)
double gx;  //dimensionless gravity in x direction  

int Ninit  = 4096; 

//IBM Vars
double r = 10.;
coord ci = {0.,0.};     // initial coordinates of cylinder

//inlet 
u.n[left]  = dirichlet(u0);
u.t[left]  = dirichlet(0);
p[left]    = neumann(0);
pf[left]   = neumann(0);

//outlet
u.n[right]  = neumann(0);
u.t[right]  = neumann(0);
p[right]    = dirichlet(0);
pf[right]   = dirichlet(0);
//#pf[top]     = neumann(0);
//u.n[top] = neumann (0);

double Rd = 1.;

int main(int argc, char * argv[])
{
  double pullx = 0.;
  if (argc > 1)
    max_level = atoi (argv[1]);
  if (argc > 2)
    L         = atof (argv[2]);
  if (argc > 3)
    u0        = atof (argv[3]);
  if (argc > 4)
    t_out     = atof (argv[4]);
  if (argc > 5)
    t_end     = atof (argv[5]);
  if (argc > 6)
    xd0       = atof (argv[6]);
  if (argc > 7)
    rhog      = atof (argv[7]);
  if (argc > 8)
    mul       = atof (argv[8]);
  if (argc > 9)
    mug       = atof (argv[9]);
  if (argc > 10)
    femax     = atof (argv[10]);
  if (argc > 11)
    uemax     = atof (argv[11]);
  if (argc > 12)
    yd0       = atof (argv[12]);
  if (argc > 13)
    pullx       = atof (argv[13]);
  if (argc > 14)
    r       = atof (argv[14]);
  if (argc > 15)
    max_level_ibm = atoi (argv[15]);
  if (argc > 16)
      Rd          = atof(argv[16]);
  if ( max_level_ibm == 0 ) max_level_ibm = max_level; 
  maxlevel = max_level_ibm;//give var to ibm
  ci.x = pullx;
  size (L);
  origin ( -(L-r*8.), 0.);
  init_grid (Ninit);
  /**
  The liquid phase is water, rho_l=1000 kg/m3, mu_l=1e-3 Ps s; 
  the gas phase is air, rho_g = 1.2 kg/m3, mu_g = 1.7e-5 Pa s; 
  surface tension is sigma=0.0688 N/m; 
  The dimensionless parameters can then be computed based on 
  rho_l, sigma, and D. */
  rho1 = 1., rho2 = rhog; 
  mu1 = mul, mu2 = mug;
  f.sigma = 1.;
  TOLERANCE = 1.e-4; 
  //printf("starting max_level_drop = %f d = %f: max_level_ibm = %f D = %f\n",max_level,Rd*2,maxlevel,r*2);
  run();
}

/**
The initial drop is spherical. */

event init (t = 0)
{
  if (!restore (file = "dump-1.360-")) {
    //refine near the IBM surface and droplet
    //refine((sq(x-ci.x) + sq(y-ci.y) < sq(r+r/10) && sq(x-ci.x) + sq(y-ci.y) > sq(r-r/10) && level < maxlevel) || (sq(x-xd0) + sq(y-yd0) < sq(Rd+Rd/10) && sq(x-xd0) + sq(y-yd0) > sq(Rd-Rd/10) && level < max_level));
    ////solid (vof, isf, - sq(x - ci.x) - sq(y - ci.y) + sq(r));
    //fraction (vof, -sq(y-yd0)-sq(x-xd0)+sq(Rd));
    //fraction (f, -sq(y-yd0)-sq(x-xd0)+sq(Rd));
    astats ss;
    int ic = 0;
    do{
        ic++;
        fraction (vof, -sq(y-ci.y)-sq(x-ci.x)+sq(r));
        fraction (f, -sq(y-yd0)-sq(x-xd0)+sq(Rd));
        ss = adapt_wavelet({vof,f}, (double[]){1.e-30,1.e-30}, maxlevel=maxlevel,minlevel=3);
    }while((ss.nf || ss.nc) && ic < 100);
    fraction (vof, -sq(y-ci.y)-sq(x-ci.x)+sq(r));
    fraction (f, -sq(y-yd0)-sq(x-xd0)+sq(Rd));
    double md = 10000000.;
    foreach( reduction(min:md) ){ 
        md = min(md,Delta);
        //OLDu.x[] = 1.e-7*u0*(1-f[])*(1-vof[]);
        u.x[] = u0*(1-vof[]);
    }
    printf("min delta = %f\n",md);
  }
}

event snapshot (t += t_out; t<=t_end ) {
  char name[80];
  sprintf (name, "dump-%05.3f", t);
  p.nodump = false;
  dump (file = name); // so that we can restart
}

/**
Adapt mesh based on the volume fraction. */
event adapt (i++) {
  adapt_wavelet ({vof,f,u}, (double[]){femax,femax,uemax,uemax}, minlevel = 3, maxlevel = max_level);
}

//log results
double dropletx = 0.;
event logdrop(i++){
    dropletx = 0.;
    double voldrop = 0.,velodrop = 0.,maxdx = 0.,mindx = 100000000.;
    foreach(reduction(+:velodrop) reduction(+:dropletx) reduction(+:voldrop) reduction(max:maxdx) reduction(min:mindx)){
        //compute droplet center and velocity.
	    double vd = 2/pi*dv()*f[];
	    voldrop += vd;
	    velodrop += vd * u.x[];
	    dropletx += vd * x;
        maxdx = max(maxdx,Delta);
        mindx = min(mindx,Delta);
    }
    dropletx /= voldrop;
    velodrop /= voldrop;
    
    //print
#if _MPI
    if(pid() == 0){
#endif
    FILE * fp = fopen("dropletData","a");
    if(t == 0)fprintf(fp,"t dt i xc vc ncells mindx maxdx\n");
    fprintf(fp,"%f %f %d %f %f %ld %f %f\n",t,dt,i,dropletx,velodrop,grid->tn,mindx,maxdx);
    fflush(fp);
    fclose(fp);
#if _MPI
    }
#endif
}

int image_width=1200;
int image_height=1200;
double dtgap = 0.01;
double fov1 = 25;
double fov2 = 1;
double fov3 = 5;
scalar vm[];
char name[80];
int plotstep = 0;
double ifmax = 0.;
double ifmin = 1.;
double vofmax = 0.;
double vofmin = 1.;
double pmax = 0.;
double pmin = 10000000000.;
double uxmax = -10000000000.;
double uxmin = 10000000000.;
double uymax = -10000000000.;
double uymin = 10000000000.;
int shiftscales = 0;
int shiftstop = 5;
event plot(t = 0.; t += t_out){
    double tx1=0.5;
    double ty1=-0.5;
    if(shiftscales == shiftstop){
        //hard reset max/min values every so often to allow color bar to find a happy scale for the current time without changing too drastically
        shiftscales = 0;
        ifmax = 0.;
        ifmin = 1.;
        vofmax = 0.;
        vofmin = 1.;
        pmax = 0.;
        pmin = 10000000000.;
        uxmax = -10000000000.;
        uxmin = 10000000000.;
        uymax = -10000000000.;
        uymin = 10000000000.;
    }
    else{
        shiftscales++;
    }
    foreach(reduction(max:ifmax) reduction(min:ifmin) reduction(max:vofmax) reduction(min:vofmin) reduction(max:pmax) reduction(min:pmin) reduction(max:uxmax) reduction(min:uxmin) reduction(max:uymax) reduction(min:uymin)){
      vm[] = sqrt(u.x[]*u.x[]+u.y[]*u.y[]);
      ifmax   =      max(ifmax,f[]);
      ifmin   =      min(ifmin,f[]);
      vofmax =      max(vofmax,vof[]);
      vofmin =      min(vofmin,vof[]);
      pmax   =      max(pmax,p[]);
      pmin   =      min(pmin,p[]);
      uxmax  =      max(uxmax,u.x[]);
      uxmin  =      min(uxmin,u.x[]);
      uymax  =      max(uymax,u.y[]);
      uymin  =      min(uymin,u.y[]);
    }
    boundary ({vm});

    view (fov = fov1, tx=tx1, ty=ty1, quat = {0,0,0,1}, width = image_width, height = image_height);

    clear();
    box();
    squares("f",max=ifmax,min=ifmin,cbar=true);
    draw_vof ("f");
    sprintf (name, "vof_%05d.png",plotstep);
    save(file=name);
    
    clear();
    box();
    squares("vof",max=vofmax,min=vofmin,map = blue_white_red,cbar=true);
    draw_vof ("vof");
    sprintf (name, "ibm_%05d.png",plotstep);
    save(file=name);
    
    clear();
    box();
    draw_vof ("f");
    draw_vof ("vof");
    squares("u.x", linear = true, max = uxmax, min = uxmin,map = blue_white_red,cbar=true);
    sprintf (name, "full-velx_%05d.png",plotstep);
    save(file=name);
    
    clear();
    box();
    cells();
    draw_vof ("f");
    draw_vof ("vof");
    sprintf (name, "cells_%05d.png",plotstep);
    save(file=name);
    
    clear();
    box();
    draw_vof ("f");
    draw_vof ("vof");
    squares("u.y", linear = true, max = uymax, min = uymin,map = blue_white_red,cbar=true);
    sprintf (name, "full-vely_%05d.png",plotstep);
    save(file=name);
    
    clear();
    box();
    draw_vof ("f");
    draw_vof ("vof");
    squares("p", linear = true, max = pmax, min = pmin,map= blue_white_red,cbar=true);
    sprintf (name, "full-p_%05d.png",plotstep);
    save(file=name);
    
    clear();
    box();
    draw_vof ("f");
    draw_vof ("vof");
    squares("vm", linear = true,map= blue_white_red,cbar=true);
    sprintf (name, "full-vm_%05d.png",plotstep);
    save(file=name);
    
    //zoom on droplet 
    double shiftx = -0.050;
    double lx = 1.0 + shiftx;
    double rx = 0.0 + shiftx;
    double tx2 = lx + (dropletx+L-r*8)/(L)*(rx-lx);//left bound = 0.9, right bound = -0.1 
    double ty2 = -0.01;
    view (fov = fov2, tx=tx2, ty=ty2, quat = {0,0,0,1}, width = image_width, height = image_height);
    clear();
    box();
    draw_vof ("f", lw = 3, lc = {0,0,0});
    draw_vof ("vof");//, lw = 5, lc = {0,0,0});
    squares("p", linear = true, max = pmax, min = pmin,map= blue_white_red,cbar=true);
    sprintf (name, "closedrop-pressure_%05d.png",plotstep);
    save(file=name);
    
    clear();
    box();
    draw_vof ("f", lw = 3, lc = {0,0,0});
    draw_vof ("vof");//, lw = 5, lc = {0,0,0});
    squares("u.x", linear = true, max = uxmax, min = uxmin,map = blue_white_red,cbar=true);
    sprintf (name, "closedrop-velx_%05d.png",plotstep);
    save(file=name);
    
    //zoom on ibm surface
    shiftx = -0.01;
    lx = 1.0 + shiftx;
    rx = 0.0 + shiftx;
    double tx3 = lx + (ci.x+L-r*8)/(L)*(rx-lx);//left bound = 1 ,right bound = 0.0
    double ty3 = -0.1;
    view (fov = fov3, tx=tx3, ty=ty3, quat = {0,0,0,1}, width = image_width, height = image_height);
    clear();
    box();
    draw_vof ("f");
    draw_vof ("vof", lw = 3, lc = {0,0,0});
    squares("p", linear = true, max = pmax, min = pmin,map= blue_white_red,cbar=true);
    sprintf (name, "closeibm-pressure_%05d.png",plotstep);
    save(file=name);
    
    clear();
    box();
    draw_vof ("f");
    draw_vof ("vof", lw = 3, lc = {0,0,0});
    squares("u.x", linear = true, max = uxmax, min = uxmin,map = blue_white_red,cbar=true);
    sprintf (name, "closeibm-velx_%05d.png",plotstep);
    save(file=name);
    
    plotstep++;
}

