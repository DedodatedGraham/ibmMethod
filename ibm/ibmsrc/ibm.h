#define BGHOSTS 2
#include "fractions.h"
#include "ibm-utils.h"
extern double max_level;

coord vc;          // object's imposed velocity
scalar ibm[];

//vector velocityGrad[];

vector forceTotal[];//Need to use in AMR statement, typically works well with error associated with velocity
face vector faceForce[];
int Ni = 10;               // # of multi-direct forcing iterations

#if AXI
#define ibmdv() (sq(Delta)*fabs(cm[]))
#else
#define ibmdv() (sq(Delta)*cm[])
#endif

event init (t = 0)
{
    ibm.refine = ibm.prolongation = fraction_refine;
    ibm.dirty = true;//set prolongation 
    foreach(){
        foreach_dimension() { 
            //velocityGrad.x[] = 0.;
            if (ibm[] == 1){
                u.x[] = vc.x;
	        }
            //set to either 0. or starting a value depending on initalization
        }
    }
}

//allows for our reflection boundary condidion to be set on a vector field
void setRefBC(vector *set){
#if AXI//assigns needed boundary conditions for axi 
    //This BC covers a simple case of a IBM method which crosses the domain on the x axis 
    vector vec;
    foreach_boundary(bottom){
    	for(vec in set){
    	    if(vec.x[] != 0. || vec.y[] != 0. || vec.x[0,1] != 0. || vec.y[0,1] != 0){
    	    	vec.x[0,-1] = vec.x[];                                                         
    	    	vec.x[0,-2] = vec.x[0,1];                                                      
    	    	vec.y[0,-1] = -1*vec.y[];                                                      
    	    	vec.y[0,-2] = -1*vec.y[0,1];                                                   
    	    }
    	}
    }
#endif
#if _MPI
    boundary(set);
#endif
}
//
//extern coord ci;
//extern double dropletx,L,r,fov1,fov2,fov3;
//extern int image_width,image_height;
extern face vector a;
//char name[80];
//#include "view.h"
event acceleration (i++)
{
    vector utemp[], cellForce[], desiredForce[], markerCoord[];
    trash({faceForce});
    // 1. Get temporary velocity (advection, diffusion, pressure)
    foreach() {
        foreach_dimension() {
            //forceTotal was old
            utemp.x[] = u.x[] + dt * (g.x[] - forceTotal.x[]);
            forceTotal.x[] = 0.;
        }
    }
    for (int counter = 0; counter < Ni; counter++) { 
        setRefBC({utemp});
        // 2. calculate the force at the marker point
        foreach() {
            coord markerVelocity = {0}, desiredVelocity, markerPoint;
            int inter = ibm[] > 1e-6 && ibm[] < 1-1e-6;//changed from 0&1 so no marker points on super small fractional cells?
	        if (inter) {
                marker_point (point, ibm, &markerPoint);
	        }
	        if (inter || empty_neighbor(point, &markerPoint, ibm)){
                // interpolate to find velocity at marker point
                double markerdv = ibmdv();
		        foreach_neighbor(){
		            coord sPoint = {x,y,z};
		            double delta_u = delta_func(sPoint,markerPoint,markerdv,Delta);
                    foreach_dimension(){
                        markerVelocity.x += utemp.x[] * delta_u * ibmdv();
		            }
                }
                // calculate the desired force at the marker point
                desiredVelocity = vc;
                foreach_dimension() {
                    desiredForce.x[] = (desiredVelocity.x - markerVelocity.x) / (dt);
                    markerCoord.x[] = markerPoint.x;
                }
            }
            else{
                foreach_dimension(){
                    desiredForce.x[] = markerCoord.x[] = 0.;
		        }
	        }
	    }
        setRefBC({desiredForce,markerCoord});
        // 3. spread the force at the marker point to the nearby cell centers
        foreach() {
            coord forceSum = {0};
            if (level == max_level) {
                double markerdv = ibmdv();
		        coord sPoint = {x,y,z};
                foreach_neighbor(){
                    if (markerCoord.x[] && level == max_level) {
			            coord mcord = {markerCoord.x[],markerCoord.y[],markerCoord.z[]}; 
		                double delta_h = delta_func(sPoint,mcord,markerdv,Delta);
                        foreach_dimension() {
                            forceSum.x += (desiredForce.x[] * delta_h * ibmdv());
                        }
                    }
		        }
            }
            foreach_dimension(){
                cellForce.x[] = forceSum.x;
	        } 
        }
        setRefBC({cellForce});
        foreach(){
            foreach_dimension() {
                forceTotal.x[] += cellForce.x[];
                utemp.x[] += dt*cellForce.x[];
            }
	    }
    }
    setRefBC({forceTotal});
    // 4. correct interfacial velocity
    foreach_face(){
        faceForce.x[] = face_value (forceTotal.x, 0) + a.x[];//apply forces to acceleration field 
    }
    a = faceForce;
}

//  g is used to find uf t+dt/2 at the next time step, so the contributions
//  from f (stored in a) should be subtracted
event end_timestep (i++)
{
    //trash({a});
    //centered_gradient (p, g);
    //trash ({velocityGrad});
    //foreach(){
    //    foreach_dimension(){
    //        velocityGrad.x[] = (u.x[1] - u.x[-1])/(2.*Delta);
	//    }
    //}
}

//coord ibm_force ()
//{
//    coord ibmForce = {0};
//    foreach(reduction(+:ibmForce)){
//        foreach_dimension(){
//            ibmForce.x += -forceTotal.x[]*ibmdv();
//	    }
//    }
//    return ibmForce;
//}
//
//double ibm_pressure (Point point, scalar ibm, scalar pressure, coord normal, coord markerPoint)
//{
//    return extrapolate_scalar (point, ibm, markerPoint, normal, pressure);
//}
//
//double ibm_vorticity (Point point, vector u, coord p, coord n) // needs improvement
//{
//    coord dudn = ibm_gradientv2 (point, u, p, n);
//
//    return dudn.y*n.x - dudn.x*n.y;
//}

