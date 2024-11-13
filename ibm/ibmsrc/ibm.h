#define BGHOSTS 2


#include "fractions.h"
#include "ibm-utils.h"


extern coord vc;          // object's imposed velocity
extern scalar vof;
extern double maxlevel;

vector desiredForce[];    // force calculated at the marker point
vector cellForce[];       // force located at the cell center (after spreading)
vector markerCoord[];     // field to store the coordinates of all marker points

vector velocityGrad[];

face vector faceForce[];  // for averaging the cell force to get the face values

vector utemp[];
vector forceTotal[];

int Ni = 10;               // # of multi-direct forcing iterations


event init (t = 0)
{
    vof.refine = vof.prolongation = fraction_refine;
    vof.dirty = true;//set prolongation 
    foreach(){
        foreach_dimension() { 
            velocityGrad.x[] = 0.;
            if (vof[] == 1){
                u.x[] = vc.x;
	    }
        }
    }
}

//allows for our reflection boundary condidion to be set on a vector field
void setRefBC(vector *set){
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
}

event acceleration (i++)
{
    trash({cellForce, desiredForce, markerCoord, faceForce, utemp, forceTotal});
    // 1. Get temporary velocity (advection, diffusion, pressure)
    foreach() {
        foreach_dimension() {
            utemp.x[] = u.x[] - dt*(p[] - p[-1])/(Delta);
            forceTotal.x[] = 0.;
        }
    }
    boundary({utemp});
    setRefBC({utemp});
    for (int counter = 0; counter < Ni; counter++) { 
        // 2. calculate the force at the marker point
        foreach() {
            coord markerVelocity = {0}, desiredVelocity, markerPoint;
            int inter = vof[] > 1e-6 && vof[] < 1-1e-6;//changed from 0&1 so no marker points on super small fractional cells?
	    if (inter) {
                marker_point (point, vof, &markerPoint);
	    }
	    if (inter || empty_neighbor(point, &markerPoint, vof)){
                // interpolate to find velocity at marker point
                double markerdv = dv();
		foreach_neighbor(){
		    coord sPoint = {x,y,z};
		    double delta_u = delta_func(sPoint,markerPoint,markerdv,Delta);
                    foreach_dimension(){
                        markerVelocity.x += utemp.x[] * delta_u * dv();
		    }
                }
                // calculate the desired force at the marker point
                desiredVelocity = vc;
                foreach_dimension() {
                    desiredForce.x[] = (desiredVelocity.x - markerVelocity.x) / dt;
                    markerCoord.x[] = markerPoint.x;
                }
            }
            else{
                foreach_dimension(){
                    desiredForce.x[] = markerCoord.x[] = 0.;
		}
	    }
	}
	boundary({desiredForce,markerCoord});
        setRefBC({desiredForce,markerCoord});
        // 3. spread the force at the marker point to the nearby cell centers
        foreach() {
            coord forceSum = {0};
            if (level == maxlevel) {
                double markerdv = dv();
		coord sPoint = {x,y,z};
                foreach_neighbor(){
                    if (markerCoord.x[] && level == maxlevel) {
			coord mcord = {markerCoord.x[],markerCoord.y[],markerCoord.z[]}; 
		        double delta_h = delta_func(sPoint,mcord,markerdv,Delta);
                        foreach_dimension() {
                            forceSum.x += (desiredForce.x[] * delta_h * dv());
                        }
                    }
		}
            }
            foreach_dimension(){
                cellForce.x[] = forceSum.x;
	    } 
        }
	boundary({cellForce});
        setRefBC({cellForce});
        foreach(){
            foreach_dimension() {
                forceTotal.x[] += cellForce.x[];
                utemp.x[] += dt*cellForce.x[];
            }
	}
    }
    setRefBC({forceTotal});
    boundary({forceTotal});
    // 4. correct interfacial velocity
    foreach_face(){
        faceForce.x[] = fm.x[]*(face_value (forceTotal.x, 0));
    }
    a = faceForce;

}

//  g is used to find uf t+dt/2 at the next time step, so the contributions
//  from f (stored in a) should be subtracted

event end_timestep (i++)
{
    trash({a});
    centered_gradient (p, g);
    trash ({velocityGrad});
    foreach(){
        foreach_dimension(){
            velocityGrad.x[] = (u.x[1] - u.x[-1])/(2.*Delta);
	}
    }
}

coord ibm_force ()
{
    coord ibmForce = {0};
    foreach(reduction(+:ibmForce)){
        foreach_dimension(){
            ibmForce.x += -forceTotal.x[]*dv();
	}
    }
    return ibmForce;
}

double ibm_pressure (Point point, scalar vof, scalar pressure, coord normal, coord markerPoint)
{
    return extrapolate_scalar (point, vof, markerPoint, normal, pressure);
}

double ibm_vorticity (Point point, vector u, coord p, coord n) // needs improvement
{
    coord dudn = ibm_gradientv2 (point, u, p, n);

    return dudn.y*n.x - dudn.x*n.y;
}

