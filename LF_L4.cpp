#include"iostream"
#include"nr3.h"
#include"fstream"
#include"svd.h"
#include"ran.h"

using namespace std; 

/***************************************************/
/*  This defines the motion via the 1st order      */
/*  leapfrog (Verlet) integration of               */
/*  a logarithmic potential in the rotating frame  */
/*  Phi = vc^2 ln (R) = 0.5 vc^2 ln (x^2 + y^2)    */
/*  Phi_eff = 0.5 vc^2 ln (x^2 + y^2)              */
/*            - vc^2/(2 CR^2) *(x^2+y^2)           */
/*  Uses (R,phi) coordinates for calculation       */
/***************************************************/

/*********************************/
/* Declare major steps necessary */
/*********************************/
  // Routine to take one step
VecDoub leapstep(int ndim, VecDoub qp, double dt, double tnow, double CR, float thetain, float epsilon, int m);
  // Acceleration for this potential
VecDoub accel(int ndim, VecDoub qp, double tnow, double CR, float thetain, float epsilon, int m);

/***********************************************/
/* Main program that controls output and steps */
/***********************************************/
int main() 
{
    
    float a, b, c, d;

    ifstream table;
    table.open("./temp_initials.txt");

    table >> a;
    table >> b;
    table >> c;
    table >> d;

    cout << a << " " << b << " " << c << " " << d << endl;

    table.close();
   
 /* Declare input variables */
   int ndim = 2; // Number of dimentions (x,y)
   VecDoub qp(ndim+ndim);
   double tnow=0.;
   
 /* Define physical constants */
   const Doub km2kpc=3.24e-17; // Conversion from km to kpc
   const Doub yr2sec=3.156e7; // Conversion from km to kpc
   const Doub PI=3.141592; // Define PI
   const Doub vc=220.0*km2kpc*yr2sec; // Circular velocity [kpc/yr]
   
 /* Initial conditions */
   tnow = 0.; // Set initial time [yrs]
 // Initial position in (x,y) coords
   float x0;
   x0 = a;
   qp[0] = x0; // x(t=0) [kpc]
   float y0;
   y0 = b;
   qp[1] = y0; // y(t=0) [kpc]
 // Initial velocity in the inertial frame
   float vx0;
   vx0 = c;
   qp[2] = vx0*km2kpc*yr2sec; // v_x0 = dx/dt|_(t=0) [kpc/yr]
   float vy0;
   vy0 = d;
   qp[3] = vy0*km2kpc*yr2sec; // v_y0 = dy/dt|_(t=0) [kpc/yr]
 // Spiral description
   double CR;
   CR = 8.0;
   float thetain;
   thetain = 25.0;
   float epsilon;
   epsilon = 0.3;
   int m;
   m = 4;

 /* Set intigration parameters */
   float ttot = 2e9; // Total time for integration [yrs]
   double dt = 1.e2; // Timestep for integration [yrs]
   int mstep = round(ttot/dt); // Number of steps to take
   int outputnumber = 2000; // Number of outputs
   int nout = floor(mstep/outputnumber); // Output density=number of steps between outputs
   
 /* Open file and write solution to it */
   stringstream ssdatfile (stringstream::in | stringstream::out);
   ssdatfile << "./qp_file/qp_(m=" << m << ")_(th=" << thetain << ")_(t=" << ttot/1e9 << ")_(CR=" << CR << ")_(eps="<< epsilon<<")_(x0=" << x0 << ")_(y0=" << y0 << ")_(vx0=" << vx0 << ")_(vy0=" << vy0 << ").txt";
   string datfile = ssdatfile.str();
   cout << datfile << '\n';
   
   ofstream file;
   file.open( datfile.c_str() );
   
 /* Integration loop */
   for (int nstep=0; nstep<mstep; nstep++) {
      /* Convert to (R,phi) coordinates in the rotating frame */
      double R, phi;
      R = pow( pow(qp[0],2) + pow(qp[1],2), 0.5);
      phi = atan(qp[1]/qp[0]) - tnow*vc/CR;
      if (qp[0] < 0  && qp[1] >= 0) phi = atan(qp[1]/qp[0]) + PI - tnow*vc/CR;
      if (qp[0] < 0  && qp[1] < 0)  phi = atan(qp[1]/qp[0]) - PI - tnow*vc/CR;
      if (qp[0] == 0 && qp[1] > 0)  phi = PI/2. - tnow*vc/CR;
      if (qp[0] == 0 && qp[1] < 0)  phi = -PI - tnow*vc/CR;
      if (qp[0] == 0 && qp[1] == 0) phi = - tnow*vc/CR;
      
      // if this is an output step print current values via the print routine
      if (nstep % nout == 0)
         file << tnow << " " << qp[0] << " " << qp[1] << " " <<
             qp[2] << " " << qp[3] << " " <<
             R << " " << phi << '\n';
   
      qp=leapstep(ndim,qp,dt,tnow,CR,thetain,epsilon,m); // take an integration step
      tnow = tnow + dt; // update the current time
   }

   // print last step
   double Rend, phiend;
   if (mstep % nout == 0)
      Rend = pow( pow(qp[0],2) + pow(qp[1],2), 0.5);
   phiend = atan(qp[1]/qp[0]) - tnow*vc/CR;
   if (qp[0] < 0  && qp[1] >= 0) phiend = atan(qp[1]/qp[0]) + PI - tnow*vc/CR;
   if (qp[0] < 0  && qp[1] < 0)  phiend = atan(qp[1]/qp[0]) - PI - tnow*vc/CR;
   if (qp[0] == 0 && qp[1] > 0)  phiend = PI/2. - tnow*vc/CR;
   if (qp[0] == 0 && qp[1] < 0)  phiend = -PI - tnow*vc/CR;
   if (qp[0] == 0 && qp[1] == 0) phiend = - tnow*vc/CR;
   
   file << tnow << " " << qp[0] << " " << qp[1] << " " <<
          qp[2] << " " << qp[3] << " " <<
          Rend << " " << phiend << '\n';


   /* Close file */
   file.close();

   return (0);
}

/***************************************************/
/* Leapfrog integrator for one step from t to t+dt */
/* NOTE: this employs the kick-drift-kick method.  */
/***************************************************/
VecDoub leapstep(int ndim, VecDoub qp, double dt, double tnow, double CR, float thetain, float epsilon, int m) {

 /* Define initial parameters */
   VecDoub a(ndim); // Define acceleration in ndim dimentions
   
   a=accel(ndim,qp,tnow,CR,thetain,epsilon,m); // call acceleration code for given position and velocity
   qp[2]=qp[2]-0.5 *dt *a[0]; // advance v_R by half step
   qp[3]=qp[3]-0.5 *dt *a[1]; // advance v_phi by half step
   qp[0]=qp[0]+dt*qp[2]; // advance x by full step
   qp[1]=qp[1]+dt*qp[3]; // advance x by full step
   a=accel(ndim,qp,tnow,CR,thetain,epsilon,m); // call ecceleration code for new coords and velocities
   qp[2]=qp[2]-0.5 *dt *a[0]; // complete v_x step
   qp[3]=qp[3]-0.5 *dt *a[1]; // complete v_y step
   
   return qp;
   
}

/*************************************************/
/* Acceleration code for a logarithmic potential */
/*  Phi = vc^2 ln (R) = 0.5 vc^2 ln (x^2 + y^2)  */
/*  a_x = vc^2 *x /(x^2 + y^2)                   */
/*  a_y = vc^2 *y /(x^2 + y^2)                   */
/*************************************************/
VecDoub accel(int ndim, VecDoub qp, double tnow, double CR, float thetain, float epsilon, int m) {
 /* Define physical constants */
   const Doub km2kpc = 3.24e-17; // Conversion from km to kpc
   const Doub yr2sec = 3.156e7; // Conversion from km to kpc
   const Doub pc2kpc = 1e-3; // Conversion from pc to kpc
   const Doub PI = 3.141592; // Define PI
   const Doub G = 4.498e-24; // Gravitational constant [kpc^3 Msun^-1 yr^-2]
   const Doub Rd = 2.5 ; // Disk scale length [kpc]
   const Doub Rsun = 8.0 ; // Distance to Sun [kpc]
   const Doub vc = 220.0*km2kpc*yr2sec; // Circular velocity [kpc/yr]
   
   /* Pitch angle enterd in degrees [rad] */
   const Doub theta = thetain *2.*PI/360. ;
   
   const Doub R= pow( pow(qp[0],2) + pow(qp[1],2), 0.5);
     
   const Doub Sig0 = 50. *exp(Rsun/Rd)/pow(pc2kpc,2) ; // Scale for surface density [Msun kpc^-2]
   const Doub SigCR = Sig0 *exp(-CR/Rd) ; // Surface brightness at CR [Msun kpc^-2]
   const Doub SigR = Sig0 *exp(-R/Rd) ; // Surface brightness at CR [Msun kpc^-2]
   const Doub alpha = m /tan(theta) ;
   const Doub A = 2.*PI*G*SigR*epsilon*R/alpha ;
   
   /* Define accelation */
   VecDoub a(ndim);
   /* Acceleration from logarithmic disk */
   double aDx = pow(vc,2)*qp[0]/(pow(qp[0],2) + pow(qp[1],2));
   double aDy = pow(vc,2)*qp[1]/(pow(qp[0],2) + pow(qp[1],2));
   /* Acceleration from spiral */
   double Sstuff = -A *sin(m*tnow*vc/CR -m*atan(qp[1]/qp[0])
                        -alpha*log( pow( pow(qp[0],2)+pow(qp[1],2), 0.5)/CR))
                   /(pow(qp[0],2) + pow(qp[1],2));
   double aSx = Sstuff *(m*qp[1] - alpha*qp[0]);
   double aSy = -Sstuff *(m*qp[0] + alpha*qp[1]);
   /* Acceleration from rotating frame */
   double aRx = -qp[0]*pow(vc/CR,2);
   double aRy = -qp[1]*pow(vc/CR,2);
   
   /* Total Acceleration */
   a[0] = aDx + aSx;
   a[1] = aDy + aSy;
   
   return a;
}