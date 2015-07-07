#define ge(x,y) (x>y)?1:0  /* indicator for [0,1)  */

#define THRESHOLD1 100   /* Not currently used */

#define THRESHOLD2 1.2  /* Not currently used */

#define THRESHOLD3 15   /* The THRESHOLD3 highest 'performing' regions will be */
                        /* indicated  -  This is the default for `threshold' in */
                        /* the program. Can be set  differently by one of the */
                        /* command line arguments. */

#define G(n,allmean,allvar,mean,var,grad,lap,totgrad,numch,thresh1,thresh2) totgrad

#define G1(n,allmean,allvar,mean,var,grad,lap,sqlap,numch,thresh1,thresh2) numch
#define G2(n,allmean,allvar,mean,var,grad,lap,sqlap,numch,thresh1,thresh2) numch
#define G3(n,allmean,allvar,mean,var,grad,lap,sqlap,numch,thresh1,thresh2) numch
#define G4(n,allmean,allvar,mean,var,grad,lap,sqlap,numch,thresh1,thresh2) numch
#define G5(n,allmean,allvar,mean,var,grad,lap,sqlap,numch,thresh1,thresh2) numch
#define G6(n,allmean,allvar,mean,var,grad,lap,sqlap,numch,thresh1,thresh2) numch
#define G7(n,allmean,allvar,mean,var,grad,lap,sqlap,numch,thresh1,thresh2) numch
#define G8(n,allmean,allvar,mean,var,grad,lap,sqlap,numch,thresh1,thresh2) numch
#define G9(n,allmean,allvar,mean,var,grad,lap,sqlap,numch,thresh1,thresh2) numch
#define G10(n,allmean,allvar,mean,var,grad,lap,sqlap,numch,thresh1,thresh2) numch

#define F(n,allmean,allvar,mean,var,grad,lap,sqlap,numch,thresh1,thresh2) ge(G(n,allmean,allvar,mean,var,grad,lap,sqlap,numch,thresh1,thresh2),thresh1)

