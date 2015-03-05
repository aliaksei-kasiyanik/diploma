#include <stdio.h>
#include <math.h>
#include "dopri5.h"
#include <iostream>
using namespace std;
#define  ndgl      8
#define  nrdens    8
#define  licont    nrdens
char format99[] = "";//x=%f  y=%12.10f %12.10f  nstep=%li\r\n

void faren(unsigned n, double x, double *arg, double *res)
{

	res[0]=-1.71*arg[0]+0.43*arg[1]+8.32*arg[2]+0.0007;
	res[1]= 1.71*arg[0]-8.75*arg[1];
	res[2]=-10.03*arg[2]+0.43*arg[3]+0.035*arg[4];
	res[3]=8.32*arg[1]+1.71*arg[2]-1.12*arg[3];
	res[4]=-1.745*arg[4]+0.43*arg[5]+0.43*arg[6];
	res[5]=-280*arg[5]*arg[7]+0.69*arg[3]+1.71*arg[4]-0.43*arg[5]+0.69*arg[6];
	res[6]=280*arg[5]*arg[7] - 1.81*arg[6];
	res[7]=-280*arg[5]*arg[7]+1.81*arg[6];
} /* faren */

void solout(long nr, double xold, double x, double* y, unsigned n, int* irtrn)
{
	static double xout;

	if (nr == 1)
	{
		printf(format99, x, y[0], y[1], nr-1);
		xout = x + 2.0;
	}
	else
		while (x >= xout)
		{
			printf(format99, xout, contd5(0, xout), contd5(1, xout), nr-1);
			xout += 2.0;
		}

} /* solout */

double get_err(double *ans)
{
	double temp[8];
	temp[0]=0.00071664185131179282753 ;
	temp[1]=0.00014005267307968342623 ;
	temp[2]=5.0123734925483684184e-05 ;
	temp[3]=0.0011169432218426586251 ;
	temp[4]=0.00064155093941671657159 ;
	temp[5]=0.0015189353975143206257 ;
	temp[6]=0.0010845275986320803818 ;
	temp[7]=0.0046154724013679230762 ; 


	/*temp[0] = 0.0007166418506676627;
	 temp[1] = 0.0001400526737015955;
	 temp[2] = 0.000050123735559464576;
	 temp[3] = 0.001116943221152604;
	 temp[4] = 0.0006415509394182438;
	 temp[5] = 0.0015189353975998006;
	 temp[6] = 0.0010845275986068634;
	 temp[7] = 0.00461547240139314;*/
//t=321
//	temp[0] =0.00073735395740377924784;
//	temp[1] = 0.00014429214338435538052;
//	temp[2] = 5.8908769930807395178e-05;
//	temp[3] = 0.001176011091929014124;
//	temp[4] = 0.0023880728643662289638;
//	temp[5] = 0.0062444271013042907878;
//	temp[6] = 0.0028512129595190277044;
//	temp[7] = 0.002848787040481005027;

//	temp[0] = 0.00067060207437807028627;
//	temp[1] = 0.00013105477331575279202;
//	temp[2] = 4.688294207925115326e-05;
//	temp[3] = 0.0010451279952393328649;
//	temp[4] = 0.00059517071282648925049;
//	temp[5] = 0.0014003610097888124792;
//	temp[6] = 0.0010149280109548290685;
//	temp[7] = 0.0046850719890452040967;  

	//      temp[0] = 0.00073735395740377924784;
	//      temp[1] = 0.00014429214338435538052;
	//      temp[2] = 5.8908769930807395178e-05;
	//      temp[3] = 0.001176011091929014124;
	//      temp[4] = 0.0023880728643662289638;
	//      temp[5] = 0.0062444271013042907878;
	//      temp[6] =  0.0028512129595190277044;
	//      temp[7] =  0.002848787040481005027;

	for (int i =0; i<8; i++)
		temp[i] = fabs(temp[i]-ans[i])/fabs(temp[i]);
	double max = 0;
	for (int i=0; i<8; i++)
		if (temp[i] > max)
			max = temp[i];
	return max;
}

int main(void)
{
	double y[ndgl];
	unsigned icont[licont], i;
	int res, iout, itoler;
	double x, xend, atoler, rtoler;

	iout = 2;
	x = 0.0;

	for (i=0; i< 7; i++)
		y[i] = 0.;
	y[0] = 1.;
	y[7] = 0.0057;

	xend = 321.8122;
	itoler = 0;
	rtoler = 1.0E-10;
	atoler = rtoler;
	icont[0] = 0;
	icont[1] = 1;
	for (int kkk = 1; kkk<13; kkk++)
	{

		iout = 2;
		x = 0.0;

		for (i=0; i< 7; i++)
			y[i] = 0.;
		y[0] = 1.;
		y[7] = 0.0057;

		itoler = 0;

		rtoler = pow(10,-kkk);///=10.;
		atoler = rtoler;
		icont[0] = 0;
		icont[1] = 1;
		cout<<"======= experiment======="<<endl;
		cout<<"       tol="<<rtoler<<endl;

		res
				= dopri5(ndgl, faren, x, y, xend, &rtoler, &atoler, itoler, solout, iout,
				stdout, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10000000, 0, 0, ndgl, NULL, licont);

		printf("x=xend  y=%12.10f %12.10f\r\n", y[0], y[1]);
		printf("rtol=%12.10f   fcn=%li   step=%li   accpt=%li   rejct=%li\r\n",
				rtoler, nfcnRead(), nstepRead(), naccptRead(), nrejctRead());

		cout.precision(20);
		cout<<nfcnRead()<<", " <<get_err(y)<<endl;
				
		cout<<"solution"<<endl;
		  for(int j = 0; j  <8 ; j++)
			  cout<<y[j]<<"  ";
		cout<<endl<<endl<<endl<<endl;
	}
	return 0;

} /* main */

