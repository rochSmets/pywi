#include <stdio.h>
#include <stdlib.h>
#include <math.h>



int main(void)
{
	char filename[512];
	int nt, i;
	int nx;
	int nz;
	int xm;
	int zm;
	float dx, dz;
	
	long sf, si;
	int nss;
	FILE *fp, *fpres;

	long offset;
	long size2d, size2ds;

	float *bx, *by, *bz;
	float *flux;

	int ix, iz;
	int ij, ij1,ij2;

	float xe, ze;

	nx = 1000;
	nz = 800;
	xm = 320;
	zm = 128;


	nss = 2;

	dx = xm/(float)(nx-2) / (sqrt(25.));
	dz = zm/(float)(nz-2) / (sqrt(25.));

	nt = 8000;

	sf = sizeof (float);
	si = sizeof(int);
	size2d = nx*nz*sf;
	size2ds = size2d*nss;
	

	bx = calloc(nx*nz, sizeof (float));
	by = calloc(nx*nz, sizeof (float));
  	bz = calloc(nx*nz, sizeof (float));	
	flux = calloc(nx*nz, sizeof (float));	

	printf("dx = %f, dz = %f\n",dx,dz);

	for (i=50; i <=nt; i+=50)
	{
		sprintf(filename, "/Volumes/drobo/nico/asymmetric/by10/data/fields-%05d.dat",i);

		fp = fopen(filename, "rb");

		if(!fp) printf("problem file\n");

		else printf("Opening %s\n", filename);

		fseek(fp, si + si + sf*4 + 2*si + 3*size2ds, SEEK_SET);
		fread(bx , nx*nz, sf, fp);
		fread(by, nx*nz, sf, fp);
		fread(bz, nx*nz, sf, fp);
		
		fclose(fp);

		printf("test1\n");

		flux[0] = 0.0;
		iz = 0;

		for (ix=1; ix < nx; ix++)
		{
			ij  = ix + 0*nx;
			ij1 = ix-1 + 0*nx;

			flux[ij] = flux[ij1] + bz[ij1]*dx*4.;
			//printf("ix = %d\n", ix);
		}
		printf("test2\n");
		ix=0;
		for (iz=1; iz < nz; iz++)
		{
			ij = ix + iz*nx;
			ij1 = ix + (iz-1)*nx;

			flux[ij] = flux[ij1] - bx[ij1]*dz*4;
		//	printf("%f \n", bx[ij]*4);
		}
		
		printf("test3\n");

		for (ix=1; ix < nx; ix++)
		{
			for(iz=1; iz < nz; iz++)
			{
				ij  = ix   +  iz*nx;
				ij1 = ix-1 +  iz*nx;
				ij2 = ix   +  (iz-1)*nx;

				flux[ij] = 0.5*(flux[ij1] + bz[ij1]*dx*4) + 0.5*(flux[ij2] - bx[ij2]*dz*4);
			}
		}

		sprintf(filename, "/Volumes/drobo/nico/asymmetric/by10/data/flux-%05d.dat",i);
		
		printf("opening res file %s\n", filename);
		fpres = fopen(filename, "wb");

		fwrite(flux, nx*nz, sizeof *flux, fp);

		fclose(fpres);

		ij = 450 + 600*nx;
		printf("fllux[450,600] = %f, bx[ij]= %f, bz[0] = %0.9f\n",flux[ij], bx[0], bz[0]);
		
	}



	free(bx);
	free(by);
	free(bz);


	return 0;
}
