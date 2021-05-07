#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>
#include<unistd.h>

/*
This program will generate a datafile of [1,4] benzothiazino [3,2-i] phenothiazine structure, which is a beta nucleating agent in isotactic polypropylene.

Structure: Moniclinic
Symmetry: P2(1)/c,
Symmetry operations:
	Molecule1 : (x, y, z)
	Molecule2 : (-x, -y, -z)
	Molecule3 : (-x, 1/2+y, 1/2-z)
	Molecule4 : (x, 1/2-y, 1/2+z)

Source for symmetry operations: http://homepage.univie.ac.at/nikos.pinotsis/spacegroup.html

Approximations: All united atoms are taken as either CH or CH2 in general, since the purpose of this crystal is to expose the cyclohexane side of the crystal to adsorbed isotactic polypropylene chains
*/

int main(int argc, char const *argv[])
{
	if(argc > 1)
	{
		// Processing argc and argv
		if(strstr(argv[1], "--help"))
		{
            printf("\n[*] Input:\n    ~~~~~~\n\n\tNONE\n\n[*] Output:\n    ~~~~~~~\n\n\t\"C18H10N2S2.data\"\n\n[*] Use this program to create TPDT-crystal (C18H10N2S2), which is in LAMMPS data file format. Output from this program can be used as input in LAMMPS script with minimal modification.\n\n");
			exit(1);
		}
		else
		{
			printf("\nError:\n~~~~~~\n\n\tUnknown command passed. Type '--help' for more information or run the executable directly\n\n");
			exit(1);
		}
	}


	int natoms, nbonds, nangles, ndihedrals, nimpropers, natomtypes, nbondtypes, nangletypes, ndihedraltypes, nimpropertypes;
	int atomid=1, moltype;
	int xnumber, ynumber, znumber;
	int fr1, fr2, fr3;

	float xmult=0, ymult=0, zmult=0;
	float a=10.472, b=5.430, c=12.356; // ( b / 2 )=2.715; ( c / 2 )=6.178;
	float halfb=2.715, halfc=6.178;
	float xlo, xhi, ylo, yhi, zlo, zhi;

	FILE *output;
	output=fopen("C18H10N2S2.data", "w");

	// Real coordinates of crystal structure based on work published by Jozef Garbarczyk
	float c1ax=0.7285*a, c2ax=0.7788*a, c3ax=0.7271*a, c4ax=0.6292*a, c5ax=0.5790*a, c6ax=0.6278*a, nx=0.8769*a, sx=0.7908*a, c1bx=0.9054*a, c2bx=0.9311*a, c3bx=0.9713*a;
	float c1ay=-0.020*b, c2ay=-0.030*b, c3ay=-0.209*b, c4ay=-0.371*b, c5ay=-0.355*b, c6ay=-0.180*b, ny=0.128*b, sy=0.1867*b, c1by=0.361*b, c2by=0.304*b, c3by=0.547*b; 
	float c1az=0.4678*c, c2az=0.3741*c, c3az=0.2905*c, c4az=0.2992*c, c5az=0.3904*c, c6az=0.4738*c, nz=0.3567*c, sz=0.5794*c, c1bz=0.5347*c, c2bz=0.4272*c, c3bz=0.6012*c;

	// Final coordinates for four molecule within the unit cell
	/*float c1ax_m1, c2ax_m1, c3ax_m1, c4ax_m1, c5ax_m1, c6ax_m1, nx_m1, sx_m1, c1bx_m1, c2bx_m1, c3bx_m1;
	float c1ay_m1, c2ay_m1, c3ay_m1, c4ay_m1, c5ay_m1, c6ay_m1, ny_m1, sy_m1, c1by_m1, c2by_m1, c3by_m1;
	float c1az_m1, c2az_m1, c3az_m1, c4az_m1, c5az_m1, c6az_m1, nz_m1, sz_m1, c1bz_m1, c2bz_m1, c3bz_m1;

	float c1ax_m2, c2ax_m2, c3ax_m2, c4ax_m2, c5ax_m2, c6ax_m2, nx_m2, sx_m2, c1bx_m2, c2bx_m2, c3bx_m2;
	float c1ay_m2, c2ay_m2, c3ay_m2, c4ay_m2, c5ay_m2, c6ay_m2, ny_m2, sy_m2, c1by_m2, c2by_m2, c3by_m2;
	float c1az_m2, c2az_m2, c3az_m2, c4az_m2, c5az_m2, c6az_m2, nz_m2, sz_m2, c1bz_m2, c2bz_m2, c3bz_m2;

	float c1ax_m3, c2ax_m3, c3ax_m3, c4ax_m3, c5ax_m3, c6ax_m3, nx_m3, sx_m3, c1bx_m3, c2bx_m3, c3bx_m3;
	float c1ay_m3, c2ay_m3, c3ay_m3, c4ay_m3, c5ay_m3, c6ay_m3, ny_m3, sy_m3, c1by_m3, c2by_m3, c3by_m3;
	float c1az_m3, c2az_m3, c3az_m3, c4az_m3, c5az_m3, c6az_m3, nz_m3, sz_m3, c1bz_m3, c2bz_m3, c3bz_m3;

	float c1ax_m4, c2ax_m4, c3ax_m4, c4ax_m4, c5ax_m4, c6ax_m4, nx_m4, sx_m4, c1bx_m4, c2bx_m4, c3bx_m4;
	float c1ay_m4, c2ay_m4, c3ay_m4, c4ay_m4, c5ay_m4, c6ay_m4, ny_m4, sy_m4, c1by_m4, c2by_m4, c3by_m4;
	float c1az_m4, c2az_m4, c3az_m4, c4az_m4, c5az_m4, c6az_m4, nz_m4, sz_m4, c1bz_m4, c2bz_m4, c3bz_m4;*/

	// Applying symmetry operations from initial fractional coordinates
	// Molecule 1
	float c1ax_m1=c1ax, c2ax_m1=c2ax, c3ax_m1=c3ax, c4ax_m1=c4ax, c5ax_m1=c5ax, c6ax_m1=c6ax, nx_m1=nx, sx_m1=sx, c1bx_m1=c1bx, c2bx_m1=c2bx, c3bx_m1=c3bx;

	float c1ay_m1=c1ay, c2ay_m1=c2ay, c3ay_m1=c3ay, c4ay_m1=c4ay, c5ay_m1=c5ay, c6ay_m1=c6ay, ny_m1=ny, sy_m1=sy, c1by_m1=c1by, c2by_m1=c2by, c3by_m1=c3by;

	float c1az_m1=c1az, c2az_m1=c2az, c3az_m1=c3az, c4az_m1=c4az, c5az_m1=c5az, c6az_m1=c6az, nz_m1=nz, sz_m1=sz, c1bz_m1=c1bz, c2bz_m1=c2bz, c3bz_m1=c3bz;

	// Molecule 2
	float c1ax_m2=-c1ax, c2ax_m2=-c2ax, c3ax_m2=-c3ax, c4ax_m2=-c4ax, c5ax_m2=-c5ax, c6ax_m2=-c6ax, nx_m2=-nx, sx_m2=-sx, c1bx_m2=-c1bx, c2bx_m2=-c2bx, c3bx_m2=-c3bx;

	float c1ay_m2=-c1ay, c2ay_m2=-c2ay, c3ay_m2=-c3ay, c4ay_m2=-c4ay, c5ay_m2=-c5ay, c6ay_m2=-c6ay, ny_m2=-ny, sy_m2=-sy, c1by_m2=-c1by, c2by_m2=-c2by, c3by_m2=-c3by;

	float c1az_m2=-c1az, c2az_m2=-c2az, c3az_m2=-c3az, c4az_m2=-c4az, c5az_m2=-c5az, c6az_m2=-c6az, nz_m2=-nz, sz_m2=-sz, c1bz_m2=-c1bz, c2bz_m2=-c2bz, c3bz_m2=-c3bz;

	// Molecule 3
	float c1ax_m3=-c1ax, c2ax_m3=-c2ax, c3ax_m3=-c3ax, c4ax_m3=-c4ax, c5ax_m3=-c5ax, c6ax_m3=-c6ax, nx_m3=-nx, sx_m3=-sx, c1bx_m3=-c1bx, c2bx_m3=-c2bx, c3bx_m3=-c3bx;

	float c1ay_m3=halfb+c1ay, c2ay_m3=halfb+c2ay, c3ay_m3=halfb+c3ay, c4ay_m3=halfb+c4ay, c5ay_m3=halfb+c5ay, c6ay_m3=halfb+c6ay, ny_m3=halfb+ny, sy_m3=halfb+sy, c1by_m3=halfb+c1by, c2by_m3=halfb+c2by, c3by_m3=halfb+c3by;

	float c1az_m3=halfc-c1az, c2az_m3=halfc-c2az, c3az_m3=halfc-c3az, c4az_m3=halfc-c4az, c5az_m3=halfc-c5az, c6az_m3=halfc-c6az, nz_m3=halfc-nz, sz_m3=halfc-sz, c1bz_m3=halfc-c1bz, c2bz_m3=halfc-c2bz, c3bz_m3=halfc-c3bz;

	// Molecule 4
	float c1ax_m4=c1ax, c2ax_m4=c2ax, c3ax_m4=c3ax, c4ax_m4=c4ax, c5ax_m4=c5ax, c6ax_m4=c6ax, nx_m4=nx, sx_m4=sx, c1bx_m4=c1bx, c2bx_m4=c2bx, c3bx_m4=c3bx;

	float c1ay_m4=halfb-c1ay, c2ay_m4=halfb-c2ay, c3ay_m4=halfb-c3ay, c4ay_m4=halfb-c4ay, c5ay_m4=halfb-c5ay, c6ay_m4=halfb-c6ay, ny_m4=halfb-ny, sy_m4=halfb-sy, c1by_m4=halfb-c1by, c2by_m4=halfb-c2by, c3by_m4=halfb-c3by;

	float c1az_m4=halfc+c1az, c2az_m4=halfc+c2az, c3az_m4=halfc+c3az, c4az_m4=halfc+c4az, c5az_m4=halfc+c5az, c6az_m4=halfc+c6az, nz_m4=halfc+nz, sz_m4=halfc+sz, c1bz_m4=halfc+c1bz, c2bz_m4=halfc+c2bz, c3bz_m4=halfc+c3bz;

	// Tilt angle for c axis (in degrees)
	float anglebeta=17.6;
	// Converting degrees to radians
	anglebeta=anglebeta/57.2958;
	float sine, cosine;
	sine=sinf(anglebeta);
	cosine=cosf(anglebeta);

	// ________________________________________________________________________________________________________________________________

	// Gathering data from user
	printf("\nEnter the size of required crystal :-");
	printf("\n\txlo: "); scanf("%f", &xlo);
	printf("\txhi: "); scanf("%f", &xhi);
	printf("\n\tylo: "); scanf("%f", &ylo);
	printf("\tyhi: "); scanf("%f", &yhi);
	printf("\n\tzlo: "); scanf("%f", &zlo);
	printf("\tzhi: "); scanf("%f", &zhi);

	printf("\nComputing crystal size...");
	/*xnumber=(xhi-xlo)/a;*/ ynumber=(yhi-ylo)/b; znumber=(zhi-zlo)/c;
	xnumber=3;// ynumber=30; znumber=30;

	printf("\nCrystal size:\n\tX-Length: %.3f\n\tY-Length: %.3f\n\tZ-Length: %.3f", (xnumber*a), (ynumber*b), (znumber*c));
	printf("\nNumber of unit cells along:\n\tX-axis: %d\n\tY-axis: %d\n\tZ-axis: %d", xnumber, ynumber, znumber);

	// Computing natom, nbonds, nangles, ndihedrals, nimpropers.
	printf("\nCreating atoms...");
	natoms=xnumber*ynumber*znumber*44; natomtypes=5;
	printf("\nIgnoring bonds, angles, dihedrals and impropers...");
	nbonds=0; nangles=0; ndihedrals=0; nimpropers=0;
	nbondtypes=0; nangletypes=0; ndihedraltypes=0; nimpropertypes=0;
	printf("\nEnter molecule ID: "); scanf("%d", &moltype);
	printf("\nEnter the starting atom index: "); scanf("%d", &atomid);

	// Printing data file
	// Header
	fprintf(output, "Created by you v1.8.1 on today, this month, this year, current time.");

	// Data file information
	fprintf(output, "\n\n\t%d\tatoms\n\t%d\tbonds\n\t%d\tangles\n\t%d\tdihedrals\n\t%d\timpropers\n\n\t%d atom types\n\t%d bond types\n\t%d angle types\n\t%d dihedral types\n\t%d improper types\n\n\t%.2f\t%.2f\txlo xhi\n\t%.2f\t%.2f\tylo yhi\n\t%.2f\t%.2f\tzlo zhi\n\nMasses\n\n\t1\t13.0907	#CG311 CH\n\t2\t14.1707	#CG321 CH2\n\t3\t32.065	# S\n\t4\t14.0067	# N\n\t5\t12.0107	# C\n\nAtoms\n\n", natoms, nbonds, nangles, ndihedrals, nimpropers, natomtypes, nbondtypes, nangletypes, ndihedraltypes, nimpropertypes, xlo, xhi, ylo, yhi, zlo, zhi);

	// Printing coordinates

	for(fr1=0; fr1<xnumber; fr1++)
	{
		for(fr2=0; fr2<ynumber; fr2++)
		{
			for(fr3=0; fr3<znumber; fr3++)
			{
				// Molecule 1
				fprintf(output, "\t%d\t%d\t1\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c1ax_m1+(xmult*a))-((c1az_m1+(zmult*c))*sine), (c1ay_m1+(ymult*b)), ((c1az_m1+(zmult*c))*cosine)); atomid++;
				fprintf(output, "\t%d\t%d\t1\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c2ax_m1+(xmult*a))-((c2az_m1+(zmult*c))*sine), (c2ay_m1+(ymult*b)), ((c2az_m1+(zmult*c))*cosine)); atomid++;
				fprintf(output, "\t%d\t%d\t2\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c3ax_m1+(xmult*a))-((c3az_m1+(zmult*c))*sine), (c3ay_m1+(ymult*b)), ((c3az_m1+(zmult*c))*cosine)); atomid++;
				fprintf(output, "\t%d\t%d\t2\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c4ax_m1+(xmult*a))-((c4az_m1+(zmult*c))*sine), (c4ay_m1+(ymult*b)), ((c4az_m1+(zmult*c))*cosine)); atomid++;
				fprintf(output, "\t%d\t%d\t2\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c5ax_m1+(xmult*a))-((c5az_m1+(zmult*c))*sine), (c5ay_m1+(ymult*b)), ((c5az_m1+(zmult*c))*cosine)); atomid++;
				fprintf(output, "\t%d\t%d\t2\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c6ax_m1+(xmult*a))-((c6az_m1+(zmult*c))*sine), (c6ay_m1+(ymult*b)), ((c6az_m1+(zmult*c))*cosine)); atomid++;

				fprintf(output, "\t%d\t%d\t4\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (nx_m1+(xmult*a))-((nz_m1+(zmult*c))*sine), (ny_m1+(ymult*b)), ((nz_m1+(zmult*c))*cosine)); atomid++;
				fprintf(output, "\t%d\t%d\t3\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (sx_m1+(xmult*a))-((sz_m1+(zmult*c))*sine), (sy_m1+(ymult*b)), ((sz_m1+(zmult*c))*cosine)); atomid++;

				fprintf(output, "\t%d\t%d\t5\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c1bx_m1+(xmult*a))-((c1bz_m1+(zmult*c))*sine), (c1by_m1+(ymult*b)), ((c1bz_m1+(zmult*c))*cosine)); atomid++;
				fprintf(output, "\t%d\t%d\t5\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c2bx_m1+(xmult*a))-((c2bz_m1+(zmult*c))*sine), (c2by_m1+(ymult*b)), ((c2bz_m1+(zmult*c))*cosine)); atomid++;
				fprintf(output, "\t%d\t%d\t1\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c3bx_m1+(xmult*a))-((c3bz_m1+(zmult*c))*sine), (c3by_m1+(ymult*b)), ((c3bz_m1+(zmult*c))*cosine)); atomid++;

				// Molecule 2
				fprintf(output, "\t%d\t%d\t1\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c1ax_m2+(xmult*a))-((c1az_m2+(zmult*c))*sine), (c1ay_m2+(ymult*b)), ((c1az_m2+(zmult*c))*cosine)); atomid++;
				fprintf(output, "\t%d\t%d\t1\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c2ax_m2+(xmult*a))-((c2az_m2+(zmult*c))*sine), (c2ay_m2+(ymult*b)), ((c2az_m2+(zmult*c))*cosine)); atomid++;
				fprintf(output, "\t%d\t%d\t2\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c3ax_m2+(xmult*a))-((c3az_m2+(zmult*c))*sine), (c3ay_m2+(ymult*b)), ((c3az_m2+(zmult*c))*cosine)); atomid++;
				fprintf(output, "\t%d\t%d\t2\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c4ax_m2+(xmult*a))-((c4az_m2+(zmult*c))*sine), (c4ay_m2+(ymult*b)), ((c4az_m2+(zmult*c))*cosine)); atomid++;
				fprintf(output, "\t%d\t%d\t2\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c5ax_m2+(xmult*a))-((c5az_m2+(zmult*c))*sine), (c5ay_m2+(ymult*b)), ((c5az_m2+(zmult*c))*cosine)); atomid++;
				fprintf(output, "\t%d\t%d\t2\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c6ax_m2+(xmult*a))-((c6az_m2+(zmult*c))*sine), (c6ay_m2+(ymult*b)), ((c6az_m2+(zmult*c))*cosine)); atomid++;

				fprintf(output, "\t%d\t%d\t4\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (nx_m2+(xmult*a))-((nz_m2+(zmult*c))*sine), (ny_m2+(ymult*b)), ((nz_m2+(zmult*c))*cosine)); atomid++;
				fprintf(output, "\t%d\t%d\t3\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (sx_m2+(xmult*a))-((sz_m2+(zmult*c))*sine), (sy_m2+(ymult*b)), ((sz_m2+(zmult*c))*cosine)); atomid++;

				fprintf(output, "\t%d\t%d\t5\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c1bx_m2+(xmult*a))-((c1bz_m2+(zmult*c))*sine), (c1by_m2+(ymult*b)), ((c1bz_m2+(zmult*c))*cosine)); atomid++;
				fprintf(output, "\t%d\t%d\t5\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c2bx_m2+(xmult*a))-((c2bz_m2+(zmult*c))*sine), (c2by_m2+(ymult*b)), ((c2bz_m2+(zmult*c))*cosine)); atomid++;
				fprintf(output, "\t%d\t%d\t1\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c3bx_m2+(xmult*a))-((c3bz_m2+(zmult*c))*sine), (c3by_m2+(ymult*b)), ((c3bz_m2+(zmult*c))*cosine)); atomid++;

				// Molecule 3
				fprintf(output, "\t%d\t%d\t1\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c1ax_m3+(xmult*a))-((c1az_m3+(zmult*c))*sine), (c1ay_m3+(ymult*b)), ((c1az_m3+(zmult*c))*cosine)); atomid++;
				fprintf(output, "\t%d\t%d\t1\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c2ax_m3+(xmult*a))-((c2az_m3+(zmult*c))*sine), (c2ay_m3+(ymult*b)), ((c2az_m3+(zmult*c))*cosine)); atomid++;
				fprintf(output, "\t%d\t%d\t2\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c3ax_m3+(xmult*a))-((c3az_m3+(zmult*c))*sine), (c3ay_m3+(ymult*b)), ((c3az_m3+(zmult*c))*cosine)); atomid++;
				fprintf(output, "\t%d\t%d\t2\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c4ax_m3+(xmult*a))-((c4az_m3+(zmult*c))*sine), (c4ay_m3+(ymult*b)), ((c4az_m3+(zmult*c))*cosine)); atomid++;
				fprintf(output, "\t%d\t%d\t2\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c5ax_m3+(xmult*a))-((c5az_m3+(zmult*c))*sine), (c5ay_m3+(ymult*b)), ((c5az_m3+(zmult*c))*cosine)); atomid++;
				fprintf(output, "\t%d\t%d\t2\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c6ax_m3+(xmult*a))-((c6az_m3+(zmult*c))*sine), (c6ay_m3+(ymult*b)), ((c6az_m3+(zmult*c))*cosine)); atomid++;

				fprintf(output, "\t%d\t%d\t4\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (nx_m3+(xmult*a))-((nz_m3+(zmult*c))*sine), (ny_m3+(ymult*b)), ((nz_m3+(zmult*c))*cosine)); atomid++;
				fprintf(output, "\t%d\t%d\t3\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (sx_m3+(xmult*a))-((sz_m3+(zmult*c))*sine), (sy_m3+(ymult*b)), ((sz_m3+(zmult*c))*cosine)); atomid++;

				fprintf(output, "\t%d\t%d\t5\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c1bx_m3+(xmult*a))-((c1bz_m3+(zmult*c))*sine), (c1by_m3+(ymult*b)), ((c1bz_m3+(zmult*c))*cosine)); atomid++;
				fprintf(output, "\t%d\t%d\t5\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c2bx_m3+(xmult*a))-((c2bz_m3+(zmult*c))*sine), (c2by_m3+(ymult*b)), ((c2bz_m3+(zmult*c))*cosine)); atomid++;
				fprintf(output, "\t%d\t%d\t1\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c3bx_m3+(xmult*a))-((c3bz_m3+(zmult*c))*sine), (c3by_m3+(ymult*b)), ((c3bz_m3+(zmult*c))*cosine)); atomid++;

				// Molecule 4
				fprintf(output, "\t%d\t%d\t1\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c1ax_m4+(xmult*a))-((c1az_m4+(zmult*c))*sine), (c1ay_m4+(ymult*b)), ((c1az_m4+(zmult*c))*cosine)); atomid++;
				fprintf(output, "\t%d\t%d\t1\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c2ax_m4+(xmult*a))-((c2az_m4+(zmult*c))*sine), (c2ay_m4+(ymult*b)), ((c2az_m4+(zmult*c))*cosine)); atomid++;
				fprintf(output, "\t%d\t%d\t2\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c3ax_m4+(xmult*a))-((c3az_m4+(zmult*c))*sine), (c3ay_m4+(ymult*b)), ((c3az_m4+(zmult*c))*cosine)); atomid++;
				fprintf(output, "\t%d\t%d\t2\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c4ax_m4+(xmult*a))-((c4az_m4+(zmult*c))*sine), (c4ay_m4+(ymult*b)), ((c4az_m4+(zmult*c))*cosine)); atomid++;
				fprintf(output, "\t%d\t%d\t2\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c5ax_m4+(xmult*a))-((c5az_m4+(zmult*c))*sine), (c5ay_m4+(ymult*b)), ((c5az_m4+(zmult*c))*cosine)); atomid++;
				fprintf(output, "\t%d\t%d\t2\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c6ax_m4+(xmult*a))-((c6az_m4+(zmult*c))*sine), (c6ay_m4+(ymult*b)), ((c6az_m4+(zmult*c))*cosine)); atomid++;

				fprintf(output, "\t%d\t%d\t4\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (nx_m4+(xmult*a))-((nz_m4+(zmult*c))*sine), (ny_m4+(ymult*b)), ((nz_m4+(zmult*c))*cosine)); atomid++;
				fprintf(output, "\t%d\t%d\t3\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (sx_m4+(xmult*a))-((sz_m4+(zmult*c))*sine), (sy_m4+(ymult*b)), ((sz_m4+(zmult*c))*cosine)); atomid++;

				fprintf(output, "\t%d\t%d\t5\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c1bx_m4+(xmult*a))-((c1bz_m4+(zmult*c))*sine), (c1by_m4+(ymult*b)), ((c1bz_m4+(zmult*c))*cosine)); atomid++;
				fprintf(output, "\t%d\t%d\t5\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c2bx_m4+(xmult*a))-((c2bz_m4+(zmult*c))*sine), (c2by_m4+(ymult*b)), ((c2bz_m4+(zmult*c))*cosine)); atomid++;
				fprintf(output, "\t%d\t%d\t1\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c3bx_m4+(xmult*a))-((c3bz_m4+(zmult*c))*sine), (c3by_m4+(ymult*b)), ((c3bz_m4+(zmult*c))*cosine)); atomid++;

				// Incrementing Z-coordinates
				zmult++;
			}
			// reset Z-coordinates; increment Y-coordinates
			zmult=0; ymult++;
		}								
		// reset Z- and Y-coordinates; increment X-coordinates
		zmult=0; ymult=0; xmult++;
	}

	printf("\n");
	fclose(output);

	return(0);

}

