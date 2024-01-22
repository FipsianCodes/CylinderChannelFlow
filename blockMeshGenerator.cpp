/*
 * =====================================================================================
 *
 *       Filename:  blockMeshGenerator.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  23.11.2020 11:47:15
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  B. Philippi (), 
 *   Organization:  CAU - Algorithmic Optimal Control
 *
 * =====================================================================================
 */

#include <iomanip>
#include <fstream>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include </media/fips/data/parareal/openFoam/CAU/libconfig-1.7.2/lib/libconfig.h++>

using namespace libconfig;
using namespace std;

#define PI 3.14159265359 

int main(int argc, char **argv)
{
  	Config cfg;

  	try {cfg.readFile("CylinderParameters.cfg");}
  	catch(const FileIOException &fioex){
	    	std::cerr << "File I/O error" << std::endl;
	    	return(EXIT_FAILURE);
  	}

    	catch(const ParseException &pex){
        	std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
          	<< " - " << pex.getError() << std::endl;
          	return(EXIT_FAILURE);
      	}
  
      	try{
		std::string content = cfg.lookup("content");
		std::printf("\nReading %s...\n", content.c_str());
      	}
      	catch(const SettingNotFoundException &nfex){printf("\nNo content found...\n");}
 
      	std::string content = cfg.lookup("content");
  
  	try{
       		double RC = cfg.lookup("RadiusCylinder");
	        printf("\n\tRC         :: Radius Cylinder                        ::     \t\t %3.2e m  \n", RC);
      	}
      	catch(const SettingNotFoundException &nfex){
              	std::cerr << "No 'RadiusCylinder' setting in configuration file." << endl;
      	}
      	double RC = cfg.lookup("RadiusCylinder");


      	try{
       		double RR = cfg.lookup("RadiusRefinement");
	        printf("\n\tRR         :: Radius Refinement                      ::     \t\t %3.2e m  \n", RR);
      	}
      	catch(const SettingNotFoundException &nfex){
        	std::cerr << "No 'RadiusRefinement' setting in configuration file." << endl;
      	}
      	double RR = cfg.lookup("RadiusRefinement");

 
      	try{
       		double H = cfg.lookup("Height");
        	printf("\n\tH          :: Domain Height                            ::     \t\t %3.2e m  \n", H);
      	}
      	catch(const SettingNotFoundException &nfex){
              	std::cerr << "No 'Height' setting in configuration file." << endl;
      	}
      	double H = cfg.lookup("Height");

      	try{
       		double W = cfg.lookup("Width");
        	printf("\n\tW          :: Domain Width                             ::     \t\t %3.2e m  \n", W);
      	}
      	catch(const SettingNotFoundException &nfex){
                std::cerr << "No 'Width' setting in configuration file." << endl;
      	}
      	double W = cfg.lookup("Width");

      	try{
       		double D = cfg.lookup("Depth");
        	printf("\n\tD          :: Domain Depth                             ::     \t\t %3.2e m  \n", D);
      	}
      	catch(const SettingNotFoundException &nfex){
                std::cerr << "No 'Depth' setting in configuration file." << endl;
      	}
      	double D = cfg.lookup("Depth");

      	try{
       		double CCx = cfg.lookup("CenterCylinder_x");
        	printf("\n\tCCx        :: Center Cylinder x-Coordinate             ::     \t\t %3.2e m  \n", CCx);
      	}
      	catch(const SettingNotFoundException &nfex){
                std::cerr << "No 'CenterClyinder_x' setting in configuration file." << endl;
      	}
      	double CCx = cfg.lookup("CenterCylinder_x");

     	try{
       		double CCy = cfg.lookup("CenterCylinder_y");
        	printf("\n\tCCy        :: Center Cylinder y-Coordinate             ::     \t\t %3.2e m  \n", CCy);
      	}
      	catch(const SettingNotFoundException &nfex){
                std::cerr << "No 'CenterCylinder_y' setting in configuration file." << endl;
      	}
      	double CCy = cfg.lookup("CenterCylinder_y");

      	try{
       		double cTM = cfg.lookup("convertToMeters");
        	printf("\n\tcTM        :: convertToMeters factor                   ::     \t\t %3.2e \n", cTM);
      	}
      	catch(const SettingNotFoundException &nfex){
                std::cerr << "No 'convertToMeters' setting in configuration file." << endl;
      	}
      	double cTM = cfg.lookup("convertToMeters");

      	try{
       		int Nx1 = cfg.lookup("CellsInDirection_x1");
        	printf("\n\tNx1        :: # Cells in x-direction, first section    ::     \t\t %d \n", Nx1);
      	}
      	catch(const SettingNotFoundException &nfex){
                std::cerr << "No 'CellsInDirection_x1' setting in configuration file." << endl;
      	}
      	int Nx1 = cfg.lookup("CellsInDirection_x1");

      	try{
       		int Nx2 = cfg.lookup("CellsInDirection_x2");
        	printf("\n\tNx2        :: # Cells in x-direction, second section   ::     \t\t %d \n", Nx2);
      	}
      	catch(const SettingNotFoundException &nfex){
                std::cerr << "No 'CellsInDirection_x2' setting in configuration file." << endl;
      	}
      	int Nx2 = cfg.lookup("CellsInDirection_x2");

      	try{
       		int Nx3 = cfg.lookup("CellsInDirection_x3");
        	printf("\n\tNx3        :: # Cells in x-direction, third section    ::     \t\t %d \n", Nx3);
      	}
      	catch(const SettingNotFoundException &nfex){
                std::cerr << "No 'CellsInDirection_x3' setting in configuration file." << endl;
      	}
      	int Nx3 = cfg.lookup("CellsInDirection_x3");
    
      	try{
       		int Ny1 = cfg.lookup("CellsInDirection_y1");
        	printf("\n\tNy1        :: # Cells in y-direction, first section    ::     \t\t %d \n", Ny1);
      	}
      	catch(const SettingNotFoundException &nfex){
               std::cerr << "No 'CellsInDirection_y1' setting in configuration file." << endl;
      	}
      	int Ny1 = cfg.lookup("CellsInDirection_y1");

      	try{
       		int Ny2 = cfg.lookup("CellsInDirection_y2");
        	printf("\n\tNy2        :: # Cells in y-direction, second section   ::     \t\t %d \n", Ny2);
      	}
      	catch(const SettingNotFoundException &nfex){
                std::cerr << "No 'CellsInDirection_y2' setting in configuration file." << endl;
      	}
      	int Ny2 = cfg.lookup("CellsInDirection_y2");

     	try{
       		int Ny3 = cfg.lookup("CellsInDirection_y3");
        	printf("\n\tNy3        :: # Cells in y-direction, third section    ::     \t\t %d \n", Ny3);
      	}
      	catch(const SettingNotFoundException &nfex){
                std::cerr << "No 'CellsInDirection_y3' setting in configuration file." << endl;
      	}
      	int Ny3 = cfg.lookup("CellsInDirection_y3");

      	try{
       		int Nz = cfg.lookup("CellsInDirection_z");
        	printf("\n\tNz         :: # Cells in z-direction                   ::     \t\t %d \n", Nz);
      	}
      	catch(const SettingNotFoundException &nfex){
                std::cerr << "No 'CellsInDirection_z' setting in configuration file." << endl;
      	}
      	int Nz = cfg.lookup("CellsInDirection_z");

      	try{
       		int Nr = cfg.lookup("CellsInRadialDirection");
        	printf("\n\tNr         :: # Cells in radial direction (towards cyl)::     \t\t %d \n", Nr);
      	}
      	catch(const SettingNotFoundException &nfex){
                std::cerr << "No 'CellsInRadialDirection' setting in configuration file." << endl;
      	}
      	int Nr = cfg.lookup("CellsInRadialDirection");

      	try{
       		double radGrad = cfg.lookup("RadialGrading");
        	printf("\n\tradGrad    :: Radial Grading                           ::     \t\t %3.2f \n", radGrad);
      	}
      	catch(const SettingNotFoundException &nfex){
                std::cerr << "No 'RadialGrading' setting in configuration file." << endl;
      	}
      	double radGrad = cfg.lookup("RadialGrading");

	printf("\nParameter successfull loaded... \n");

//	Start blockMeshDict Generation         	
	
	int n_vert = 20;

	double vertices_z0[n_vert][3] = {0.0};
	double vertices_z1[n_vert][3] = {0.0};

	double RC_sin = RC*sin(PI*.25);
	double RR_sin = RR*sin(PI*.25);

	// compute length based on the cylinder diameter 

	double Cx = CCx*RC;
	double Cy = CCy*RC;
	
	H = H*RC;
	W = W*RC;
	D = D*RC;

	// x-segments
	
	double Dx_1 = Cx - RR_sin;
	double Dx_2 = Cx + RR_sin;
	
	// y-segments
	
	double Dy_1 = Cy - RR_sin;
	double Dy_2 = Cy + RR_sin;

	// x-segments cylinder
	
	double dx_1 = Cx - RC_sin;
	double dx_2 = Cx + RC_sin;

	// y-segemnts cylinder 

	double dy_1 = Cy - RC_sin;
	double dy_2 = Cy + RC_sin;

	// radial refinement towards the cylinder surface 
	
   	double radGrad_inv = 1.0/radGrad;

	// boxes x_direction

	vertices_z0[1][0] = Dx_1;
	vertices_z0[2][0] = Dx_2;
	vertices_z0[3][0] = W;

	vertices_z0[4][0] = 0.0;
	vertices_z0[5][0] = Dx_1;
	vertices_z0[6][0] = Dx_2;
	vertices_z0[7][0] = W;

	vertices_z0[8][0] = 0.0;
	vertices_z0[9][0] = Dx_1;
	vertices_z0[10][0] = Dx_2;
	vertices_z0[11][0] = W;

	vertices_z0[12][0] = 0.0;
	vertices_z0[13][0] = Dx_1;
	vertices_z0[14][0] = Dx_2;
	vertices_z0[15][0] = W;

	// box y-direction

	vertices_z0[4][1] = Dy_1;
	vertices_z0[5][1] = Dy_1;
	vertices_z0[6][1] = Dy_1;
	vertices_z0[7][1] = Dy_1;

	vertices_z0[8][1] = Dy_2;
	vertices_z0[9][1] = Dy_2;
	vertices_z0[10][1] = Dy_2;
	vertices_z0[11][1] = Dy_2;

	vertices_z0[12][1] = H;
	vertices_z0[13][1] = H;
	vertices_z0[14][1] = H;
	vertices_z0[15][1] = H;

	// cylinder 
	
	vertices_z0[16][0] = dx_1;
	vertices_z0[17][0] = dx_2;
	vertices_z0[18][0] = dx_1;
	vertices_z0[19][0] = dx_2;

	vertices_z0[16][1] = dy_1;
	vertices_z0[17][1] = dy_1;
	vertices_z0[18][1] = dy_2;
	vertices_z0[19][1] = dy_2;

	// vertices in z-direction

	for (int i=0;i<n_vert;i++){
		vertices_z1[i][0] = vertices_z0[i][0];
		vertices_z1[i][1] = vertices_z0[i][1];
		vertices_z1[i][2] = D;
	}
	
	// write blockMeshDict to file 

	ofstream blockMeshDict;

	blockMeshDict.open("blockMeshDict");

	// set writing precision 
	//
	//	precision == 8 digits
	//

	blockMeshDict << std::fixed << std::setprecision(8);

	blockMeshDict <<  "/*--------------------------------*- C++ -*----------------------------------*|\n";
	blockMeshDict <<  "  ========                 |\n";
    	blockMeshDict <<  "  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n";
     	blockMeshDict <<  "   \\    /   O peration     | Website:  https://openfoam.org\n";
	blockMeshDict <<  "    \\  /    A nd           | Version:  8\n";
	blockMeshDict <<  "     \\/     M anipulation  |\n";
	blockMeshDict <<  "|*---------------------------------------------------------------------------*/\n";
	blockMeshDict <<  "FoamFile\n";
	blockMeshDict <<  "{\n";
	blockMeshDict <<  "     version     2.0;\n";
	blockMeshDict <<  "     format      ascii;\n";
	blockMeshDict <<  "     class       dictionary;\n";
	blockMeshDict <<  "     object      blockMeshDict;\n";
	blockMeshDict <<  "}\n";
	blockMeshDict <<  "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n";

	// add vertices 
	
	blockMeshDict << "vertices\n";
	blockMeshDict << "(\n";

	for (int i=0;i<n_vert;i++){
		blockMeshDict << "\t( " << vertices_z0[i][0] << "  " << vertices_z0[i][1] << "  " << vertices_z0[i][2] << " )  // " << i << "\n";
	}
	for (int i=0;i<n_vert;i++){
		blockMeshDict << "\t( " << vertices_z1[i][0] << "  " << vertices_z1[i][1] << "  " << vertices_z1[i][2] << " )  // " << i+20 << "\n";
	}

	blockMeshDict << ");\n";	

	// add blocks

	blockMeshDict << "blocks\n";
	blockMeshDict << "(\n";

	blockMeshDict << "\t hex (0 1 5 4 20 21 25 24) (" << Nx1 << " " << Ny1 << " " << Nz << ") simpleGrading (1 1 1)\n";
	blockMeshDict << "\t hex (1 2 6 5 21 22 26 25) (" << Nx2 << " " << Ny1 << " " << Nz << ") simpleGrading (1 1 1)\n";
	blockMeshDict << "\t hex (2 3 7 6 22 23 27 26) (" << Nx3 << " " << Ny1 << " " << Nz << ") simpleGrading (1 1 1)\n";
  
	blockMeshDict << "\t hex (4  5  9  8  24 25 29 28) (" << Nx1 << " " << Ny2 << " " << Nz << ") simpleGrading (1 1 1)\n";
	blockMeshDict << "\t hex (6  7  11 10 26 27 31 30) (" << Nx3 << " " << Ny2 << " " << Nz << ") simpleGrading (1 1 1)\n";

	blockMeshDict << "\t hex (8  9  13 12 28 29 33 32) (" << Nx1 << " " << Ny3 << " " << Nz << ") simpleGrading (1 1 1)\n";
	blockMeshDict << "\t hex (9  10 14 13 29 30 34 33) (" << Nx2 << " " << Ny3 << " " << Nz << ") simpleGrading (1 1 1)\n";
	blockMeshDict << "\t hex (10 11 15 14 30 31 35 34) (" << Nx3 << " " << Ny3 << " " << Nz << ") simpleGrading (1 1 1)\n";
  
        // clyinder 
        blockMeshDict << "\t hex (5  6  17 16 25 26 37 36) (" << Nx2 << " " << Nr << " " << Nz << ") simpleGrading (1 " << radGrad << " 1)\n";
	blockMeshDict << "\t hex (17 6  10 19 37 26 30 39) (" << Nr << " " << Ny2 << " " << Nz << ") simpleGrading (" << radGrad_inv << " 1 1)\n";
	blockMeshDict << "\t hex (18 19 10 9  38 39 30 29) (" << Nx2 << " " << Nr << " " << Nz << ") simpleGrading (1 " << radGrad_inv <<" 1)\n";
	blockMeshDict << "\t hex (5  16 18 9  25 36 38 29) (" << Nr << " " << Ny2 << " " << Nz << ") simpleGrading (" << radGrad << " 1 1)\n";

	blockMeshDict << ");\n";	

	// add edges
	
	blockMeshDict << "edges\n";
	blockMeshDict << "(\n";       

	blockMeshDict << "\t arc 5  6  (" << Cx << " " << Cy-RR << " 0)\n";
	blockMeshDict << "\t arc 25 26 (" << Cx << " " << Cy-RR << " " << D << ")\n";
  
        blockMeshDict << "\t arc 5  9  (" << Cx-RR << " " << Cy << "  0)\n";      
	blockMeshDict << "\t arc 25 29 (" << Cx-RR << " " << Cy << " " << D << ")\n";

	blockMeshDict << "\t arc 6  10 (" << Cx+RR << " " << Cy << " 0)\n";      
	blockMeshDict << "\t arc 26 30 (" << Cx+RR << " " << Cy << " " << D << ")\n";
  
	blockMeshDict << "\t arc 9  10 (" << Cx << " " << Cy+RR << " 0)\n";
	blockMeshDict << "\t arc 29 30 (" << Cx << " " << Cy+RR << " " << D << ")\n";
  
          // cylinder wall
	blockMeshDict << "\t arc 16 17 (" << Cx << " " << Cy-RC << " 0)\n";
	blockMeshDict << "\t arc 36 37 (" << Cx << " " << Cy-RC << " " << D << ")\n";
  
	blockMeshDict << "\t arc 17 19 (" << Cx+RC << " " << Cy << " 0)\n";
	blockMeshDict << "\t arc 37 39 (" << Cx+RC << " " << Cy << " " << D << ")\n";;

        blockMeshDict << "\t arc 18 19 (" << Cx << " " << Cy+RC << " 0)\n";
 	blockMeshDict << "\t arc 38 39 (" << Cx << " " << Cy+RC << " " << D << ")\n";
  
	blockMeshDict << "\t arc 16 18 (" << Cx-RC << " " << Cy << " 0)\n";
 	blockMeshDict << "\t arc 36 38 (" << Cx-RC << " " << Cy << " " << D << ")\n";
  
 	blockMeshDict  << ");\n";
	
	// define boundary patches || normal vectors have to point outside of the domain!!!
	
	blockMeshDict << "boundary\n";
	blockMeshDict << "(\n";

	// inlet faces
	blockMeshDict << "\tinlet\n";
	blockMeshDict << "\t{\n";

	blockMeshDict << "\t\ttype \t patch;\n";
	blockMeshDict << "\t\tfaces\n";
	blockMeshDict << "\t\t(\n";

	blockMeshDict << "\t\t(0 4 24 20)\n";
	blockMeshDict << "\t\t(4 8 24 28)\n";
	blockMeshDict << "\t\t(8 12 32 28)\n";

	blockMeshDict << "\t\t);\n";
	blockMeshDict << "\t}\n";

	// outlet faces 
	
	blockMeshDict << "\toutlet\n";
	blockMeshDict << "\t{\n";

	blockMeshDict << "\t\ttype patch;\n";
	blockMeshDict << "\t\tfaces\n";
	blockMeshDict << "\t\t(\n";

	blockMeshDict << "\t\t(3 23 27 7)\n";
	blockMeshDict << "\t\t(7 27 31 11)\n";
	blockMeshDict << "\t\t(11 31 35 15)\n";

	blockMeshDict << "\t\t);\n";
	blockMeshDict << "\t}\n";

	// cylinder wall faces

	blockMeshDict << "\tcylinder\n";
	blockMeshDict << "\t{\n";

	blockMeshDict << "\t\ttype wall;\n";
	blockMeshDict << "\t\tfaces\n";
	blockMeshDict << "\t\t(\n";

	blockMeshDict << "\t\t(17 19 39 37)\n";
	blockMeshDict << "\t\t(18 38 39 19)\n";
	blockMeshDict << "\t\t(16 36 38 18)\n";
	blockMeshDict << "\t\t(16 17 36 37)\n";

	blockMeshDict << "\t\t);\n";
	blockMeshDict << "\t}\n";

	// top and bottom symmetryPlane
	

//	blockMeshDict << "\ttop\n";
//	blockMeshDict << "\t{\n";

//	blockMeshDict << "\t\ttype symmetryPlane;\n";
//	blockMeshDict << "\t\tfaces\n";
//	blockMeshDict << "\t\t(\n";

	blockMeshDict << "\ttop0\n";
	blockMeshDict << "\t{\n";

	blockMeshDict << "\t\ttype cyclic;\n";
    blockMeshDict << "\t\tneighbourPatch bottom0;\n";
	blockMeshDict << "\t\tfaces\n";
	blockMeshDict << "\t\t(\n";

    blockMeshDict << "\t\t(12 32 33 13)\n";

	blockMeshDict << "\t\t);\n";
	blockMeshDict << "\t}\n";

	blockMeshDict << "\ttop1\n";
	blockMeshDict << "\t{\n";

	blockMeshDict << "\t\ttype cyclic;\n";
    blockMeshDict << "\t\tneighbourPatch bottom1;\n";
	blockMeshDict << "\t\tfaces\n";
	blockMeshDict << "\t\t(\n";

	blockMeshDict << "\t\t(13 33 34 14)\n";

	blockMeshDict << "\t\t);\n";
	blockMeshDict << "\t}\n";

	blockMeshDict << "\ttop2\n";
	blockMeshDict << "\t{\n";

	blockMeshDict << "\t\ttype cyclic;\n";
    blockMeshDict << "\t\tneighbourPatch bottom2;\n";
	blockMeshDict << "\t\tfaces\n";
	blockMeshDict << "\t\t(\n";

	blockMeshDict << "\t\t(14 34 35 15)\n";

	blockMeshDict << "\t\t);\n";
	blockMeshDict << "\t}\n";

//    blockMeshDict << "\t\t(12 32 33 13)\n";
//	blockMeshDict << "\t\t(13 33 34 14)\n";
//	blockMeshDict << "\t\t(14 34 35 15)\n";

//	blockMeshDict << "\t\t);\n";
//	blockMeshDict << "\t}\n";

//	blockMeshDict << "\tbottom\n";
//	blockMeshDict << "\t{\n";

//	blockMeshDict << "\t\ttype symmetryPlane;\n";
//	blockMeshDict << "\t\tfaces\n";
//	blockMeshDict << "\t\t(\n";
	
//	blockMeshDict << "\t\t(0 1 21 20)\n";
//	blockMeshDict << "\t\t(1 2 22 21)\n";
//	blockMeshDict << "\t\t(2 3 23 22)\n";

//	blockMeshDict << "\t\t);\n";
//	blockMeshDict << "\t}\n";

	blockMeshDict << "\tbottom0\n";
	blockMeshDict << "\t{\n";

	blockMeshDict << "\t\ttype cyclic;\n";
    blockMeshDict << "\t\tneighbourPatch top0;\n";
	blockMeshDict << "\t\tfaces\n";
	blockMeshDict << "\t\t(\n";
	
	blockMeshDict << "\t\t(0 1 21 20)\n";

	blockMeshDict << "\t\t);\n";
	blockMeshDict << "\t}\n";

	blockMeshDict << "\tbottom1\n";
	blockMeshDict << "\t{\n";

	blockMeshDict << "\t\ttype cyclic;\n";
    blockMeshDict << "\t\tneighbourPatch top1;\n";
	blockMeshDict << "\t\tfaces\n";
	blockMeshDict << "\t\t(\n";
	
	blockMeshDict << "\t\t(1 2 22 21)\n";

	blockMeshDict << "\t\t);\n";
	blockMeshDict << "\t}\n";

	blockMeshDict << "\tbottom2\n";
	blockMeshDict << "\t{\n";

	blockMeshDict << "\t\ttype cyclic;\n";
    blockMeshDict << "\t\tneighbourPatch top2;\n";
	blockMeshDict << "\t\tfaces\n";
	blockMeshDict << "\t\t(\n";
	
	blockMeshDict << "\t\t(2 3 23 22)\n";

	blockMeshDict << "\t\t);\n";
	blockMeshDict << "\t}\n";

	// front an back faces 
	//
	// 	if Nz > 1 the boundary condition symmetryPlane is applied, otherwise 
	// 	the faces are treated as empty.
	//

	if (Nz > 1) // symmetryPlane
	{

	blockMeshDict << "\tfront\n";
	blockMeshDict << "\t{\n";

	blockMeshDict << "\t\ttype symmetryPlane;\n";
	blockMeshDict << "\t\tfaces\n";
	blockMeshDict << "\t\t(\n";

	blockMeshDict << "\t\t(0 4 5 1)\n";
	blockMeshDict << "\t\t(1 5 6 2)\n";
	blockMeshDict << "\t\t(2 6 7 3)\n";
	blockMeshDict << "\t\t(4 8 9 5)\n";
	blockMeshDict << "\t\t(5 9 18 16)\n";
	blockMeshDict << "\t\t(5 16 17 6)\n";
	blockMeshDict << "\t\t(17 19 10 6)\n";
	blockMeshDict << "\t\t(18 9 10 19)\n";
	blockMeshDict << "\t\t(6 10 11 7)\n";
	blockMeshDict << "\t\t(8 12 13 9)\n";
	blockMeshDict << "\t\t(9 13 14 10)\n";
	blockMeshDict << "\t\t(10 14 15 11)\n";

	blockMeshDict << "\t\t);\n";
	blockMeshDict << "\t}\n";

	blockMeshDict << "\tback\n";
	blockMeshDict << "\t{\n";

	blockMeshDict << "\t\ttype symmetryPlane;\n";
	blockMeshDict << "\t\tfaces\n";
	blockMeshDict << "\t\t(\n";
	
	blockMeshDict << "\t\t(20 21 25 24)\n";
	blockMeshDict << "\t\t(21 22 26 25)\n";
	blockMeshDict << "\t\t(22 23 27 26)\n";
	blockMeshDict << "\t\t(24 25 29 28)\n";
	blockMeshDict << "\t\t(25 36 38 29)\n";
	blockMeshDict << "\t\t(25 26 37 36)\n";
	blockMeshDict << "\t\t(37 26 30 39)\n";
	blockMeshDict << "\t\t(38 39 30 29)\n";
	blockMeshDict << "\t\t(26 27 31 30)\n";
	blockMeshDict << "\t\t(28 29 32 33)\n";
	blockMeshDict << "\t\t(29 30 34 33)\n";
	blockMeshDict << "\t\t(30 31 35 34)\n";

	blockMeshDict << "\t\t);\n";
	blockMeshDict << "\t}\n";

	}
	else // empty faces (no solving in z-direction)
	{

	blockMeshDict << "\tFrontAndBack\n";
	blockMeshDict << "\t{\n";

	blockMeshDict << "\t\ttype empty;\n";
	blockMeshDict << "\t\tfaces\n";
	blockMeshDict << "\t\t(\n";

	blockMeshDict << "\t\t(0 4 5 1)\n";
	blockMeshDict << "\t\t(1 5 6 2)\n";
	blockMeshDict << "\t\t(2 6 7 3)\n";
	blockMeshDict << "\t\t(4 8 9 5)\n";
	blockMeshDict << "\t\t(5 9 18 16)\n";
	blockMeshDict << "\t\t(5 16 17 6)\n";
	blockMeshDict << "\t\t(17 19 10 6)\n";
	blockMeshDict << "\t\t(18 9 10 19)\n";
	blockMeshDict << "\t\t(6 10 11 7)\n";
	blockMeshDict << "\t\t(8 12 13 9)\n";
	blockMeshDict << "\t\t(9 13 14 10)\n";
	blockMeshDict << "\t\t(10 14 15 11)\n";

	blockMeshDict << "\t\t(20 21 25 24)\n";
	blockMeshDict << "\t\t(21 22 26 25)\n";
	blockMeshDict << "\t\t(22 23 27 26)\n";
	blockMeshDict << "\t\t(24 25 29 28)\n";
	blockMeshDict << "\t\t(25 36 38 29)\n";
	blockMeshDict << "\t\t(25 26 37 36)\n";
	blockMeshDict << "\t\t(37 26 30 39)\n";
	blockMeshDict << "\t\t(38 39 30 29)\n";
	blockMeshDict << "\t\t(26 27 31 30)\n";
	blockMeshDict << "\t\t(28 29 32 33)\n";
	blockMeshDict << "\t\t(29 30 34 33)\n";
	blockMeshDict << "\t\t(30 31 35 34)\n";

	blockMeshDict << "\t\t);\n";
	blockMeshDict << "\t}\n";

	}


	blockMeshDict << ");\n"; // boundary bracket 
	

	// occasional mergePatchPairs statement 
	
	blockMeshDict << "mergePatchPairs\n";
	blockMeshDict << "();\n";


	blockMeshDict <<  "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //";

	blockMeshDict.close();


      	return(EXIT_SUCCESS);
}






