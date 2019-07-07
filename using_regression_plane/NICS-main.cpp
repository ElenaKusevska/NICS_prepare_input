#include "NICS-module.h"
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <math.h>
#include <Eigen/Dense>

//-------------------------------------------------------------
// Program to find the centers of rings fo NICS(0) and NICS(1) 
// calculation and write some quazi-.gjf file
//-------------------------------------------------------------

int main() {

   //---------------------------------------------
   // Get the information out of the input file:
   //--------------------------------------------- 

   // Variables to get out of input file:
   int nmols;
   std::vector<int> nrings_per_molecule;
   std::vector<std::string> filenames;
   std::vector<std::vector<int>> rings_atoms; // atoms making up rings
   
   // Read input file:
   read_input_file(nmols, filenames, nrings_per_molecule, rings_atoms);
   std::cout << "nmols: " << nmols << std::endl;
   print_vector_string(filenames);
   std::cout << "nrings_per_molecule: " ;
   print_vector_int(nrings_per_molecule);
   print_matrix_int(rings_atoms);

   //--------------------------------------------------------
   // Get the coordinates of the molecules and rings from 
   // the .xyz files, move and reorient them, and print the 
   // final coordinates:
   //--------------------------------------------------------

   // Variables:
   std::string in_xyz_file_name, out_ring_xyz_file_name, out_xyz_file_name;
   std::string line, temp;
   std::vector<std::string> lables_ring, lables;
   std::vector<double> mass, xring, yring, zring,  x, y, z;
   double temp_double;
   double xsum, ysum, zsum, xcenter, ycenter, zcenter;
   int i, j, k, kk, line_num, ring_num;
   Eigen::Vector3d vec1, vec2, normal1, normal2;
   Eigen::Matrix3d inert;
   Eigen::MatrixXd ring_coordinates(1,1), coordinates(1,1);

   ring_num = 0;
   for (i=0; i<nmols; i=i+1) { // for each of the molecules:
      for (j=0; j<nrings_per_molecule[i]; j=j+1) { // for each ring in mol
         
         // Reset lables and coordinates arrays:
         mass.resize(0);
         lables.resize(0);
         x.resize(0);
         y.resize(0);
         z.resize(0);
         lables_ring.resize(0);
         xring.resize(0);
         yring.resize(0);
         zring.resize(0);

         // Open the input file for that molecule and ring:
         in_xyz_file_name = filenames[i] + ".xyz";
         std::ifstream in_xyz_file;
         in_xyz_file.open(in_xyz_file_name);

         // get the coordinates of the rings and the molecule:
         line_num = 1;
         while (!in_xyz_file.eof()) {
            getline(in_xyz_file,line);
            if (line_num > 2) {
               std::stringstream ssin(line);
               std::vector<std::string> split_line;
               split_line.resize(0);
               while (ssin.good()){
                  ssin >> temp;
                  split_line.push_back(temp);
               }
               // uncomment if segmentation fault stod
               //std::cout << print_vector_string(split_line) << std::endl;
               if (split_line.size() > 1) {
                  // get lables and coordinates for whole molecule
                  lables.push_back(split_line[0]);
                  temp_double = std::stod(split_line[1]);
                  x.push_back(temp_double);
                  temp_double = std::stod(split_line[2]);
                  y.push_back(temp_double);
                  temp_double = std::stod(split_line[3]);
                  z.push_back(temp_double);
                  for (k=0; k<rings_atoms[ring_num].size(); k=k+1){
                     if (line_num == rings_atoms[ring_num][k] + 2) {
                        // get lables and coordinates for ring:
                        lables_ring.push_back(split_line[0]);
                        temp_double = std::stod(split_line[1]);
                        xring.push_back(temp_double);
                        temp_double = std::stod(split_line[2]);
                        yring.push_back(temp_double);
                        temp_double = std::stod(split_line[3]);
                        zring.push_back(temp_double);
                        continue;
                     }
                  }
               }
            }
            line_num = line_num + 1;
         }

         // Find center of ring:
         xsum = 0.0;
         ysum = 0.0;
         zsum = 0.0;
         for (k=0; k<xring.size(); k=k+1) {
            xsum = xsum + xring[k];
            ysum = ysum + yring[k];
            zsum = zsum + zring[k];
         }
         xcenter = xsum/xring.size();
         ycenter = ysum/yring.size();
         zcenter = zsum/zring.size();

         // Move the coordinates vectors for the ring
         // so that the center of the ring is positioned
         // at (0,0,0):
         for (k=0; k<xring.size(); k=k+1) {
            xring[k] = xring[k] - xcenter;
            yring[k] = yring[k] - ycenter;
            zring[k] = zring[k] - zcenter;
         }

         // Move the coordinates vectors for the 
         // entire molecule so that the center of the 
         // ring is positioned at (0,0,0):
         for (k=0; k<x.size(); k=k+1) {
            x[k] = x[k] - xcenter;
            y[k] = y[k] - ycenter;
            z[k] = z[k] - zcenter;
         }

         // find plane passing through all the points of
         // the ring using the equation 
         // C=(XY.transpose*XY).inverse*XY.transpose*Z
         Eigen::Vector3d Cvec;
         Eigen::VectorXd Zvec;
         Eigen::MatrixXd XYmat;

         Zvec.resize(zring.size());
         XYmat.resize(xring.size(),3);
         for (k=0; k<xring.size(); k=k+1) {
            XYmat(k,0) = 1;
            XYmat(k,1) = xring[k];
            XYmat(k,2) = yring[k];
            Zvec(k) = zring[k];
         }

         Cvec = (XYmat.transpose()*XYmat).inverse()*XYmat.transpose()*Zvec;

         // Define two vectors in this plane:
         vec1(0) = 0.8;
         vec1(1) = 1.5;
         vec1(2) = Cvec(0) + Cvec(1)*vec1(0) + Cvec(2)*vec1(1);
         vec2(0) = -1.2;
         vec2(1) = 0.6;
         vec2(2) = Cvec(0) + Cvec(1)*vec2(0) + Cvec(2)*vec2(1);


         // find a normal to this plane at (0,0,0)
         normal1 = vec1.cross(vec2);
         std::cout << "normal:" << std::endl;
         std::cout << normal1 << std::endl;

         double norm_of_normal;
         norm_of_normal = (sqrt(normal1(0)*normal1(0) +
                  normal1(1)*normal1(1) + normal1(2)*normal1(2)));
         std::cout << "norm of normal: " << normal1.norm() << std::endl;

         // normalize normal 1 to length 1:
         for (k=0; k<normal1.size(); k=k+1) {
            normal1(k) = normal1(k) / norm_of_normal;
         }

         norm_of_normal = (sqrt(normal1(0)*normal1(0) +
                  normal1(1)*normal1(1) + normal1(2)*normal1(2)));
         std::cout << "norm of normal: " << norm_of_normal << std::endl;
     
         // convert vectors xring, yring, and zring to a
         // ring_coordinates matrix:
         ring_coordinates.resize(xring.size(),3);
         for(k=0; k<xring.size(); k=k+1) {
            ring_coordinates(k,0) = xring[k];
            ring_coordinates(k,1) = yring[k];
            ring_coordinates(k,2) = zring[k];
         }
 
         // Convert vector x, y, and z to a coordinates matrix
         // and add dummy atoms:
         coordinates.resize(x.size() + 3,3);
         for (k=0; k<x.size(); k=k+1) {
            coordinates(k,0) = x[k];
            coordinates(k,1) = y[k];
            coordinates(k,2) = z[k];
         }
         coordinates(x.size(),0) = normal1(0);
         coordinates(x.size(),1) = normal1(1);
         coordinates(x.size(),2) = normal1(2);
         coordinates(x.size()+1,0) = 0.0;
         coordinates(x.size()+1,1) = 0.0;
         coordinates(x.size()+1,2) = 0.0;
         coordinates(x.size()+2,0) = -normal1(0);
         coordinates(x.size()+2,1) = -normal1(1);
         coordinates(x.size()+2,2) = -normal1(2);

      
         // Print final coordinates for the ring
         out_ring_xyz_file_name = filenames[i] + "_ring_"
            + std::to_string(j) + ".xyz";
         if_file_exist_delete(out_ring_xyz_file_name);
         std::ofstream out_ring_xyz_file;
         out_ring_xyz_file.open(out_ring_xyz_file_name);

         out_ring_xyz_file << rings_atoms[ring_num].size() + 3 << std::endl;
         out_ring_xyz_file << "sldgjksjg" << std::endl;
         for (kk=0; kk<rings_atoms[ring_num].size(); kk=kk+1) {
            out_ring_xyz_file << lables_ring[kk] << "         "
                         << ring_coordinates(kk,0) << "         "
                         << ring_coordinates(kk,1) << "         "
                         << ring_coordinates(kk,2) << std::endl;
         }
         out_ring_xyz_file << "N          " << normal1(0) << "         "
                           << normal1(1) << "         " << normal1(2)
                           << std::endl;
         out_ring_xyz_file << "O          " << vec1(0) << "         "
                           << vec1(1) << "         " << vec1(2)
                           << std::endl;
         out_ring_xyz_file << "O          " << vec2(0) << "         "
                           << vec2(1) << "         " << vec2(2)
                           << std::endl;
         
         // Align normal2 to z-axis, so that
         // the ring is in the xy plane"
         align(normal1, coordinates);
         align(normal1, ring_coordinates);

         // Print final coordinates for the molecule,
         // adding the dummy atoms:
         out_xyz_file_name = filenames[i] + "_" + 
                                  std::to_string(j) + ".xyz";
         if_file_exist_delete(out_xyz_file_name);
         std::ofstream out_xyz_file;
         out_xyz_file.open(out_xyz_file_name);

         // printing:
         out_xyz_file << x.size() + 3 << std::endl;
         out_xyz_file << "0 1" << std::endl;
         for (kk=0; kk<x.size(); kk=kk+1) {
            out_xyz_file << std::left << std::setw(2) << lables[kk] 
                         << std::setw(18) << std::right
                         << std::fixed << std::setprecision(7)
                         << coordinates(kk,0) 
                         << std::setw(19) << std::right
                         << std::fixed << std::setprecision(7)
                         << coordinates(kk,1) 
                         << std::setw(19) << std::right
                         << std::fixed << std::setprecision(7)
                         << coordinates(kk,2) << std::endl;
         }
         for (kk=x.size(); kk<x.size()+3; kk=kk+1) {
         out_xyz_file << std::left << std::setw(2) << "Bq"
                      << std::setw(18) << std::right
                      << std::fixed << std::setprecision(7)
                      << coordinates(kk,0) 
                      << std::setw(19) << std::right
                      << std::fixed << std::setprecision(7)
                      << coordinates(kk,1)
                      << std::setw(19) << std::right
                      << std::fixed << std::setprecision(7)
                      << coordinates(kk,2) << std::endl;
         }

         // End of loop, reset for next ring:
         ring_num = ring_num + 1;
         in_xyz_file.close();
         out_ring_xyz_file.close();
         out_xyz_file.close();
      }

   }

   return 0;
}
