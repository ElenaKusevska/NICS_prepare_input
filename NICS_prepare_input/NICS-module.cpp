#include <array>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <math.h>
#include <Eigen/Dense>

int print_vector_int(std::vector<int> A) {
   int i;
   for (i=0; i<A.size(); i=i+1) {
      std::cout << A[i] <<  " " ;
   }
   std::cout << " " << std::endl;
   std::cout << " " << std::endl;
   return 0;
}

int print_vector_double(std::vector<double> A) {
   int i;
   for (i=0; i<A.size(); i=i+1) {
      std::cout << A[i] <<  " " ;
   }
   std::cout << " " << std::endl;
   std::cout << " " << std::endl;
   return 0;
}

int print_vector_string(std::vector<std::string> A) {
   int i;
   for (i=0; i<A.size(); i=i+1) {
      std::cout << A[i] << " ";
   }
   std::cout << " " << std::endl;
   std::cout << " " << std::endl;
   return 0;
}

int print_matrix_int(std::vector<std::vector<int> > A) {
   int i, j;
   for (i=0; i<A.size(); i=i+1) {
      for (j=0; j<A[i].size(); j=j+1) {
         std::cout << A[i][j] << " ";
      }
      std::cout << " " << std::endl;
   }
   std::cout << " " << std::endl;
   return 0;
}

bool fexists(const char *filename) {
   std::ifstream ifile(filename);
   return ifile.good();
}

int if_file_exist_delete (std::string filename) {
   if (fexists(filename.c_str())) {
      if (std::remove(filename.c_str()) != 0) {
          std::cout << "failed to remove " << filename << std::endl;
          exit(1);
       }
       else {
          std::cout << filename << " found and deleted " << std::endl;
       }
   }
   return 0;
}

int x_rotation_matrix (double cosx, double sinx, Eigen::Matrix3d& Rx) {

   Rx(0,0)=1;
   Rx(0,1)=0;
   Rx(0,2)=0;
   Rx(1,0)=0;
   Rx(1,1)=cosx;
   Rx(1,2)=-sinx;
   Rx(2,0)=0;
   Rx(2,1)=sinx;
   Rx(2,2)=cosx;

   std::cout << "Rx " << Rx << std::endl;

   return 0;
}

int y_rotation_matrix (double cosy, double siny, Eigen::Matrix3d& Ry) {

   Ry(0,0)=cosy;
   Ry(0,1)=0;
   Ry(0,2)=-siny;
   Ry(1,0)=0;
   Ry(1,1)=1;
   Ry(1,2)=0;
   Ry(2,0)=siny;
   Ry(2,1)=0;
   Ry(2,2)=cosy;

   std::cout << "Ry " << Ry << std::endl;

   return 0;
}

int align (Eigen::Vector3d vec_in, Eigen::MatrixXd& coords) {

   double cosx, sinx, cosy, siny;
   Eigen::Matrix3d Rx, Ry;
   Eigen::Vector3d vecz;
   Eigen::MatrixXd coords_t(coords.cols(),coords.rows());
   coords_t = coords.transpose();

   if (vec_in(0,0) == 0) {
      if (vec_in(1,0) != 0) { // if vec_in = [0,y,z]
         
         // find matrix for rotation about the x-axis:
         cosx = vec_in(2,0)/sqrt( pow(vec_in(1,0),2) + pow(vec_in(2,0),2) );
         sinx = vec_in(1,0)/sqrt( pow(vec_in(1,0),2) + pow(vec_in(2,0),2) );
         x_rotation_matrix(cosx, sinx, Rx);

         // rotate space about the x - axis to bring vec_in
         // in the xz-plane:
         coords_t = Rx*coords_t;
         vecz = Rx*vec_in;
      }
      else if (vec_in(1,0) == 0) {// [0,0,z]
         // all is good
         vecz = vec_in;
      }
   }
   if (vec_in(1,0) == 0) {
      if (vec_in(0,0) != 0) { // [x, 0, z]
         
         // find matrix for rotation about the y-axis:
         cosy = vec_in(2,0)/sqrt( pow(vec_in(0,0),2) + pow(vec_in(2,0),2) );
         siny = vec_in(0,0)/sqrt( pow(vec_in(0,0),2) + pow(vec_in(2,0),2) );
         y_rotation_matrix (cosy, siny, Ry);

         // rotate space around the y - axis, so that the 
         // rotation axis lies along the positive z - axis:
         coords_t = Ry*coords_t;
         vecz = Ry*vec_in; // should have the form [0,0,z]
      }
   }
   if (vec_in(0,0) != 0 ) {
      if (vec_in(1,0) != 0) { // [x, y, z]

         Eigen::Vector3d vecxz;

         // find matrix for rotation about the x-axis:
         cosx = vec_in(2,0)/sqrt( pow(vec_in(1,0),2) + pow(vec_in(2,0),2) );
         sinx = vec_in(1,0)/sqrt( pow(vec_in(1,0),2) + pow(vec_in(2,0),2) );
         x_rotation_matrix (cosx, sinx, Rx);
        
         // rotate space about the x - axis to bring vec_in
         // in the xz-plane:
         std::cout << Ry.rows() << Ry.cols() << coords.rows() << coords.cols()    ;
         coords_t = Rx*coords_t;
         vecxz = Rx*vec_in;
        
         //find matrix for rotation about the y-axis
         cosy = vecxz(2,0)/sqrt( pow(vecxz(0,0),2) + pow(vecxz(2,0),2) );
         siny = vecxz(0,0)/sqrt( pow(vecxz(0,0),2) + pow(vecxz(2,0),2) );
         y_rotation_matrix (cosy, siny, Ry);

         // rotate space around the y - axis, so that the rotation axis lies
         // along the positive z - axis:
         std::cout << Ry.rows() << Ry.cols() << coords.rows() << coords.cols();
         coords_t = Ry*coords_t;
         vecz = Ry*vecxz; //should have the form [0,0,z]
      }
   }
   std::cout << "vecz: " << vecz << std::endl;
   return 0;
}

int read_input_file(int& nmols, std::vector<std::string>& filenames, \
      std::vector<int>& nrings_per_molecule, \
      std::vector<std::vector<int>>& rings_atoms) {

   std::string line, temp;
   char test_char;
   int nrings_total, nfile, counts, finds_nl, rr, temp_int;

   //----------------------------------
   // count number of molecules,
   // and the total number of rings:
   //----------------------------------

   std::ifstream inputfile;
   inputfile.open("input");
   nmols = 0;
   finds_nl = 0; // finds new line
   nrings_total = 0;
   while (!inputfile.eof()) {
      inputfile.get(test_char);
      if (finds_nl == 0) { // if this is the first character in the line
         if (test_char == 'm') { // and the first character is 'm'
            nmols = nmols + 1;
         }
         if (test_char == 'r') {
            nrings_total = nrings_total + 1;
         }
      }
      finds_nl = finds_nl + 1; // because it's probably not the
                                       // end of the line
      if (test_char == '\n') { // though, if it is the end of the line:
         finds_nl = 0;
      }
   }
   std::cout << "nrings total: " << nrings_total << std::endl;

   //---------------------------------------------
   // get the filenames of the molecules,
   // count the total numbers of rings, and the
   // number of rings per molecule:
   //---------------------------------------------

   filenames.resize(nmols);
   nrings_per_molecule.resize(nmols);
   rings_atoms.resize(nrings_total);
   nfile = -1; // because c++ array numbering starts from 0
   rr = 0; // the current ring being read
   finds_nl = 0; // finds new line
   inputfile.clear(); // To clear the EOF from previously
   inputfile.seekg(0, std::ios::beg); // set to beggining of file
   while (!inputfile.eof()) {
      inputfile.get(test_char);
      if (finds_nl == 0) { // if this is the first character in the line
         if (test_char == 'm') { // we have a new molecule
            nfile = nfile + 1;
            nrings_per_molecule[nfile] = 0; // initialize to 0 rings
         }
         if (test_char == 'f') { // the filename of the new molecule is
                                 // written here
            getline(inputfile,line); // get this line
            std::stringstream ssin(line); // break up line in string stream
            counts = 1;
            while (ssin.good()){
               ssin >> temp;
               if (counts = 2) {
                  filenames[nfile] = temp;
               }
               counts = counts + 1;
            }
            // go back by 1 character so that get reads '\n' from
            // this line which was read by getline:
            inputfile.seekg(-1, std::ios::cur);
         }
         if (test_char == 'r') { 
            // there is another ring for this molecule:
            nrings_per_molecule[nfile] = nrings_per_molecule[nfile] + 1;
            // get the atoms that make up this ring
            getline(inputfile,line); // get this line
            std::stringstream ssin(line); // break up line in string stream
            counts = 1;
            while (ssin.good()){
               ssin >> temp;
               std::cout << temp << " ";
               if (counts > 1) {
                  temp_int = std::atoi(temp.c_str());
                  rings_atoms[rr].push_back(temp_int);
               }
               counts = counts + 1;
            }
            std::cout << " " << std::endl;
            // move on to next ring
            rr = rr + 1;
            // go back by 1 character so that get reads '\n' from
            // this line which was read by getline:
            inputfile.seekg(-1, std::ios::cur);
         }
      }
      finds_nl = finds_nl + 1; // because it's probably not the
                                       // end of the line
      if (test_char == '\n') { // though, if it is the end of the line:
         finds_nl = 0;
      }
   }

   inputfile.close();
   return 0;
}
