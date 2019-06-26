#ifndef GET_FROM_GAUSSIAN_OUTPUT_H
#define GET_FROM_GAUSSIAN_OUTPUT_H
#include <vector>
#include <string>
#include <Eigen/Dense>

int print_vector_int(std::vector<int> A);

int print_vector_double(std::vector<double> A);

int print_vector_string(std::vector<std::string> A);

int print_matrix_int(std::vector<std::vector<int> > A);

bool fexists(const char *filename);

int if_file_exist_delete (std::string filename);

int x_rotation_matrix (double cosx, double sinx, Eigen::Matrix3d& Rx);

int y_rotation_matrix (double cosy, double siny, Eigen::Matrix3d& Ry);

int align (Eigen::Vector3d vec_in, Eigen::MatrixXd& coords);

int read_input_file(int& nmols, std::vector<std::string>& filenames, \
      std::vector<int>& nrings_per_molecule, \
      std::vector<std::vector<int>>& rings_atoms);

#endif

