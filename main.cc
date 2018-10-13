#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <random>

#include "./dm_decomp.h"

int main(int argc,char **argv){
  int rows = 5 , cols = 5;
  if(argc<2){
    std::cout << argv[0] << " [matrix filename]" << std::endl;
    return 0;
  }
  std::ifstream ifs(argv[1]);
  ifs >> rows >> cols;
  Eigen::SparseMatrix<double> A(rows,cols);
  std::vector<Eigen::Triplet<double> > t;

  for(int i=0;i<rows;++i){
    for(int j=0;j<cols;++j){
      double a;
      ifs >> a;
      if(a!=0){
        t.emplace_back(i,j,a);
      }
    }
  }

  A.setFromTriplets(t.begin(),t.end());

  Eigen::MatrixXd dmat = Eigen::MatrixXd(A);
  std::cout << "original matrix:" << std::endl;
  for(int i=0;i<rows;++i){
    for(int j=0;j<cols;++j){
      std::cout << std::setw(3) << dmat(i,j) ;
    }
    std::cout << std::endl;
  }

  dm_decomposition dm;
  dm(A);
  return 0;
}
