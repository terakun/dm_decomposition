#include <iostream>
#include <fstream>
#include <algorithm>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <random>

#include "./dm_decomp.h"

int main(){
  int rows = 5 , cols = 5;
  std::ifstream ifs("./data.txt");
  ifs >> rows >> cols;
  Eigen::SparseMatrix<double> A(rows,cols);
  std::vector<Eigen::Triplet<double> > t;

  std::random_device rnd;
  std::mt19937 mt(rnd());
  std::vector<int> shuffled_rows(rows),shuffled_cols(cols);
  for(int r=0;r<rows;++r) shuffled_rows[r]=r;
  std::shuffle(shuffled_rows.begin(),shuffled_rows.end(),mt);
  for(int c=0;c<cols;++c) shuffled_cols[c]=c;
  std::shuffle(shuffled_cols.begin(),shuffled_cols.end(),mt);

  for(int i=0;i<rows;++i){
    for(int j=0;j<cols;++j){
      double a;
      ifs >> a;
      if(a!=0){
        t.emplace_back(shuffled_rows[i],shuffled_cols[j],a);
      }
    }
  }

  A.setFromTriplets(t.begin(),t.end());

  Eigen::MatrixXd dmat = Eigen::MatrixXd(A);
  std::cout << "original matrix:" << std::endl;
  std::cout << dmat << std::endl;

  dm_decomposition dm;
  dm(A);
  return 0;
}
