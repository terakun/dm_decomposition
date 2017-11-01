#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <random>

#include "./dm_decomp.h"

int main(){
  int rows = 5 , cols = 5;
  Eigen::SparseMatrix<double> A(rows,cols);
  int elemsize = 10; // 成分の数
  std::vector<Eigen::Triplet<double> > t(elemsize);

  std::random_device rnd;
  std::mt19937 mt(rnd());

  std::uniform_int_distribution<> row_dist(0,rows-1);
  std::uniform_int_distribution<> col_dist(0,cols-1);
  std::uniform_real_distribution<> dist(1,2);
  std::uniform_int_distribution<> int_dist(1,5);
  for(int i=0;i<elemsize;++i){
    t[i] = Eigen::Triplet<double>(row_dist(mt),col_dist(mt),int_dist(mt));
  }

  A.setFromTriplets(t.begin(),t.end());

  for(int i=0;i<A.outerSize();++i){
    for(Eigen::SparseMatrix<double>
        ::InnerIterator it(A,i);it;++it){
      std::cout << it.value() <<",";
      std::cout << it.row() <<",";//SVector にはない
      std::cout << it.col() <<",";//SVector にはない
      std::cout << it.index() << std::endl;
    }
  }

  Eigen::MatrixXd dmat = Eigen::MatrixXd(A);
  std::cout << dmat << std::endl;

  dm_decomposition dm;
  dm(A);
  return 0;
}
