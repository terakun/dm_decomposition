#ifndef DM_DECOMP_H
#define DM_DECOMP_H
#include <vector>
#include <utility>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using Graph = std::vector<std::vector<int>>;
using Edge = std::pair<int,int>;

class bipartite_matching{
  std::vector<bool> used_;
  std::vector<int> match_;
  bool dfs(const Graph &,int);
  public:
  int operator()(const Graph &g,int,std::vector<Edge> &edges);
};

class strong_connection_components{
  std::vector<int> vs_;
  std::vector<bool> visited_;
  public:
  void operator()(const Graph &,std::vector<int> &labels);
  void dfs(const Graph &,int);
  void rdfs(const Graph &,int,int,std::vector<int> &labels);
};

class dm_decomposition{
  Graph graph_;
  Eigen::SparseMatrix<double> mat_;
  std::vector<int> row_perm , col_perm;

  public:
  void operator()(const Eigen::SparseMatrix<double> &A);
};

#endif
