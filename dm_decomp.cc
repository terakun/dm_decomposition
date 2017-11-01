#include <iostream>
#include <iomanip>
#include <map>
#include "./dm_decomp.h"

void dm_decomposition::operator()(const Eigen::SparseMatrix<double> &A){
  mat_ = A;
  if(mat_.rows() != mat_.cols()){
    std::cout << "the matrix doesn't have inverse mat" << std::endl;
    return ;
  }

  int n = mat_.rows();
  std::map<int,int> row_nodes,col_nodes;
  Graph graph(2*n);
  for(int i=0;i<mat_.outerSize();++i){
    for(Eigen::SparseMatrix<double>::InnerIterator it(mat_,i);it;++it){
      if(it.value()!=0){
        graph[it.row()].push_back(it.col()+n);
        row_nodes[it.row()]++;
        col_nodes[it.col()]++;
      }
    }
  }

  // 最大マッチングを求める
  // [0,n) 行のノード
  // [n,2n) 列のノード
  // 2n ソース
  // 2n+1 シンク
  int src_node = 2*n; 
  int sink_node = src_node+1;
  
  Graph flow = graph;
  flow.emplace_back();
  flow.emplace_back();
  // ダミーのソース ( Ford Fulkerson法を使うために設定 )
  for(const auto &r:row_nodes){
    flow[src_node].push_back(r.first);
  }
  // ダミーのシンク
  for(const auto &c:col_nodes){
    flow[c.first+n].push_back(sink_node);
  }

  bipartite_matching bm;
  std::vector<Edge> edges;
  bm(flow,n,edges);

   // 完全マッチングでなければ終了（解が一意に定まらない）
  if(edges.size() != row_nodes.size()){
    std::cout << "Incomplete matching!" << std::endl;
    return ;
  }
 
  // 完全マッチングの部分は無向辺にする
  for(const auto &e : edges){
    if(e.first!=src_node&&e.second!=sink_node){
      std::cout << e.first << " " << e.second-n << std::endl;
      graph[e.second].push_back(e.first);
    }
  }

  // 強連結成分分解
  std::vector<int> labels;
  strong_connection_components scc;
  scc(graph,labels);
 
  Eigen::MatrixXd mat_dense = Eigen::MatrixXd(mat_);

  std::cout << "   " ;
  for(int j=0;j<n;++j){
    std::cout << std::setw(3) << labels[j+n];
  }
  std::cout << std::endl;

  for(int i=0;i<n;++i){
    std::cout << std::setw(3) << labels[i];
    for(int j=0;j<n;++j){
      std::cout << std::setw(3) << mat_dense(i,j);
    }
    std::cout << std::endl;
  }

/*
  std::vector<std::pair<int,int>> idx_and_labels;
  for(int i=0;i<labels.size();++i){
    idx_and_labels.emplace_back(i,labels[i]);
  }
  std::sort(idx_and_labels.begin(),idx_and_labels.end());

  std::cout << "Blocked upper triangle matrix : " << std::endl;*/
  // ブロック上三角行列に分解する
}

bool bipartite_matching::dfs(const Graph &g,int v){
  used_[v] = true;
  for(const auto &u:g[v]){
    int w = match_[u];
    if(w<0||!used_[w]&&dfs(g,w)){
      match_[v] = u;
      match_[u] = v;
      return true;
    }
  }
  return false;
}

int bipartite_matching::operator()(const Graph &g,int n,std::vector<Edge> &edges){
  int res = 0;
  match_.resize(g.size());
  used_.resize(g.size());
  std::fill(match_.begin(),match_.end(),-1);
  for(int v = 0; v < g.size();++v){
    std::fill(used_.begin(),used_.end(),false);
    if(dfs(g,v)){
      res++;
    }
  }

  for(int i=0;i<n;++i){
    if(match_[i]>=0) edges.emplace_back(i,match_[i]);
  }

  return res;
}

void strong_connection_components::operator()(const Graph &g,std::vector<int> &labels){
  int n = g.size();

  Graph rev_g(n);
  // reverse graphを作成
  for(int i=0;i<n;++i){
    for(const auto &n:g[i]){
      rev_g[n].push_back(i);
    }
  }
  
  visited_.resize(n);
  // gに対して深さ優先探索
  std::fill(visited_.begin(),visited_.end(),false);
  for(int i=0;i<n;++i){
    if(!visited_[i]) dfs(g,i);
  }
  std::fill(visited_.begin(),visited_.end(),false);
  // rev_gに対して深さ優先探索
  int l = 0;
  labels.resize(n);
  for(int i=vs_.size()-1;i>=0;--i){
    if(!visited_[vs_[i]]) rdfs(rev_g,vs_[i],l++,labels);
  }
}

void strong_connection_components::dfs(const Graph &g,int i){
  visited_[i] = true;
  for(const auto &n:g[i]){
    if(!visited_[n]) dfs(g,n);
  }
  vs_.push_back(i);
}

void strong_connection_components::rdfs(const Graph &rev_g,int i,int l,std::vector<int> &labels){
  visited_[i] = true;
  labels[i] = l;
  for(const auto &n:rev_g[i]){
    if(!visited_[n]) rdfs(rev_g,n,l,labels);
  }
}

void dm_decomposition::solve(const Eigen::VectorXd &b,Eigen::VectorXd &x){

}
