#include "RcppEigen.h"


// [[Rcpp::depends(RcppEigen)]]

// distribution for the dfcop
//
// [[Rcpp::export]]
Eigen::VectorXd dfcop_pdf_cpp(const Eigen::VectorXd& w, 
                              const Eigen::MatrixXd& pv,
                              const Eigen::MatrixXd& x)
{
  int n = x.rows();
  int d = x.cols();
  int m = w.size();
  
  Eigen::VectorXd result(n);
  Eigen::VectorXd ones = Eigen::VectorXd::Ones(m);
  Eigen::RowVectorXd wt = w.transpose();

  Eigen::MatrixXd tmp(pv);
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < d; j++) {
      if(x(i,j) == 1) {
        tmp.col(j) = pv.col(j);
      } else {
        tmp.col(j) = ones - pv.col(j);
      }
    }
    result(i) = wt * tmp.rowwise().prod();
  }
 
  return result;
}
