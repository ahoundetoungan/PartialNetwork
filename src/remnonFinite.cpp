// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#define NDEBUG 1

using namespace Rcpp;
using namespace arma;
using namespace std;

//[[Rcpp::export]]
List rem_non_fin (const arma::mat& net) {
  arma::mat out = net;
  int n         = net.n_rows;
  arma::uvec id = arma::regspace<arma::uvec>(0, n - 1);
  arma::umat R  = arma::repmat(id, 1, n);
  arma::umat C  = arma::repmat(id.t(), n, 1);
  arma::uvec tp = arma::find_nan(out);
  while(tp.n_elem > 0){
    R             = R.elem(tp);
    C             = C.elem(tp);
    arma::uvec z  = arma::join_cols(R, C);
    arma::uvec zu = arma::sort(arma::unique(z));
    int nzu       = zu.n_elem;
    arma::vec czu(nzu);
    for(int s(0); s < nzu; ++ s){
      czu(s)      = sum(z == zu(s));
    }
    arma::uvec x  = zu.elem(arma::find(czu == czu.max()));
    id.shed_rows(x);
    out           = net.rows(id);
    out           = out.cols(id);
    tp            = arma::find_nan(out);
    n             = out.n_rows;
    arma::umat y  = arma::regspace<arma::uvec>(0, n - 1);
    R             = arma::repmat(y, 1, n);
    C             = arma::repmat(y.t(), n, 1);
  }
  return List::create(Named("net") = out, Named("id") = id + 1);
}
