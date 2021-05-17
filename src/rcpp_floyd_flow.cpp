#include <Rcpp.h>
using namespace Rcpp;

//' Floyd-Warshall algorithm for shortest-path distances and total flows
//'
//' @param C numeric matrix of edge costs
//' @return D matrix of shortest-path distances, i.e. sum of costs along minimising paths
//' @return F matrix of total flows, i.e. sum of inverse costs, along shortest-paths
// [[Rcpp::export]]
List rcpp_floyd_flow(NumericMatrix C) {
  List res;
  int n = C.nrow();
  NumericMatrix D(n, n);
  NumericMatrix F(n, n);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      D[i + n * j] = C[i + n * j];
      F[i + n * j] = 1 / C[i + n * j];
    }
  }
  for (int i = 0; i < n; i++) {
    D[i + n * i] = 0;              /* no self cycle */
    F[i + n * i] = 0;              /* no self cycle */
  }
  for (int k = 0; k <n; k++) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        if (D[i + n * k] + D[k + n * j] < D[i + n * j]) {
          D[i + n * j] = D[i + n * k] + D[k + n * j];
          F[i + n * j] = F[i + n * k] + F[k + n * j];
        }
      }
    }
  }
  res["D"] = D;
  res["F"] = F;
  return(res);
}
