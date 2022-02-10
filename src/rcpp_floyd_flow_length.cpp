#include <Rcpp.h>
using namespace Rcpp;

// Implementation of Floyd-Warshall algorithm for shortest-path distances (sum of costs)
// Additionally, compute also the total flow (sum of strengths) along those paths and
// the Euclidean length of shortest-paths (sum of Euclidean lengths) and
// the total Euclidean closeness (sum of reciprocals Euclidean lengths) along sh.-paths.

//' Floyd-Warshall algorithm for shortest-path distances and total flows
//'
//' Additionally, it computes also the Euclidean length of shortest-paths (sum of Euclidean lengths)
//' and the total Euclidean closeness (sum of reciprocals Euclidean lengths) along shortest-paths.
//' @param C numeric matrix of edge costs
//' @param L numeric matrix of edge (Euclidean) lenghts
//' @return D matrix of shortest-path distances, i.e. sum of costs along minimising paths
//' @return F matrix of total flows, i.e. sum of inverse costs, along shortest-paths
//' @return D_eucl matrix of Euclidean distances along shortest-paths (paths minimising sum of costs)
//' @return F_eucl matrix of total Eucl. flows, i.e. sum of inverse lengths, along shortest-paths
//' @export
// [[Rcpp::export]]
List rcpp_floyd_flow_length(NumericMatrix C, NumericMatrix L) {
  List res;
  int n = C.nrow();
  NumericMatrix D(n, n);
  NumericMatrix F(n, n);
  NumericMatrix D_eucl(n, n);
  NumericMatrix F_eucl(n, n);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      D[i + n * j] = C[i + n * j];
      F[i + n * j] = 1 / C[i + n * j];
      D_eucl[i + n * j] = L[i + n * j];
      F_eucl[i + n * j] = 1 / L[i + n * j];
    }
  }
  for (int i = 0; i < n; i++) {
    D[i + n * i] = 0;              /* no self cycle */
    F[i + n * i] = 0;              /* no self cycle */
    D_eucl[i + n * i] = 0;         /* no self cycle */
    F_eucl[i + n * i] = 0;         /* no self cycle */
  }
  for (int k = 0; k <n; k++) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        if (D[i + n * k] + D[k + n * j] < D[i + n * j]) {
          D[i + n * j] = D[i + n * k] + D[k + n * j];
          F[i + n * j] = F[i + n * k] + F[k + n * j];
          D_eucl[i + n * j] = D_eucl[i + n * k] + D_eucl[k + n * j];
          F_eucl[i + n * j] = F_eucl[i + n * k] + F_eucl[k + n * j];
        }
      }
    }
  }
  res["D"] = D;
  res["F"] = F;
  res["D_eucl"] = D_eucl;
  res["F_eucl"] = F_eucl;
  return(res);
}
