#include <vector>

double assemblyRHSLoad(std::vector<int>& e, std::vector<int>& A){
   int n = e.size();
   double RHS;

   for(int i = 0; i < n; i++){
      //getQuadrature(e[i], "Name") = X, Y, Z, weigth;
      get_quadrature_xpoints(e[i], X, h);
      int nqp = X.size();

      for(int q = 0; q<nqp; q++){
         RHS += evaluate_Basis(e[i], A[i], X[q], Y[q], Z[q]) * f.eval(X[q], Y[q], Z[q]) * weight[q];
      }
   }
   return RHS;
}
