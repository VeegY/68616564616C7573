#include<vector>

double assemblyRHSNeumann(std::vector<int>& e, std::vector<int>& A, int Ebene){
   int n = e.size();
   double RHS;
   for(int i = 0; i < n; i++){
      getQuadrature(e[i], "Name") = X, Y, Z, weigth;
      int nqp = X.size();
         for(int q = 0; q<nqp; q++){
            //X-Y-Ebene
            if(Ebene == 1){
               RHS += evaluate_2D_Basis(e[i], A[i], Ebene, X[q], Y[q]) * g.eval(X[q], Y[q], Z[q]) * weight[q];
            }
            //X-Z-Ebene
            if(Ebene == 2){
               RHS += evaluate_2D_Basis(e[i], A[i], Ebene, X[q], Z[q]) * g.eval(X[q], Y[q], Z[q]) * weight[q];
            }
            //Y-Z-Ebene
            if(Ebene == 3){
               RHS += evaluate_2D_Basis(e[i], A[i], Ebene, Y[q], Z[q]) * g.eval(X[q], Y[q], Z[q]) * weight[q];
            }
         }
   }
   return RHS;
}
