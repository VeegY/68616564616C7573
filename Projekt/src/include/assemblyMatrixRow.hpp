#include <vector>

void assemblyMatrixRow(std::vector<int>& e, std::vector<int>& A)
{
   int n = e.size();
   std::vector<int> a(8);
   std::vector<double> value(27);
   std::vector<int> column(27);
   std::vector<double> grad_Basis1(3);
   std::vector<double> grad_Basis2(3);
   //std::vector<double> Zeile(Nx*Ny*Nz);

   double zwsp;

   a[0]=0; a[1]=1; a[2]= 1+y; a[3]=y; a[4]= z; a[5]=1+z; a[6]= 1+y+z; a[7]=y+z;

   //Schleife ueber alle Elemente
   for(int i = 0; i < n; i++){
      //getQuadrature(e[i], "Name") = X, Y, Z, weigth;
      get_quadrature_xpoints(e[i], X, h);
      int nqp = X.size();
      //Schleife ueber alle Ecken von e
      for(int B = 0; B < 8; B++){
         //Quadratur
         zwsp = 0;
         for(int q = 0; q<nqp; q++){
            grad_Basis1 = evaluate_gradient_Basis(e[i], A[i], X[q], Y[q], Z[q]);
            grad_Basis2 = evaluate_gradient_Basis(e[i], B, X[q], Y[q], Z[q]);
            zwsp += grad_Basis1.dot(grad_Basis2) * weigth[q];
         }

         int j=0;
         bool abbr = false;
         while(abbr==false||j<27){
            if(e[i] + a[B] == 0){
                column[j] = e[i] + a[B];
                value[j] += zwsp;
                abbr = true;
            }
            else if(e[i] + a[B] == column[j]){
                value[j] += zwsp;
                abbr = true;
            }
            j++;
         }
         //Zeile[e[i] + a[B]] += zwsp;
      }
   }

   //Zeile zu value und column

}
