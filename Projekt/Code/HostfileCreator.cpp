#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
using namespace std;

int main()
{    
     //remove("txtFiles/network.txt");      
     string line,command="ls";     //sudo apt-get scan >> network.txt
     system(command.c_str());      //Nur Linux-freundlich!
     
     ifstream network; 
     ofstream hostfile;
     hostfile.open("txtFiles/hostfile.txt",ios::trunc);
     network.open("txtFiles/network.txt");
     if(network.is_open() && hostfile.is_open())
     {
       
          while(getline(network,line))
          {    
               
               if(line.compare(0,3,"192")==0)
               {
                    hostfile << line << " slots=" << "n\n";  
                    //n Prozesse, Prozesse als Array angeben?
                    //Der Hostf.Cr. muss eigentlich wissen was auf 
                    //wie vielen Boards berrechnet werden soll, danach
                    //kann eine effektive Aufteilung erfolgen, Angestrebt 
                    //aber 5 Prozesse pro Board, 4 GPU 1CPU.
               }
               
          }
          hostfile.close();
          network.close();
     }
     else{cout << "Critical Error regarding hostfile.txt" << endl;}
     return 0;
}

