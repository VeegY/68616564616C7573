#include "vtkFileWriter.hpp"
#include <cstdlib>
int main()
{
    vtkFileWriter file ("ergebnisse/writertest", "test", 3, 3, 3, 2);  //Ordner ergebnisse muss existieren 
    float* data;
    float* data2;
    float* data6;
    data=new float[27];
    data2=new float[27];
    data6=new float[27];
    Vector<double> data3(27);
    Vector<double> data4(27);
    Vector<float> data5(8);
    Vector<float> vecx(8);
    Vector<float> vecy(8);
    Vector<float> vecz(8);
    for (int i=0; i<27; ++i){
        data[i]=static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        data3[i]=static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
        data2[i]=static_cast <float> (i);
        data4[i]=static_cast<double> (i);
    }
    for (int i=0; i<8; ++i){
        vecx[i]=0;
        vecy[i]=1;
        vecz[i]=2;
    }
    string name("rand");
    file.addPointDataToTimestep(data, 27, 1, name);
    file.addPointDataToTimestep(data2, 27, 1, "point2");
    file.addCellDataToAll(data5, "celldata");
    file.addPointDataToTimestep(data4, 0, "point2");
    file.addPointDataToTimestep(data3, 0, name);
    file.addPointVecToAll(data,data2,data6,27,"pointvec");
    file.addCellVecToAll(vecx,vecy,vecz,"cellvec");

    delete data;
    delete data2;
    delete data6;
    return 0;

}
