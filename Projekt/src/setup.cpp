#include "include/assemblefem.hpp"

namespace Icarus
{

int assembleFem::setup_A(int row, std::vector<char>& disc_points)
{
    int length(0);
    int numairelements(0);
    std::vector<int> ubf(27, 0); // used basis functions
    _A.resize(8);
    if (disc_points[row + z + y + 1] == 'a')
    {
        _A[numairelements] = 0;
        ++numairelements;
        ubf[13] = 1;
        ubf[14] = 1;
        ubf[16] = 1;
        ubf[17] = 1;
        ubf[22] = 1;
        ubf[23] = 1;
        ubf[25] = 1;
        ubf[26] = 1;
    }
    if (disc_points[row + z + y - 1] == 'a')
    {
        _A[numairelements] = 1;
        ++numairelements;
        ubf[12] = 1;
        ubf[13] = 1;
        ubf[15] = 1;
        ubf[16] = 1;
        ubf[21] = 1;
        ubf[22] = 1;
        ubf[24] = 1;
        ubf[25] = 1;
    }
    if (disc_points[row + z - y + 1] == 'a')
    {
        _A[numairelements] = 2;
        ++numairelements;
        ubf[10] = 1;
        ubf[11] = 1;
        ubf[13] = 1;
        ubf[14] = 1;
        ubf[19] = 1;
        ubf[20] = 1;
        ubf[22] = 1;
        ubf[23] = 1;
    }
    if (disc_points[row + z - y - 1] == 'a')
    {
        _A[numairelements] = 3;
        ++numairelements;
        ubf[9] = 1;
        ubf[10] = 1;
        ubf[12] = 1;
        ubf[13] = 1;
        ubf[18] = 1;
        ubf[19] = 1;
        ubf[21] = 1;
        ubf[22] = 1;
    }
    if (disc_points[row - z + y + 1] == 'a')
    {
        _A[numairelements] = 4;
        ++numairelements;
        ubf[4] = 1;
        ubf[5] = 1;
        ubf[7] = 1;
        ubf[8] = 1;
        ubf[13] = 1;
        ubf[14] = 1;
        ubf[16] = 1;
        ubf[17] = 1;
    }
    if (disc_points[row - z + y - 1] == 'a')
    {
        _A[numairelements] = 5;
        ++numairelements;
        ubf[3] = 1;
        ubf[4] = 1;
        ubf[6] = 1;
        ubf[7] = 1;
        ubf[12] = 1;
        ubf[13] = 1;
        ubf[15] = 1;
        ubf[16] = 1;
    }
    if (disc_points[row - z - y + 1] == 'a')
    {
        _A[numairelements] = 6;
        ++numairelements;
        ubf[1] = 1;
        ubf[2] = 1;
        ubf[4] = 1;
        ubf[5] = 1;
        ubf[10] = 1;
        ubf[11] = 1;
        ubf[13] = 1;
        ubf[14] = 1;
    }
    if (disc_points[row - z - y - 1] == 'a')
    {
        _A[numairelements] = 7;
        ++numairelements;
        ubf[0] = 1;
        ubf[1] = 1;
        ubf[3] = 1;
        ubf[4] = 1;
        ubf[9] = 1;
        ubf[10] = 1;
        ubf[12] = 1;
        ubf[13] = 1;
    }
    _A.resize(numairelements);
    for (int i(0); i < 27; ++i)
        length += ubf[i];
/*
   ohl-------ohr
   /         /
ovl-------ovr


   uhl-------uhr
   /         /
uvl-------uvr
*/
    return length;
}

void assembleFem::setup_e(int row)
{
    int n(_A.size());
    _e.resize(n);

    for(int i(0); i < n; i++)
    {
        switch(_A[i])
        {
        case 0: _e[i] = row;
            break;
        case 1: _e[i] = row-1;
            break;
        case 2: _e[i] = row-y;
            break;
        case 3: _e[i] = row-y-1;
            break;
        case 4: _e[i] = row-z;
            break;
        case 5: _e[i] = row-z-1;
            break;
        case 6: _e[i] = row-z-y;
            break;
        case 7: _e[i] = row-z-y-1;
            break;
        }
    }
}

void assembleFem::setup_neumann(int row, int plane, bool rigthbacktop, std::vector<char>& disc_points)
{
    _A.clear();
    _e.clear();

    if(!rigthbacktop) //postive normal gradient
    {
        switch(plane)
        {
        case 1:
            if(disc_points[row+z+y+1] == 'a' && disc_points[row+y+1] == 'b')
            {
                _e.push_back(row);
                _A.push_back(0);
            }
            if(disc_points[row+z+y-1] == 'a' && disc_points[row+y] == 'b')
            {
                _e.push_back(row -1);
                _A.push_back(1);
            }
            if(disc_points[row+z-y+1] == 'a' && disc_points[row+1] == 'b')
            {
                _e.push_back(row -y);
                _A.push_back(2);
            }
            if(disc_points[row+z-y-1] == 'a' && disc_points[row] == 'b')
            {
                _e.push_back(row -y-1);
                _A.push_back(3);
            }
            break;
        case 2:
            if(disc_points[row+z+y+1] == 'a' && disc_points[row+z+1] == 'b')
            {
                _e.push_back(row);
                _A.push_back(0);
            }
            if(disc_points[row+z+y-1] == 'a' && disc_points[row+z] == 'b')
            {
                _e.push_back(row -1);
                _A.push_back(1);
            }
            if(disc_points[row-z+y+1] == 'a' && disc_points[row+1] == 'b')
            {
                _e.push_back(row -z);
                _A.push_back(4);
            }
            if(disc_points[row-z+y-1] == 'a' && disc_points[row] == 'b')
            {
                _e.push_back(row -z-1);
                _A.push_back(5);
            }
            break;
        case 3:
            if(disc_points[row+z+y+1] == 'a' && disc_points[row+z+y] == 'b')
            {
                _e.push_back(row);
                _A.push_back(0);
            }
            if(disc_points[row+z-y+1] == 'a' && disc_points[row+z] == 'b')
            {
                _e.push_back(row -y);
                _A.push_back(2);
            }
            if(disc_points[row-z+y+1] == 'a' && disc_points[row+y] == 'b')
            {
                _e.push_back(row -z);
                _A.push_back(4);
            }
            if(disc_points[row-z-y+1] == 'a' && disc_points[row] == 'b')
            {
                _e.push_back(row -z-y);
                _A.push_back(6);
            }
            break;
        default: assert(false);
        }
    }
    else //negativ normal gradient
    {
        switch(plane)
        {
        case 1:
            if(disc_points[row-z+y+1] == 'a' && disc_points[row+y+1] == 'b')
            {
                _e.push_back(row);
                _A.push_back(4);
            }
            if(disc_points[row-z+y-1] == 'a' && disc_points[row+y] == 'b')
            {
                _e.push_back(row -1);
                _A.push_back(5);
            }
            if(disc_points[row-z-y+1] == 'a' && disc_points[row+1] == 'b')
            {
                _e.push_back(row -y);
                _A.push_back(6);
            }
            if(disc_points[row-z-y-1] == 'a' && disc_points[row] == 'b')
            {
                _e.push_back(row -y-1);
                _A.push_back(7);
            }
            break;
        case 2:
            if(disc_points[row+z-y+1] == 'a' && disc_points[row+z+1] == 'b')
            {
                _e.push_back(row);
                _A.push_back(2);
            }
            if(disc_points[row+z-y-1] == 'a' && disc_points[row+z] == 'b')
            {
                _e.push_back(row -1);
                _A.push_back(3);
            }
            if(disc_points[row-z-y+1] == 'a' && disc_points[row+1] == 'b')
            {
                _e.push_back(row -z);
                _A.push_back(6);
            }
            if(disc_points[row-z-y-1] == 'a' && disc_points[row] == 'b')
            {
                _e.push_back(row -z-1);
                _A.push_back(7);
            }
            break;
        case 3:
            if(disc_points[row+z+y-1] == 'a' && disc_points[row+z+y] == 'b')
            {
                _e.push_back(row);
                _A.push_back(1);
            }
            if(disc_points[row+z-y-1] == 'a' && disc_points[row+z] == 'b')
            {
                _e.push_back(row -y);
                _A.push_back(3);
            }
            if(disc_points[row-z+y-1] == 'a' && disc_points[row+y] == 'b')
            {
                _e.push_back(row -z);
                _A.push_back(5);
            }
            if(disc_points[row-z-y-1] == 'a' && disc_points[row] == 'b')
            {
                _e.push_back(row -z-y);
                _A.push_back(7);
            }
            break;
        default: assert(false);
        }
    }
}

void assembleFem::setup_plane_rigthbacktop(int row, std::vector<char>& disc_points, std::vector<int>& planes, std::vector<bool>& rightbacktops)
{
    planes.clear();
    rightbacktops.clear();
    if ((disc_points[row-z-y-1] != 'a' && disc_points[row+z-y-1] == 'a')
     || (disc_points[row-z-y+1] != 'a' && disc_points[row+z-y+1] == 'a')
     || (disc_points[row-z+y-1] != 'a' && disc_points[row+z+y-1] == 'a')
     || (disc_points[row-z+y+1] != 'a' && disc_points[row+z+y+1] == 'a'))
    {
        planes.push_back(1);
        rightbacktops.push_back(false);
    }
    if ((disc_points[row-z-y-1] == 'a' && disc_points[row+z-y-1] != 'a')
     || (disc_points[row-z-y+1] == 'a' && disc_points[row+z-y+1] != 'a')
     || (disc_points[row-z+y-1] == 'a' && disc_points[row+z+y-1] != 'a')
     || (disc_points[row-z+y+1] == 'a' && disc_points[row+z+y+1] != 'a'))
    {
        planes.push_back(1);
        rightbacktops.push_back(true);
    }

    if ((disc_points[row-z-y-1] != 'a' && disc_points[row-z+y-1] == 'a')
     || (disc_points[row-z-y+1] != 'a' && disc_points[row-z+y+1] == 'a')
     || (disc_points[row+z-y-1] != 'a' && disc_points[row+z+y-1] == 'a')
     || (disc_points[row+z-y+1] != 'a' && disc_points[row+z+y+1] == 'a'))
    {
        planes.push_back(2);
        rightbacktops.push_back(false);
    }
    if ((disc_points[row-z-y-1] == 'a' && disc_points[row-z+y-1] != 'a')
     || (disc_points[row-z-y+1] == 'a' && disc_points[row-z+y+1] != 'a')
     || (disc_points[row+z-y-1] == 'a' && disc_points[row+z+y-1] != 'a')
     || (disc_points[row+z-y+1] == 'a' && disc_points[row+z+y+1] != 'a'))
    {
        planes.push_back(2);
        rightbacktops.push_back(true);
    }

    if ((disc_points[row-z-y-1] != 'a' && disc_points[row-z-y+1] == 'a')
     || (disc_points[row-z+y-1] != 'a' && disc_points[row-z+y+1] == 'a')
     || (disc_points[row+z-y-1] != 'a' && disc_points[row+z-y+1] == 'a')
     || (disc_points[row+z+y-1] != 'a' && disc_points[row+z+y+1] == 'a'))
    {
        planes.push_back(3);
        rightbacktops.push_back(false);
    }
    if ((disc_points[row-z-y-1] == 'a' && disc_points[row-z-y+1] != 'a')
     || (disc_points[row-z+y-1] == 'a' && disc_points[row-z+y+1] != 'a')
     || (disc_points[row+z-y-1] == 'a' && disc_points[row+z-y+1] != 'a')
     || (disc_points[row+z+y-1] == 'a' && disc_points[row+z+y+1] != 'a'))
    {
        planes.push_back(3);
        rightbacktops.push_back(true);
    }
}

}//namespace Icarus
