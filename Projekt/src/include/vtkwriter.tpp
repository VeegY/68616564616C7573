
#ifndef __VTKWRITER_TPP_
#define __VTKWRITER_TPP_

// nur für intellisense
#include "vtkwriter.hpp"

namespace Icarus
{

vtkWriter::vtkWriter(std::string filename, std::string title, int xdim, int ydim, int zdim, unsigned int timesteps):
    _filename(filename),
    _title(title),
    _tsteps(timesteps),
    _xdim(xdim),
    _ydim(ydim),
    _zdim(zdim)
    {
        _cells=((_xdim-1)*(_ydim-1)*(_zdim-1));
        _points=(_xdim*_ydim*_zdim);
        _cell_data_written_last=new bool[_tsteps]();
        _point_data_written_last=new bool[_tsteps]();
        std::ofstream file;
        std::string hfname(_filename);
        for (unsigned i=0; i<_tsteps; ++i)
        {
            _point_data_written_last[i]=false;
            _cell_data_written_last[i]=false;
            hfname=_filename;
            hfname.append(".vtk.");
            hfname.append(std::to_string(i));
            file.open(hfname.c_str(), std::ios::out | std::ios::trunc );
            if (file.is_open())
            {
                file << "# vtk DataFile Version 3.0\n";
                file << title << std::endl;
                file << "ASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS ";
                file << xdim << " " << ydim << " " << zdim << std::endl;
                file << "ASPECT_RATIO 1 1 1\nORIGIN 0 0 0\n\n";
                file.close();
            }
            else LOG_ERROR("Unable to open file.");
        }
    }
vtkWriter::~vtkWriter()
{
    delete[] _cell_data_written_last;
    delete[] _point_data_written_last;
}

template<typename type>
void vtkWriter::addPointDataToTimestep(const type data[], const int length, const int timestep, std::string name)
{
    if (0 > timestep || timestep >= _tsteps) LOG_ERROR("invalid timestep")
    else
    {
        if (length==_points)
        {
            std::ofstream file;
            std::string hfname(_filename);
            hfname.append(".vtk.");
            hfname.append(std::to_string(timestep));
            file.open(hfname.c_str(), ios::out | ios::app);
            if (file.is_open())
            {
                if (point_data_written_last[timestep]==false)
                {
                    file << "POINT_DATA " << _points <<std::endl;
                }
                file << "SCALARS " << name <<" float"<<std::endl;
                file << "LOOKUP_TABLE default" << std::endl;
                for (int i=0; i<_points; ++i)
                {
                    file << static_cast<float>(data[i])<<std::endl;
                }
                file << std::endl;
                file.close();
                point_data_written_last[timestep]=true;
                cell_data_written_last[timestep]=false;
            }
            else LOG_ERROR("Unable to open file.");
        }
        else
        {
            LOG_ERROR("invalid length of data");
        }
    }
}

template<typename type>
void vtkWriter::addPointDataToTimestep(const FullVector<type>& data,const unsigned timestep, const std::string name)
{
	if (0 > timestep || timestep >= _tsteps) LOG_ERROR("invalid timestep");
    else
    {
        if (data.get_dim()==_points)
        {
            std::ofstream file;
            std::string hfname(_filename);
            hfname.append(".vtk.");
            hfname.append(std::to_string(timestep));
            file.open(hfname.c_str(), std::ios::out | std::ios::app);
            if (file.is_open())
            {
                if (_point_data_written_last[timestep]==false)
                {
                    file << "POINT_DATA " << _points <<std::endl;
                }
                file << "SCALARS " << name <<" float"<<std::endl;
                file << "LOOKUP_TABLE default" << std::endl;
                for (int i=0; i<_points; ++i)
                {
                    file << static_cast<float>(data[i])<<std::endl;
                }
                file << std::endl;
                file.close();
                _point_data_written_last[timestep]=true;
                _cell_data_written_last[timestep]=false;
            }
            else LOG_ERROR("Unable to open file.");
        }
        else
        {
            LOG_ERROR("invalid length of data");
        }
    }
}

template<typename type>
void vtkWriter::addPointVecToTimestep(const FullVector<type>& datax, const FullVector<type>& datay, const FullVector<type>& dataz, const int timestep, const std::string name)
{
	if (0 > timestep || timestep >= _tsteps) LOG_ERROR("invalid timestep")
	else
	{
		if (datax.get_dim() == _points || datax.get_dim() == _points || datax.get_dim() == _points)
		{
			std::ofstream file;
			std::string hfname(_filename);
			hfname.append(".vtk.");
			hfname.append(std::to_string(timestep));
			file.open(hfname.c_str(), ios::out | ios::app);
			if (file.is_open())
			{
				if (point_data_written_last[timestep] == false)
				{
					file << "POINT_DATA " << _points << std::endl;
				}
				file << "SCALARS " << name << " float" << std::endl;
				file << "LOOKUP_TABLE default" << std::endl;
				for (int i = 0; i < _points; ++i)
				{
					file << static_cast<float>datax[i] << " " << datay[i] << " " << dataz[i] << std::endl;
				}
				file << std::endl;
				file.close();
				point_data_written_last[timestep] = true;
				cell_data_written_last[timestep] = false;
			}
			else LOG_ERROR("Unable to open file.");
		}
		else
		{
			LOG_ERROR("invalid length of data");
		}
	}
}

}
#endif // __VTKWRITER_HPP_
