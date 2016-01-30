#include "include/vtkwriter.hpp"

namespace Icarus
{

	vtkWriter::vtkWriter(std::string filename, std::string title, int xdim, int ydim, int zdim, unsigned int timesteps) :
		_filename(filename),
		_title(title),
		_tsteps(timesteps),
		_xdim(xdim),
		_ydim(ydim),
		_zdim(zdim)
	{
		_cells = ((_xdim - 1)*(_ydim - 1)*(_zdim - 1));
		_points = (_xdim*_ydim*_zdim);
		_cell_data_written_last = new bool[_tsteps]();
		_point_data_written_last = new bool[_tsteps]();
		std::ofstream file;
		std::string hfname(_filename);
		for (unsigned i = 0; i < _tsteps; ++i)
		{
			_point_data_written_last[i] = false;
			_cell_data_written_last[i] = false;
			hfname = _filename;
			hfname.append(".vtk.");
			hfname.append(std::to_string(i));
			file.open(hfname.c_str(), std::ios::out | std::ios::trunc);
			if (file.is_open())
			{
				file << "# vtk DataFile Version 3.0\n";
				file << title << std::endl;
				file << "ASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS ";
				file << xdim << " " << ydim << " " << zdim << std::endl;
				file << "SPACING 1 1 1\nORIGIN 0 0 0\n\n";
				file.close();
			}
			else LOG_ERROR("Unable to open file.");
		}
	}

	vtkWriter::vtkWriter(std::string filename, std::string title, int xdim, int ydim, int zdim, double h, unsigned int timesteps):
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
                file << "SPACING" << h << " " << h << " " << h << "\nORIGIN 0 0 0\n\n";
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
}

