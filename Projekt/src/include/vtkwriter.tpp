namespace Icarus
{

template<typename type>
void vtkWriter::addPointDataToTimestep(const type data[], const int length, const unsigned timestep, std::string name)
{
    if (0 > timestep || timestep >= _tsteps) LOG_ERROR("invalid timestep");
    else
    {
        if (length==_points)
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
void vtkWriter::addPointDataToTimestep(const FullVector<type>& data, const unsigned timestep, const std::string name)
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
void vtkWriter::addPointVecToTimestep(const FullVector<type>& datax, const FullVector<type>& datay, const FullVector<type>& dataz, const unsigned timestep, const std::string name)
{
	if (0 > timestep || timestep >= _tsteps) LOG_ERROR("invalid timestep");
	else
	{
		if (datax.get_dim() == _points && datay.get_dim() == _points && dataz.get_dim() == _points)
		{
			std::ofstream file;
			std::string hfname(_filename);
			hfname.append(".vtk.");
			hfname.append(std::to_string(timestep));
			file.open(hfname.c_str(), std::ios::out | std::ios::app);
			if (file.is_open())
			{
				if (_point_data_written_last[timestep] == false)
				{
					file << "POINT_DATA " << _points << std::endl;
				}
				file << "VECTORS " << name << " float" << std::endl;
				file << "LOOKUP_TABLE default" << std::endl;
				for (int i = 0; i < _points; ++i)
				{
					file << static_cast<float>(datax[i]) << " " << static_cast<float>(datay[i]) << " " << static_cast<float>(dataz[i]) << std::endl;
				}
				file << std::endl;
				file.close();
				_point_data_written_last[timestep] = true;
				_cell_data_written_last[timestep] = false;
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
void vtkWriter::addPointVecToTimestep(const type datax[], const type datay[],
                                      const type dataz[], const int length,
                                      const unsigned timestep, std::string name)
{
	if (0 > timestep || timestep >= _tsteps) LOG_ERROR("invalid timestep");
	else
	{
		if (length == _points)
		{
			std::ofstream file;
			std::string hfname(_filename);
			hfname.append(".vtk.");
			hfname.append(std::to_string(timestep));
			file.open(hfname.c_str(), std::ios::out | std::ios::app);
			if (file.is_open())
			{
				if (_point_data_written_last[timestep] == false)
				{
					file << "POINT_DATA " << _points << std::endl;
				}
				file << "VECTORS " << name << " float" << std::endl;
				file << "LOOKUP_TABLE default" << std::endl;
				for (int i = 0; i < _points; ++i)
				{
					file << static_cast<float>(datax[i]) << " " << static_cast<float>(datay[i]) << " " << static_cast<float>(dataz[i]) << std::endl;
				}
				file << std::endl;
				file.close();
				_point_data_written_last[timestep] = true;
				_cell_data_written_last[timestep] = false;
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
void vtkWriter::addCellDataToTimestep(const FullVector<type>& data,
                                       const unsigned timestep, const std::string name)
{
	if (0 > timestep || timestep >= _tsteps) LOG_ERROR("invalid timestep");
    else
    {
        if (data.get_dim()==_cells)
        {
            std::ofstream file;
            std::string hfname(_filename);
            hfname.append(".vtk.");
            hfname.append(std::to_string(timestep));
            file.open(hfname.c_str(), std::ios::out | std::ios::app);
            if (file.is_open())
            {
                if (_cell_data_written_last[timestep]==false)
                {
                    file << "CELL_DATA " << _cells <<std::endl;
                }
                file << "SCALARS " << name <<" float"<<std::endl;
                file << "LOOKUP_TABLE default" << std::endl;
                for (int i=0; i<_cells; ++i)
                {
                    file << static_cast<float>(data[i])<<std::endl;
                }
                file << std::endl;
                file.close();
                _point_data_written_last[timestep]=false;
                _cell_data_written_last[timestep]=true;
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
void vtkWriter::addCellDataToTimestep(const type data[], const int length,
                                   const unsigned timestep, const std::string name)
{
	if (0 > timestep || timestep >= _tsteps) LOG_ERROR("invalid timestep");
    else
    {
        if (length==_cells)
        {
            std::ofstream file;
            std::string hfname(_filename);
            hfname.append(".vtk.");
            hfname.append(std::to_string(timestep));
            file.open(hfname.c_str(), std::ios::out | std::ios::app);
            if (file.is_open())
            {
                if (_cell_data_written_last[timestep]==false)
                {
                    file << "CELL_DATA " << _cells <<std::endl;
                }
                file << "SCALARS " << name <<" float"<<std::endl;
                file << "LOOKUP_TABLE default" << std::endl;
                for (int i=0; i<_cells; ++i)
                {
                    file << static_cast<float>(data[i])<<std::endl;
                }
                file << std::endl;
                file.close();
                _point_data_written_last[timestep]=false;
                _cell_data_written_last[timestep]=true;
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
void vtkWriter::addCellVecToTimestep(const FullVector<type>& datax,
                                      const FullVector<type>& datay,
                                      const FullVector<type>& dataz,
                                      const unsigned timestep,
                                      const std::string name)
{
	if (0 > timestep || timestep >= _tsteps) LOG_ERROR("invalid timestep");
	else
	{
		if (datax.get_dim() == _cells && datay.get_dim() == _cells && dataz.get_dim() == _cells)
		{
			std::ofstream file;
			std::string hfname(_filename);
			hfname.append(".vtk.");
			hfname.append(std::to_string(timestep));
			file.open(hfname.c_str(), std::ios::out | std::ios::app);
			if (file.is_open())
			{
				if (_cell_data_written_last[timestep] == false)
				{
					file << "CELL_DATA " << _cells << std::endl;
				}
				file << "VECTORS " << name << " float" << std::endl;
				file << "LOOKUP_TABLE default" << std::endl;
				for (int i = 0; i < _cells; ++i)
				{
					file << static_cast<float>(datax[i]) << " " << static_cast<float>(datay[i]) << " " << static_cast<float>(dataz[i]) << std::endl;
				}
				file << std::endl;
				file.close();
				_point_data_written_last[timestep] = false;
				_cell_data_written_last[timestep] = true;
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
void vtkWriter::addCellVecToTimestep(const type datax[], const type datay[],
                                      const type dataz[], const int length,
                                      const unsigned timestep, std::string name)
{
	if (0 > timestep || timestep >= _tsteps) LOG_ERROR("invalid timestep");
	else
	{
		if (length == _cells)
		{
			std::ofstream file;
			std::string hfname(_filename);
			hfname.append(".vtk.");
			hfname.append(std::to_string(timestep));
			file.open(hfname.c_str(), std::ios::out | std::ios::app);
			if (file.is_open())
			{
				if (_cell_data_written_last[timestep] == false)
				{
					file << "CELL_DATA " << _cells << std::endl;
				}
				file << "VECTORS " << name << " float" << std::endl;
				file << "LOOKUP_TABLE default" << std::endl;
				for (int i = 0; i < _cells; ++i)
				{
					file << static_cast<float>(datax[i]) << " " << static_cast<float>(datay[i]) << " " << static_cast<float>(dataz[i]) << std::endl;
				}
				file << std::endl;
				file.close();
				_point_data_written_last[timestep] = false;
				_cell_data_written_last[timestep] = true;
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
void vtkWriter::addPointDataToAll(const type data[],
                                  const int length, std::string name)
{
    for (unsigned i(0); i<_tsteps; ++i )
    {
        addPointDataToTimestep(data, length, i, name);
    }
}


template<typename type>
void vtkWriter::addPointDataToAll(const FullVector<type>& data,
                                  const std::string name)
{
    for (unsigned i(0); i<_tsteps; ++i)
    {
        addPointDataToTimestep(data, i, name);
    }
}


template<typename type>
void vtkWriter::addCellDataToAll(const type data[],
                                 const int length,
                                 std::string name)
{
    for (unsigned i(0); i<_tsteps; ++i )
    {
        addCellDataToTimestep(data, length, i, name);
    }
}


template<typename type>
void vtkWriter::addCellDataToAll(const FullVector<type>& data,
                                 const std::string name)
{
    for (unsigned i(0); i<_tsteps; ++i)
    {
        addCellDataToTimestep(data, i, name);
    }
}


template<typename type>
void vtkWriter::addPointVecToAll(const FullVector<type>& datax,
                                const FullVector<type>& datay,
                                const FullVector<type>& dataz,
                                const std::string name)

{
    for (unsigned i(0); i<=_tsteps; i++)
    {
        addPointVecToTimestep(datax, datay, dataz, i, name);
    }
}


template<typename type>
void vtkWriter::addPointVecToAll(const type datax[],
                                const type datay[],
                                const type dataz[],
                                int length,
                                std::string name)
{
    for (unsigned i(0); i<=_tsteps; i++)
    {
        addPointVecToTimestep(datax, datay, dataz, length, i, name);
    }
}


template<typename type>
void vtkWriter::addCellVecToAll(const FullVector<type>& datax,
                                const FullVector<type>& datay,
                                const FullVector<type>& dataz,
                                const std::string name)

{
    for (unsigned i(0); i<=_tsteps; i++)
    {
        addCellVecToTimestep(datax, datay, dataz, i, name);
    }
}


template<typename type>
void vtkWriter::addCellVecToAll(const type datax[],
                                const type datay[],
                                const type dataz[],
                                int length,
                                std::string name)
{
    for (unsigned i(0); i<=_tsteps; i++)
    {
        addCellVecToTimestep(datax, datay, dataz, length, i, name);
    }
}

}
