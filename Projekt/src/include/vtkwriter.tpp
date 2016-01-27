namespace Icarus
{

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
