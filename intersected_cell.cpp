#include "intersected_cell.h"
#include <iostream>

namespace Reservoir {
    namespace WellIndexCalculation {

        Vector3d IntersectedCell::xvec() const {
            return corners()[5] - corners()[4];
        }

        Vector3d IntersectedCell::yvec() const {
            return corners()[6] - corners()[4];
        }

        Vector3d IntersectedCell::zvec() const {
            return corners()[0] - corners()[4];
        }

        double IntersectedCell::dx() const {
            return xvec().norm();
        }

        double IntersectedCell::dy() const {
            return yvec().norm();
        }

        double IntersectedCell::dz() const {
            return zvec().norm();
        }

        Vector3d IntersectedCell::get_segment_entry_point(int segment_index) const{
        	return entry_points_[segment_index];
        }

        Vector3d IntersectedCell::get_segment_exit_point(int segment_index) const{
        	return exit_points_[segment_index];
        }

        double IntersectedCell::get_segment_radius(int segment_index) const{
        	return segment_radius_[segment_index];
        }

        int IntersectedCell::num_segments() const{
        	return entry_points_.size();
        }

        void IntersectedCell::add_new_segment(Vector3d entry_point, Vector3d exit_point, double radius){
        	entry_points_.push_back(entry_point);
        	exit_points_.push_back(exit_point);
        	segment_radius_.push_back(radius);
        }

        double IntersectedCell::cell_well_index() const {
            return well_index_;
        }

        void IntersectedCell::set_cell_well_index(double well_index) {
            well_index_ = well_index;
        }

        void IntersectedCell::set_segment_calculation_data(int segment_index, std::string name, double value)
        {
        	// Check if this name already exists
        	std::map<std::string, std::vector<double>>::iterator it = calculation_data_.find(name);
        	if(it != calculation_data_.end())
        	{
        		if (segment_index >= 0 && segment_index < calculation_data_[name].size())
        		{
        			calculation_data_[name].at(segment_index) = value;
        		}
            	else if(segment_index == calculation_data_[name].size()){
            		calculation_data_[name].push_back(value);
            	}
            	else std::runtime_error("This segment index is out of bounds.");
        	}
        	else
        	{
        		calculation_data_[name].push_back(value);
        	}
        }

        std::map<std::string, std::vector<double>>& IntersectedCell::get_calculation_data()
		{
        	return calculation_data_;
		}

        int IntersectedCell::GetIntersectedCellIndex( std::vector<IntersectedCell> &cells, Grid::Cell grdcell){
        	if (cells.size() == 0)
        	{
        		cells.push_back(IntersectedCell(grdcell));
        		return 0;
        	}
        	else
        	{
        		for(int cell_index = 0 ; cell_index < cells.size(); cell_index++)
        		{
        			if (cells.at(cell_index).global_index() == grdcell.global_index())
        			{
        				return cell_index;
        			}
        		}

        		cells.push_back(IntersectedCell(grdcell));
        		return cells.size() - 1;
        	}
        }
    }
}
