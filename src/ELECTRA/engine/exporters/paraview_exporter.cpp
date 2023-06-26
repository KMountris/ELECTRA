/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */


#include "ELECTRA/exporters/paraview_exporter.hpp"


namespace ELECTRA {

ParaviewExporter::ParaviewExporter()
{}


ParaviewExporter::~ParaviewExporter()
{}


void ParaviewExporter::CreateVtu(const WeakModel3D &weak_model_3d)
{
    // Store in the exporter the 3D weak formulation model for for output.
    this->weak_model_3d_ = weak_model_3d;
    
    // Create xml declaration for the document.
    tinyxml2::XMLDeclaration *declaration = this->output_.NewDeclaration();
    this->output_.InsertFirstChild(declaration);
    
    // Create the main xml branch.
    tinyxml2::XMLElement *vtk_file = this->output_.NewElement("VTKFile");
    vtk_file->SetAttribute("type", "UnstructuredGrid");
    vtk_file->SetAttribute("version", "0.1");
    vtk_file->SetAttribute("byte_order", "LittleEndian");
    this->output_.InsertEndChild(vtk_file);
    
    // Create UnstructuredGrid branch.
    tinyxml2::XMLElement *unstruct_grid = this->output_.NewElement("UnstructuredGrid");
    vtk_file->InsertEndChild(unstruct_grid);

    // Create Piece branch.
    tinyxml2::XMLElement *piece = this->output_.NewElement("Piece");
    piece->SetAttribute("NumberOfPoints", this->weak_model_3d_.Grid().NodesNum());
    piece->SetAttribute("NumberOfCells", this->weak_model_3d_.TetrahedraMeshElemsNum());
    unstruct_grid->InsertEndChild(piece);
    
    // Create Points branch.
    tinyxml2::XMLElement *points = this->output_.NewElement("Points");
    piece->InsertEndChild(points);
    
    // Create DataArray branch for points coordinates.
    tinyxml2::XMLElement *points_coords = this->output_.NewElement("DataArray");
    points_coords->SetAttribute("type", "Float64");
    points_coords->SetAttribute("Name", "Points");
    points_coords->SetAttribute("NumberOfComponents", 3);
    points_coords->SetAttribute("format", "ascii");
    
    // Generate string with the nodes coordinates.
    std::string coordinates = "\n\t\t\t\t\t\t\t\t\t";
    for (auto &node : this->weak_model_3d_.Grid().Nodes()) {
        coordinates += std::to_string(node.Coordinates().X());
        coordinates += " ";
        coordinates += std::to_string(node.Coordinates().Y());
        coordinates += " ";
        coordinates += std::to_string(node.Coordinates().Z());
        coordinates += "\n\t\t\t\t\t\t\t\t\t";
    }
    
    // Write the nodes coordinates and insert the branch in the xml tree.
    points_coords->SetText(coordinates.c_str());
    points->InsertEndChild(points_coords);
    
    // Create Cells branch.
    tinyxml2::XMLElement *cells = this->output_.NewElement("Cells");
    piece->InsertEndChild(cells);
    
    // Create Cells connectivity, offsets, and type DataArray branches.
    tinyxml2::XMLElement *cells_connectivity = this->output_.NewElement("DataArray");
    cells_connectivity->SetAttribute("type", "Int32");
    cells_connectivity->SetAttribute("Name", "connectivity");
    cells_connectivity->SetAttribute("format", "ascii");
    
    tinyxml2::XMLElement *cells_offsets = this->output_.NewElement("DataArray");
    cells_offsets->SetAttribute("type", "Int32");
    cells_offsets->SetAttribute("Name", "offsets");
    cells_offsets->SetAttribute("format", "ascii");
    
    tinyxml2::XMLElement *cells_types = this->output_.NewElement("DataArray");
    cells_types->SetAttribute("type", "UInt8");
    cells_types->SetAttribute("Name", "types");
    cells_types->SetAttribute("format", "ascii");

    // Generate strings with the connectivity, offsets, and types of the cells.
    std::string connectivity = "\n\t\t\t\t\t\t\t\t\t";
    std::string offsets = "\n\t\t\t\t\t\t\t\t\t";
    std::string types = "\n\t\t\t\t\t\t\t\t\t";
    
    for (auto &elem : this->weak_model_3d_.TetrahedralMesh().Elements()) {

        // The id of the current element.
        auto id = &elem - &this->weak_model_3d_.TetrahedralMesh().Elements()[0];
        
        connectivity += std::to_string(elem.N1()) + " " + std::to_string(elem.N2()) + " ";
        connectivity += std::to_string(elem.N3()) + " " + std::to_string(elem.N4()) + " ";
        offsets += std::to_string(4*id+4) + " ";
        types += std::to_string(10) + " ";
        
        // Add a new line every 10 elements.
        if ((id !=0) && !(id % 10)) {
            connectivity += "\n\t\t\t\t\t\t\t\t\t";
            offsets += "\n\t\t\t\t\t\t\t\t\t";
            types += "\n\t\t\t\t\t\t\t\t\t";
        }
    }

    // Add a new line in the end of the string if has not been added during the loop.
    if (!(this->weak_model_3d_.TetrahedraMeshElemsNum() % 5)) {
        connectivity += "\n\t\t\t\t\t\t\t\t\t";
        offsets += "\n\t\t\t\t\t\t\t\t\t";
        types += "\n\t\t\t\t\t\t\t\t\t";
    }
    
    //Add attributes to connectivity, offsets, and types branches.
    cells_connectivity->SetText(connectivity.c_str());
    cells_offsets->SetText(offsets.c_str());
    cells_types->SetText(types.c_str());
    
    //Assign connectivity, offsets, and types branches.
    cells->InsertEndChild(cells_connectivity);
    cells->InsertEndChild(cells_offsets);
    cells->InsertEndChild(cells_types);

}


//void ParaviewExporter::AddPointNormals()
//{
//    //Get the first node of the vtu file.
//    tinyxml2::XMLNode *parent_node = this->vtu_doc_.FirstChildElement("VTKFile");
//    if (parent_node == nullptr) {
//        std::string error = "Error: Unknown XML file given. First node name expected to be VTKFile ...";
//        throw std::runtime_error(error.c_str());
//    }
    
    
//    //Find the branch with name Piece.
//    tinyxml2::XMLElement *piece = parent_node->FirstChildElement("UnstructuredGrid")->FirstChildElement("Piece");
//    if (piece == nullptr) {
//        std::string error = "Error: Searched for XML node named 'Piece' failed. Expected hierarchy is: VTKFile->UnstructuredGrid->Piece ...";
//        throw std::runtime_error(error.c_str());
//    }

    
//    //Find the branch with name PointData.
//    tinyxml2::XMLElement *point_data = piece->FirstChildElement("PointData");
//    if (point_data == nullptr) {
//        std::string error = "Error: Searched for XML node named 'PointData' failed. Expected to be child of node 'Piece' ...";
//        throw std::runtime_error(error.c_str());
//    }
    
//    //Set normals attribute in point_data branch.
//    point_data->SetAttribute("Normals", "normals");
    
//    //Create point normals branch.
//    tinyxml2::XMLElement *point_normals = this->vtu_doc_.NewElement("DataArray");
//    point_normals->SetAttribute("type", "Float32");
//    point_normals->SetAttribute("Name", "normals");
//    point_normals->SetAttribute("NumberOfComponents", 3);
//    point_normals->SetAttribute("format", "ascii");
    
//    //String to write point normals information
//    std::string point_normals_text = "\n\t\t\t\t\t\t\t\t\t";
//    for (int i = 0; i != this->grid_.GetNumNodes(); ++i) {
//        point_normals_text += std::to_string(this->grid_.Nodes().at(i).Nx()) + " ";
//        point_normals_text += std::to_string(this->grid_.Nodes().at(i).Ny()) + " ";
//        point_normals_text += std::to_string(this->grid_.Nodes().at(i).Nz());
//        point_normals_text += "\n\t\t\t\t\t\t\t\t\t";
//    }
    
//    //Assign point normals text to the point normals branch.
//    point_normals->SetText(point_normals_text.c_str());

//    //Assign point normals branch as child of point data branch.
//    point_data->InsertEndChild(point_normals);

//}


void ParaviewExporter::AddScalarField(const Eigen::VectorXd& scalar_field, const std::string &scalar_field_name)
{
    // Check consistency of given vector field's size with model's grid.
    if (scalar_field.rows() != this->weak_model_3d_.Grid().NodesNum()) {
        std::string error = "[CLOUDEA Error] The size of the given vector field is different from the model's grid size.";
        throw std::invalid_argument(error.c_str());
    }
    
    // Get the parent branch of the exporting xml file.
    tinyxml2::XMLNode *parent_branch = this->output_.FirstChildElement("VTKFile");

    // Check if parent branch is the expected one (VTKFile).
    if (parent_branch == nullptr) {
        std::string error = "[CLOUDEA Error] Unknown XML file. Parent branch expected to be VTKFile.";
        throw std::runtime_error(error.c_str());
    }
    
    // Get the Piece branch.
    tinyxml2::XMLElement *piece = parent_branch->FirstChildElement("UnstructuredGrid")->FirstChildElement("Piece");

    // Check if Piece branch was found.
    if (piece == nullptr) {
        std::string error = "[CLOUDEA Error] 'Piece' branch was not found. Expected XML tree is: VTKFile->UnstructuredGrid->Piece.";
        throw std::runtime_error(error.c_str());
    }

    // Get the PointData branch.
    tinyxml2::XMLElement *point_data = piece->FirstChildElement("PointData");

    // Create PointData branch if does not exist.
    if (point_data == nullptr) {
        point_data = this->output_.NewElement("PointData");
        piece->InsertEndChild(point_data);
    }
    
    // Add vector field attribute.
    if (point_data->Attribute("Scalars") != nullptr) {
        // Update existing attribute.
        std::string old_attribute = point_data->FirstAttribute()->Value();
        std::string new_attribute = old_attribute + " " + scalar_field_name;
        point_data->SetAttribute("Scalars", new_attribute.c_str());
    }
    else {
        point_data->SetAttribute("Scalars", scalar_field_name.c_str());
    }
    
    // Create scalar field branch.
    tinyxml2::XMLElement *point_scalars = this->output_.NewElement("DataArray");
    point_scalars->SetAttribute("type", "Float64");
    point_scalars->SetAttribute("Name", scalar_field_name.c_str());
    point_scalars->SetAttribute("NumberOfComponents", 1);
    point_scalars->SetAttribute("format", "ascii");
    
    //Generate string with the values of the scalar field.
    std::string scalars = "\n\t\t\t\t\t\t\t\t\t";
    for (int i = 0; i != scalar_field.rows(); ++i) {
        scalars += std::to_string(scalar_field.coeff(i));
        scalars += "\n\t\t\t\t\t\t\t\t\t";
    }
    
    //Assign vector values and insert vector field branch in the xml tree.
    point_scalars->SetText(scalars.c_str());
    point_data->InsertEndChild(point_scalars);
}


void ParaviewExporter::AddVectorField(const Eigen::MatrixXd &vector_field, const std::string &vector_field_name)
{
    // Check consistency of given vector field's size with model's grid.
    if (vector_field.rows() != this->weak_model_3d_.Grid().NodesNum()) {
        std::string error = "[CLOUDEA Error] The size of the given vector field is different from the model's grid size.";
        throw std::invalid_argument(error.c_str());
    }
    
    // Get the parent branch of the exporting xml file.
    tinyxml2::XMLNode *parent_branch = this->output_.FirstChildElement("VTKFile");

    // Check if parent branch is the expected one (VTKFile).
    if (parent_branch == nullptr) {
        std::string error = "[CLOUDEA Error] Unknown XML file. Parent branch expected to be VTKFile.";
        throw std::runtime_error(error.c_str());
    }
    
    // Get the Piece branch.
    tinyxml2::XMLElement *piece = parent_branch->FirstChildElement("UnstructuredGrid")->FirstChildElement("Piece");

    // Check if Piece branch was found.
    if (piece == nullptr) {
        std::string error = "[CLOUDEA Error] 'Piece' branch was not found. Expected XML tree is: VTKFile->UnstructuredGrid->Piece.";
        throw std::runtime_error(error.c_str());
    }

    // Get the PointData branch.
    tinyxml2::XMLElement *point_data = piece->FirstChildElement("PointData");

    // Create PointData branch if does not exist.
    if (point_data == nullptr) {
        point_data = this->output_.NewElement("PointData");
        piece->InsertEndChild(point_data);
    }
    
    // Add vector field attribute.
    if (point_data->Attribute("Vectors") != nullptr) {
        // Update existing attribute.
        std::string old_attribute = point_data->FirstAttribute()->Value();
        std::string new_attribute = old_attribute + " " + vector_field_name;
        point_data->SetAttribute("Vectors", new_attribute.c_str());
    }
    else {
        point_data->SetAttribute("Vectors", vector_field_name.c_str());
    }
    
    // Create vector field branch.
    tinyxml2::XMLElement *point_vectors = this->output_.NewElement("DataArray");
    point_vectors->SetAttribute("type", "Float64");
    point_vectors->SetAttribute("Name", vector_field_name.c_str());
    point_vectors->SetAttribute("NumberOfComponents", 3);
    point_vectors->SetAttribute("format", "ascii");
    
    //Generate string with the values of the vector field.
    std::string vectors = "\n\t\t\t\t\t\t\t\t\t";
    for (int i = 0; i != vector_field.rows(); ++i) {
        vectors += std::to_string(vector_field.coeff(i, 0)) + " ";
        vectors += std::to_string(vector_field.coeff(i, 1)) + " ";
        vectors += std::to_string(vector_field.coeff(i, 2));
        vectors += "\n\t\t\t\t\t\t\t\t\t";
    }
    
    //Assign vector values and insert vector field branch in the xml tree.
    point_vectors->SetText(vectors.c_str());
    point_data->InsertEndChild(point_vectors);
    

    /////////////NOTE////////////////////
    //    //cell_data should be created too...
}


void ParaviewExporter::ClearVectorFields()
{
    // Get the PointData branch.
    tinyxml2::XMLElement *point_data = this->output_.FirstChildElement("VTKFile")->FirstChildElement("UnstructuredGrid")->
                                        FirstChildElement("Piece")->FirstChildElement("PointData");

    //Delete PointData branch if it exists.
    if (point_data != nullptr) {
        this->output_.DeleteNode(point_data);
    }

}


void ParaviewExporter::Export(const std::string &export_filename)
{    

    // Initialize the path of the exporting file.
    std::string path = "";

    // Position of the last slash in the exporting file's name.
    std::size_t last_slash = export_filename.find_last_of("/\\");

    // Get the path directory of the exporting file name.
    if (last_slash != std::string::npos) {
        path = export_filename.substr(0, last_slash);
    }

    // Create the path's directory if it doesn't exist.
    boost::filesystem::path dir(path);
    if (!path.empty() && !boost::filesystem::exists(dir)) {
        boost::filesystem::create_directories(dir);
    }

    // Initialize the exporting file name's extension.
    std::string ext = "";

    // Search for the extension.
    if (export_filename.find_last_of(".") != std::string::npos) {
        ext = export_filename.substr(export_filename.find_last_of("."));
    }

    //Output result with error check.
    tinyxml2::XMLError out_result;

    // Add .vtu extension before exporting if it's missing from the exporting file's name.
    if (ext != ".vtu") {
        std::string complete_filename = export_filename + ".vtu";
        out_result = this->output_.SaveFile(complete_filename.c_str());
    }
    else {
        out_result = this->output_.SaveFile(export_filename.c_str());
    }

    // Check if export was sucessful.
    if (out_result != tinyxml2::XML_SUCCESS) {
        std::string error = "[CLOUDEA ERROR] Unable to export the file: " + export_filename;
        throw std::runtime_error(error.c_str());
    }
}


void ParaviewExporter::CreatePvdAnimation(const std::string &vtu_files_dir, const int &files_number, const std::string &pvd_filename, double time_inc)
{
    // Check if the requested directory is available.
    if(!boost::filesystem::exists(vtu_files_dir) || !boost::filesystem::is_directory(vtu_files_dir)) {
        std::string error = "[CLOUDEA ERROR] Cannot generate ParaView animation file (.pvd). Expected path to an existing directory.";
        throw std::invalid_argument(error.c_str());
    }

    // Container to store the vtu files paths.
    std::map<int, std::string> vtu_files_paths;

    // Search vtu files in the directory.
    std::string file_id = "";

    for (auto &file : boost::make_iterator_range(boost::filesystem::directory_iterator(vtu_files_dir),{})){

        if(boost::filesystem::is_regular_file(file.path()) && (file.path().extension() == ".vtu")) {
            // Set the id of the currently readed file.
            file_id = file.path().string().substr(file.path().string().find_first_of("0123456789"));
            file_id = file_id.substr(0, file_id.find_first_of("."));

            // Check if file_id string is numerical.
            if (file_id.empty() || file_id.find_first_not_of("0123456789") != std::string::npos) {
                std::string error = "[CLOUDEA ERROR] Cannot generate Paraview animation file (.pvd). Expected file format filenameXX.vtu."
                                    " Where XX is a number.";
                throw std::runtime_error(error.c_str());
            }

            // Insert and auto-sort file paths in the container.
            vtu_files_paths.emplace(std::stoi(file_id), file.path().string());
        }
    }

    // Create pvd animation xml file.
    tinyxml2::XMLDocument animation;

    // Create xml declaration for the animation.
    tinyxml2::XMLDeclaration *declaration = animation.NewDeclaration();
    animation.InsertFirstChild(declaration);

    // Create the main xml branch.
    tinyxml2::XMLElement *vtk_file = animation.NewElement("VTKFile");
    vtk_file->SetAttribute("type", "Collection");
    vtk_file->SetAttribute("version", "0.1");
    vtk_file->SetAttribute("byte_order", "LittleEndian");
    animation.InsertEndChild(vtk_file);

    // Create Collection branch.
    tinyxml2::XMLElement *collection = animation.NewElement("Collection");

    // Generate vtu file paths string.
    double time_step = 0.;
    for (auto &file_path : vtu_files_paths) {
        // Add file in collection if its index is smaller than the requested number of files.
        if (file_path.first < files_number) {
            
            tinyxml2::XMLElement *dataset = animation.NewElement("DataSet");

            if (time_inc > 2*std::numeric_limits<double>::epsilon()) {
               time_step += time_inc; 
               dataset->SetAttribute("timestep", time_step);
            } else {
                 dataset->SetAttribute("timestep", file_path.first);
            }

            dataset->SetAttribute("group","");
            dataset->SetAttribute("part",0);
            dataset->SetAttribute("file",file_path.second.c_str());
            collection->InsertEndChild(dataset);
        }
    }
    vtk_file->InsertEndChild(collection);

    // Initialize the path of the animation file.
    std::string path = "";

    // Position of the last slash in the exporting file's name.
    std::size_t last_slash = pvd_filename.find_last_of("/\\");

    // Get the path directory of the exporting file name.
    if (last_slash != std::string::npos) {
        path = pvd_filename.substr(0, last_slash);
    }

    // Create the path's directory if it doesn't exist.
    boost::filesystem::path dir(path);
    if (!path.empty() && !boost::filesystem::exists(dir)) {
        boost::filesystem::create_directories(dir);
    }

    //Export result with error check.
    tinyxml2::XMLError out_result;

    // Add .pvd extension before exporting animation if it's missing from the file's name.
    std::string ext = pvd_filename.substr(pvd_filename.find_last_of("."));
    if (ext != ".pvd") {
        std::string complete_filename = pvd_filename + ".pvd";
        out_result = animation.SaveFile(complete_filename.c_str());
    }
    else {
        out_result = animation.SaveFile(pvd_filename.c_str());
    }

    // Check if export was sucessful.
    if (out_result != tinyxml2::XML_SUCCESS) {
        std::string error = "[CLOUDEA ERROR] Unable to export the ParaView animation file (.pvd).";
        throw std::runtime_error(error.c_str());
    }


}


} //end of namespace CLOUDEA

