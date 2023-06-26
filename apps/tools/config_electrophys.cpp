/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

// #ifndef ELECTRA_APPS_TOOLS_CONFIG_GEO_TPP_
// #define ELECTRA_APPS_TOOLS_CONFIG_GEO_TPP_

#include "config_electrophys.hpp"

namespace APP_ELECTRA
{

ConfigElectrophys::ConfigElectrophys()
{
    // Set map of entries in configuration file to ap model types.
    this->ep_model_map_["bueno2008"] = ELECTRA::EpModelType::Bueno;
    this->ep_model_map_["ohara2011m"] = ELECTRA::EpModelType::OHara;
    this->ep_model_map_["gaur2021"] = ELECTRA::EpModelType::Gaur2021;
    this->ep_model_map_["gong2020"] = ELECTRA::EpModelType::Gong2020;
    this->ep_model_map_["gong2020m"] = ELECTRA::EpModelType::Gong2020m;
    this->ep_model_map_["paci2013v"] = ELECTRA::EpModelType::PaciVentri;
    this->ep_model_map_["tentusscher2006"] = ELECTRA::EpModelType::TenTusscher2006;
    this->ep_model_map_["courtemanche1998"] = ELECTRA::EpModelType::Courtemanche;
    this->ep_model_map_["grandi2011a"] = ELECTRA::EpModelType::GrandiAtri;
    this->ep_model_map_["maleckar2009"] = ELECTRA::EpModelType::Maleckar2009;
    this->ep_model_map_["maccannell2007"] = ELECTRA::EpModelType::MacCannell;
    this->ep_model_map_["stewart2009"] = ELECTRA::EpModelType::Stewart;

    // Set map of entries in configuration file to cell types.
    this->cell_type_map_["ventricular"] = ELECTRA::CellType::ventricular;
    this->cell_type_map_["endo"] = ELECTRA::CellType::endo;
    this->cell_type_map_["mid"] = ELECTRA::CellType::mid;
    this->cell_type_map_["epi"] = ELECTRA::CellType::epi;
    this->cell_type_map_["fibro"] = ELECTRA::CellType::fibro;
    this->cell_type_map_["atrial"] = ELECTRA::CellType::atrial;
    this->cell_type_map_["left atrial"] = ELECTRA::CellType::left_atrial;
    this->cell_type_map_["right atrial"] = ELECTRA::CellType::right_atrial;
    this->cell_type_map_["purkinje"] = ELECTRA::CellType::purkinje;
}


ConfigElectrophys::~ConfigElectrophys()
{}


void ConfigElectrophys::SetCellElectrophysiology(const Parser &parser, const std::unordered_map<std::string, IMP::NodeSet> &node_sets,
                                                 const std::string &body_type, int nodal_cells_num, std::vector<std::unique_ptr<ELECTRA::EpBasic>> &nodal_cells,
                                                 std::vector<ELECTRA::EpVaryingParams> &cell_varying_param_groups, std::ostream &stream) const
{
    // Check that body type is valid.
    if (body_type != "tissue" && body_type != "conduction system") {
        std::string error_msg = "Could not set cell electrophysiology. Body type must be: tissue or conduction system.";
        throw std::invalid_argument(ELECTRA::Logger::Error(error_msg));
    }

    // Reset nodal cells vector.
    nodal_cells.clear();
    nodal_cells.resize(nodal_cells_num);
    if (parser.HasAttribute(body_type+".electrophysiology.models number")) {
        // Set electrophysiology model at tissue cells from node sets.
        this->SetCellModelsFromNodeSets(parser, node_sets, body_type, nodal_cells_num, nodal_cells, cell_varying_param_groups, stream);
    } else {
        // Set electrophysiology model at tissue cells individually.
        this->SetCellModelsIndividually(parser, body_type, nodal_cells_num, nodal_cells, cell_varying_param_groups, stream);
    }

}


void ConfigElectrophys::SetCellModelsFromNodeSets(const Parser &parser, const std::unordered_map<std::string, IMP::NodeSet> &node_sets,
        const std::string &body_type, int nodal_cells_num, std::vector<std::unique_ptr<ELECTRA::EpBasic>> &nodal_cells,
        std::vector<ELECTRA::EpVaryingParams> &cell_varying_param_groups, std::ostream &stream) const
{
    // Get number of electrophysiology model nodesets.
    unsigned short models_num = parser.GetValue<unsigned short>(body_type+".electrophysiology.models number");

    // Assign electrophysiology model to the nodes of each nodeset.
    cell_varying_param_groups.clear();
    std::set<std::string> used_ep_models;
    std::vector<bool> initialized_nodes_flag(nodal_cells.size(), false);
    for (unsigned short i = 1; i <= models_num; ++i) {
        std::string ep_model_path = body_type+".electrophysiology.model-"+std::to_string(i);

        // Get electrophysiology model name.
        std::string ep_model_name = parser.GetValue<std::string>(ep_model_path+".model type");
        used_ep_models.insert(ep_model_name);
        std::transform(std::begin(ep_model_name), std::end(ep_model_name), std::begin(ep_model_name), ::tolower);
        if (this->ep_model_map_.find(ep_model_name) == this->ep_model_map_.end()) {
            std::string error_msg = "Electrophysiology model type: [" + ep_model_name + "] is not supported.";
            throw std::invalid_argument(ELECTRA::Logger::Error(error_msg));
        }

        // Get cell type.
        std::string cell_type = parser.GetValue<std::string>(ep_model_path+".cell type");
        std::transform(std::begin(cell_type), std::end(cell_type), std::begin(cell_type), ::tolower);
        if (this->cell_type_map_.find(cell_type) == this->cell_type_map_.end()) {
            std::string error_msg = "Electrophysiology cell type: [" + cell_type + "] is not supported.";
            throw std::invalid_argument(ELECTRA::Logger::Error(error_msg));
        }

        // Get nodeset name.
        std::string nset_name = parser.GetValue<std::string>(ep_model_path+".nodeset");
        if (node_sets.find(nset_name) == node_sets.end()) {
            std::string error_msg = "Could not set cell models from nodesets. The nodeset [" + nset_name + "] does not exist.";
            throw std::runtime_error(ELECTRA::Logger::Error(error_msg));
        }

        // Get manual cell initialization file if given.
        std::string manual_init_file = "";
        if (parser.HasAttribute(ep_model_path+".manual init file"))
            manual_init_file = parser.GetValue<std::string>(ep_model_path+".manual init file");

        // Set nodal cells of the specified node set.
        for (const auto &id : node_sets.at(nset_name).NodeIds()) {
            initialized_nodes_flag[id] = true;
            nodal_cells[id] = ELECTRA::EpFactory::Create(this->ep_model_map_.at(ep_model_name));
            nodal_cells[id]->Initialize(this->cell_type_map_.at(cell_type));
            if (manual_init_file != "") this->ManualCellInitialization(parser, manual_init_file, nodal_cells[id]);
        }

        // Get loading curve file if given.
        std::string load_curve_file = "";
        ELECTRA::EpVaryingParams ep_var_params;
        if (parser.HasAttribute(ep_model_path+".load curve file"))
            load_curve_file = parser.GetValue<std::string>(ep_model_path+".load curve file");

        // Read the load curves.
        if (load_curve_file != "") {
            this->SetEpVaryingParams(parser, load_curve_file, ep_model_name, ep_var_params);
            ep_var_params.SetCellIds(node_sets.at(nset_name).NodeIds());
            cell_varying_param_groups.emplace_back(ep_var_params);
        }
    }

    int initialized_counter = 0;
    for (const auto &initialized : initialized_nodes_flag) {
        if (initialized) initialized_counter++;
    }
    if (initialized_counter != nodal_cells_num) {
        std::string error_msg = "Could not set cell models from nodesets. The nodesets contain [" + std::to_string(initialized_counter) +
                                "] nodes but the model contains [" + std::to_string(nodal_cells_num) + "] nodes.";
        throw std::runtime_error(ELECTRA::Logger::Error(error_msg));
    }

    // Report used electrophysiology model types.
    stream << ELECTRA::Logger::Message("Used electrophysiology models\n");
    for (const auto &it : used_ep_models)
        stream << "                 - " + it + "\n";

}


void ConfigElectrophys::SetCellModelsIndividually(const Parser &parser, const std::string &body_type, int nodal_cells_num,
        std::vector<std::unique_ptr<ELECTRA::EpBasic>> &nodal_cells, std::vector<ELECTRA::EpVaryingParams> &cell_varying_param_groups, std::ostream &stream) const
{
    // Get file for manual nodal cell parameters initialization if available.
    std::string manual_init_file = "";
    bool with_manual_init = false;
    if (parser.HasAttribute(body_type+".electrophysiology.manual init file")) {
        manual_init_file = parser.GetValue<std::string>(body_type+".electrophysiology.manual init file");
        with_manual_init = true;
    }

    // Allow manual initialization only if ep model and cell type are in single value format.
    if (parser.IsArray(body_type+".electrophysiology.model type") || parser.IsArray(body_type+".electrophysiology.cell type")) {
        with_manual_init = false;
        stream << termcolor::yellow << ELECTRA::Logger::Warning("Manual cell initialization is not supported when either electrophysiology model or cell type is in array format. Default initialization will be used instead.") << termcolor::reset;
    }

    // Assign ep model types.
    std::set<std::string> used_ep_models;
    if (parser.IsSingleValue(body_type+".electrophysiology.model type")) {
        std::string ep_model_name = parser.GetValue<std::string>(body_type+".electrophysiology.model type");
        used_ep_models.insert(ep_model_name);
        std::transform(std::begin(ep_model_name), std::end(ep_model_name), std::begin(ep_model_name), ::tolower);

        // Same ep model for all nodal cells in the range.
        for (int id = 0; id != nodal_cells_num; ++id) {
            if (this->ep_model_map_.find(ep_model_name) != this->ep_model_map_.end()) {
                nodal_cells[id] = ELECTRA::EpFactory::Create(this->ep_model_map_.at(ep_model_name));
            } else {
                std::string error_msg = "Electrophysiology model type: " + ep_model_name + " is not supported.";
                throw std::runtime_error(ELECTRA::Logger::Error(error_msg));
            }
        }
    } else if (parser.IsArray(body_type+".electrophysiology.model type")) {
        // For the case that a different type of ep model is declared for each nodal cell.
        auto object = parser.GetObject(body_type+".electrophysiology.model type");
        auto ep_model_names = object.get<std::vector<std::string>>();

        if (ep_model_names.size() != static_cast<std::size_t>(nodal_cells_num)) {
            std::string error_msg = "Could not read electrophysiology model. Model type is not defined for all the nodal cells.";
            throw std::runtime_error(ELECTRA::Logger::Error(error_msg));
        }

        // Different ep model for each nodal cell in the range.
        int i = 0;
        for (int id = 0; id != nodal_cells_num; ++id) {
            used_ep_models.insert(ep_model_names[i]);
            std::transform(std::begin(ep_model_names[i]), std::end(ep_model_names[i]), std::begin(ep_model_names[i]), ::tolower);
            if (this->ep_model_map_.find(ep_model_names[i]) != this->ep_model_map_.end()) {
                nodal_cells[id] = ELECTRA::EpFactory::Create(this->ep_model_map_.at(ep_model_names[i]));
            } else {
                std::string error_msg = "Electrophysiology model type: " + ep_model_names[i] + " is not supported.";
                throw std::runtime_error(ELECTRA::Logger::Error(error_msg));
            }
            i++;
        }
    }

    // Assign cell types.
    if (parser.IsSingleValue(body_type+".electrophysiology.cell type")) {
        std::string cell_type_name = parser.GetValue<std::string>(body_type+".electrophysiology.cell type");
        std::transform(std::begin(cell_type_name), std::end(cell_type_name), std::begin(cell_type_name), ::tolower);

        // Same cell type for all nodal cells in the range.
        for (int id = 0; id != nodal_cells_num; ++id) {
            if (this->cell_type_map_.find(cell_type_name) != this->cell_type_map_.end()) {
                nodal_cells[id]->Initialize(this->cell_type_map_.at(cell_type_name));
                if (with_manual_init) {
                    this->ManualCellInitialization(parser, manual_init_file, nodal_cells[id]);
                }
            } else {
                std::string error_msg = "Electrophysiology cell type: " + cell_type_name + "is not supported.";
                throw std::runtime_error(ELECTRA::Logger::Error(error_msg));
            }
        }
    } else if (parser.IsArray(body_type+".electrophysiology.cell type")) {
        auto object = parser.GetObject(body_type+".electrophysiology.cell type");
        auto cell_type_names = object.get<std::vector<std::string>>();

        if (cell_type_names.size() != static_cast<std::size_t>(nodal_cells_num)) {
            std::string error_msg = "Could not read electrophysiology cell type. Cell type is not defined for all the nodal cells.";
            throw std::runtime_error(ELECTRA::Logger::Error(error_msg));
        }

        // Different cell type for each nodal cell in the range.
        int i = 0;
        for (int id = 0; id != nodal_cells_num; ++id) {
            std::transform(std::begin(cell_type_names[i]), std::end(cell_type_names[i]), std::begin(cell_type_names[i]), ::tolower);
            if (this->cell_type_map_.find(cell_type_names[i]) != this->cell_type_map_.end()) {
                nodal_cells[id]->Initialize(this->cell_type_map_.at(cell_type_names[i]));
            } else {
                std::string error_msg = "Electrophysiology cell type: " + cell_type_names[i] + " is not supported.";
                throw std::runtime_error(ELECTRA::Logger::Error(error_msg));
            }
            i++;
        }
    }

    // Get file for nodal cell varying parameters if available.
    std::string var_params_file = "";
    bool with_var_params = false;
    if (parser.HasAttribute(body_type+".electrophysiology.load curve file")) {
        var_params_file = parser.GetValue<std::string>(body_type+".electrophysiology.load curve file");
        with_var_params = true;
    }

    // Allow varying parameters initialization only if ep model and cell type are in single value format.
    if (parser.IsArray(body_type+".electrophysiology.model type") || parser.IsArray(body_type+".electrophysiology.cell type")) {
        with_var_params = false;
        stream << termcolor::yellow << ELECTRA::Logger::Warning("Varying parameters initialization is not supported when either electrophysiology model or cell type is in array format. Default initialization will be used instead.") << termcolor::reset;
    }

    if (with_var_params) {
        ELECTRA::EpVaryingParams ep_var_params;
        std::string ep_model_name = parser.GetValue<std::string>(body_type+".electrophysiology.model type");
        std::transform(std::begin(ep_model_name), std::end(ep_model_name), std::begin(ep_model_name), ::tolower);
        this->SetEpVaryingParams(parser, var_params_file, ep_model_name, ep_var_params);
        cell_varying_param_groups.clear();
        cell_varying_param_groups.emplace_back(ep_var_params);
    }

    // Report used electrophysiology model types.
    stream << ELECTRA::Logger::Message("Used electrophysiology models:\n");
    for (const auto &it : used_ep_models) {
        stream << "                 - " + it + "\n";
    }
}


void ConfigElectrophys::ManualCellInitialization(const Parser &parser, const std::string &init_file, std::unique_ptr<ELECTRA::EpBasic> &cell) const
{

    // Resolve manual initialization file path.
    std::string filename = parser.ResolvePath(init_file);

    // Open manual initialization file.
    std::ifstream init(filename, std::ios::in);
    if (!init.is_open()) {
        std::string error_msg = "Could not open file for manual electrophysiology model initialization. Check given file path: "+filename;
        throw std::invalid_argument(ELECTRA::Logger::Error(error_msg));
    }

    // Initialize line, key and val entries for file processing.
    std::stringstream ss;
    std::string key = "";
    std::string line = "";
    double val = 0.;

    // Initialize boolean to check the type of the data in process.
    bool is_var = false;
    bool is_prm = false;
    bool is_block = false;

    // Iterate over the file.
    while (std::getline(init, line)) {

        // Search for a data category.
        if (line.find("[") != std::string::npos) {
            // Make letters small to search for exact match.
            std::transform(line.begin(), line.end(), line.begin(), ::tolower);
            line = line.substr(line.find_first_of("["), line.find_first_of("]")+1);

            // Set data boolean corresponding to the current category.
            if (line == "[variables]") { is_var = true; is_prm = false; is_block = false; }
            else if (line == "[parameters]") { is_var = false; is_prm = true; is_block = false; }
            else if (line == "[current blocks]") { is_var = false; is_prm = false; is_block = true; }

        }

        // Try to set data if a category was found.
        if (is_var || is_prm || is_block) {

            // Check that the line in neither category header nor comment nor empty.
            if (line.find("[") == std::string::npos && line.find("#") == std::string::npos && line.length() > 1) {

                // Remove equal sign if necessary.
                std::replace(line.begin(), line.end(), '=', ' ');

                // Extract key and value of the data.
                ss << line;
                ss >> key >> val;
                // Reset stringstream for next iteration.
                ss.str(std::string()); ss.clear();

                // Check if key exist in the electrophysiology model's mapped data.
                if (!cell->HasDataEntry(key)) {
                    std::string error_msg = "Could not initialize electrophysiology model manually. Data with key: " + key + " does not belong to the model.";
                    throw std::out_of_range(ELECTRA::Logger::Error(error_msg));
                }

                //Set the data to the corresponding entry of the ap model of the cell.
                if (is_var) { cell->SetVar(cell->MappedDataEntry(key), val); }
                else if (is_prm) { cell->SetPrm(cell->MappedDataEntry(key), val); }
                else if (is_block) { cell->SetBlockCoeff(cell->MappedDataEntry(key), val); }
            }
        } // End of Try to set data if a category was found.

    } // End of Iterate over the file.

}


void ConfigElectrophys::SetEpVaryingParams(const Parser &parser, const std::string &load_curve_file, const std::string &ep_model_name,
        ELECTRA::EpVaryingParams &ep_var_params) const
{
    // Resolve manual initialization file path.
    std::string filename = load_curve_file;
    std::filesystem::path filepath(filename);
    std::string root_path_string{filepath.root_path()};
    if (!filepath.has_root_path() || root_path_string.size() == 1) {
        if (filename[0] == '/' || filename[0] == '\\') {
            filename = parser.ParentPath() + filename;
        } else {
            filename = parser.ParentPath() + "/" + filename;
        }
    }

    // Open loading curves file.
    std::ifstream file(filename, std::ios_base::in);
    if (!file.is_open()) {
        std::string error_msg = "Could not open file to read load curves for electrophysiology model. Check given path.";
        throw std::invalid_argument(ELECTRA::Logger::Error(error_msg));
    }

    // Containers to store the parameters names, time, and data.
    std::vector<std::string> params_names;
    std::vector<double> time;
    std::vector<ELECTRA::LoadCurve> params_load_curves;

    // Read load curves file.
    std::string line = "";
    std::stringstream ss;
    while (std::getline(file, line)) {
        // Get the time-varying parameters names.
        if (ELECTRA::ALGORITHM::ExistWord(line, "data")) {
            // Read the "#data:" part of the line.
            std::replace(std::begin(line), std::end(line), ':', ' ');
            ss << line;

            // Read the parameters names. First is always the time.
            std::string data = "";
            ss >> data; ss >> data; // Skip the part "Data time" from the header.
            while (ss >> data) { params_names.emplace_back(data); }
            ss.str(""); ss.clear();

            // Initialize the parameters data containers. First params_name refers to time.
            params_load_curves.resize(params_names.size());
        }

        // Get the load curves data.
        if (line.find("#") == std::string::npos && !ELECTRA::ALGORITHM::ExistWord(line, "d")) {
            double val = 0.;
            ss << line;

            // Store time entry.
            ss >> val;
            time.emplace_back(val);

            // Store parameters data.
            int id = 0;
            while (ss >> val) { params_load_curves[id++].AddData(val); }
            ss.str(""); ss.clear();
        }

    } // End of Read file.

    // Create a cell of the required type for data mapping.
    std::unique_ptr<ELECTRA::EpBasic> cell = ELECTRA::EpFactory::Create(this->ep_model_map_.at(ep_model_name));

    // Store the load curves in an unordered map.
    std::unordered_map<std::size_t, ELECTRA::LoadCurve> loads;
    for (std::size_t i = 0; i != params_load_curves.size(); ++i) {
        loads[cell->MappedDataEntry(params_names[i])] = params_load_curves[i];
    }

    // Set the electrophysiology varying parameters.
    ep_var_params.Clear();
    ep_var_params.SetTime(time);
    ep_var_params.SetLoads(loads);

}

} // end of namespace APP_ELECTRA