/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */
#ifndef ELECTRA_APPS_TOOLS_CONFIG_ELECTRIC_TPP_
#define ELECTRA_APPS_TOOLS_CONFIG_ELECTRIC_TPP_

#include "config_electric.hpp"

namespace APP_ELECTRA
{

template<short DIM, short CELL_NODES>
ConfigElectric<DIM, CELL_NODES>::ConfigElectric() : electric_mat_map_()
{
    this->electric_mat_map_["transversal"] = ELECTRA::ElectricType::transversal;
    this->electric_mat_map_["orthotropic"] = ELECTRA::ElectricType::orthotropic;
}


template<short DIM, short CELL_NODES>
ConfigElectric<DIM, CELL_NODES>::~ConfigElectric()
{}


template<short DIM, short CELL_NODES>
void ConfigElectric<DIM, CELL_NODES>::InitializeMaterial(const Parser &parser, std::shared_ptr<ELECTRA::ElectricBasic<DIM>> &electric_mat) const
{
    // Initialize electric material.
    std::string mat_name = parser.GetValue<std::string>("tissue.material.electric.type");
    std::transform(std::begin(mat_name), std::end(mat_name), std::begin(mat_name), ::tolower);
    electric_mat = ELECTRA::ElectricFactory<DIM>::Create(this->electric_mat_map_.at(mat_name));
}


template<short DIM, short CELL_NODES>
void ConfigElectric<DIM, CELL_NODES>::SetMaterialProperties(const Parser &parser, const std::unordered_map<std::string, IMP::NodeSet> &model_nodesets,
        const ELECTRA::MeasureUnits &units, std::shared_ptr<ELECTRA::ElectricBasic<DIM>> &electric_mat, std::ostream &stream) const
{
    // Assign transmembrane diffusivity values to the material. Used for monodomain model.
    if (parser.HasAttribute("tissue.material.electric.transmembrane diffusivity")) {
        std::vector<double> diff_values;
        this->AssignDiffusivity(parser, model_nodesets, "transmembrane diffusivity", units, electric_mat->NodesNum(), diff_values);
        electric_mat->SetTransDiffusivity(diff_values);
    }

    // Assign internal & extracellular diffusivity values to the material. Used for bidomain model.
    if (parser.HasAttribute("tissue.material.electric.intracellular diffusivity") ||
        parser.HasAttribute("tissue.material.electric.extracellular diffusivity")) {
        std::vector<double> diff_values;
        // Set intracellular diffusivity.
        this->AssignDiffusivity(parser, model_nodesets, "intracellular diffusivity", units, electric_mat->NodesNum(), diff_values);
        electric_mat->SetInterDiffusivity(diff_values);

        // Set extracellular diffusivity.
        this->AssignDiffusivity(parser, model_nodesets, "extracellular diffusivity", units, electric_mat->NodesNum(), diff_values);
        electric_mat->SetExterDiffusivity(diff_values);
    }

    // Assign transversal ratio values to the material.
    std::vector<double> ratio_values;
    this->AssignTransversalRatio(parser, model_nodesets, electric_mat->NodesNum(), ratio_values);
    electric_mat->SetTransversalRatio(ratio_values);

    // Assign capacitance to the material.
    std::vector<double> cap_values;
    this->AssignCapacitance(parser, model_nodesets, units, electric_mat->NodesNum(), cap_values);
    electric_mat->SetCapacitance(cap_values);

    // Assign fibers direction to the material.
    Eigen::MatrixXd fibers;
    this->AssignFibers(parser, electric_mat->NodesNum(), fibers);
    electric_mat->SetLongFibers(fibers);

    stream << ELECTRA::Logger::Message("Material parameters assigned successfully\n");

}


template<short DIM, short CELL_NODES>
void ConfigElectric<DIM, CELL_NODES>::ComputeDiffusionTensors(std::shared_ptr<ELECTRA::ElectricBasic<DIM>> &electric_mat, std::ostream &stream) const
{
    (DIM > 1 && CELL_NODES == 2) ? electric_mat->ComputeDiffusionTensors1D() : electric_mat->ComputeDiffusionTensors();
    stream << ELECTRA::Logger::Message("Material nodal diffusion tensors computed successfully.\n");
}


template<short DIM, short CELL_NODES>
void ConfigElectric<DIM, CELL_NODES>::AssignDiffusivity(const Parser &parser, const std::unordered_map<std::string, IMP::NodeSet> &model_nodesets,
        const std::string &diff_type, const ELECTRA::MeasureUnits &units, int values_num, std::vector<double> &diff_values) const
{
    if (parser.IsSingleValue("tissue.material.electric."+diff_type)) {
        diff_values.clear();
        double value = parser.GetValue<double>("tissue.material.electric."+diff_type);
        diff_values.assign(values_num, value);
    } else if (parser.IsArray("tissue.material.electric."+diff_type)) {
        auto object = parser.GetObject("tissue.material.electric."+diff_type);
        diff_values = object.get<std::vector<double>>();

        // Check if the number of the parsed diffusivity values is the expected.
        if (static_cast<int>(diff_values.size()) != values_num) {
            std::string error_str = "Could not assign "+diff_type+" to electric material. The material has [" + std::to_string(values_num) +
                                    "] and the parsed diffusivity values are [" + std::to_string(diff_values.size()) + "]";
            throw std::runtime_error(ELECTRA::Logger::Error(error_str));
        }
    } else if (parser.IsMultiArray("tissue.material.electric."+diff_type)) {
        this->ObtainValuesFromNodesets(parser, "tissue.material.electric."+diff_type, model_nodesets, values_num, diff_values);

    } else {
        std::string error_str = "Could not assign "+diff_type+" values to electric material. Diffusivity attribute has unexpected format.";
        throw std::invalid_argument(ELECTRA::Logger::Error(error_str));
    }

    // Get units of diffusivity.
    std::string nominator = parser.GetValue<std::string>("tissue.material.electric.diffusivity unit");
    std::string denominator = "";
    std::size_t pos = nominator.find('/');
    if (pos != std::string::npos) {
        denominator = nominator.substr(pos+1);
        nominator = nominator.substr(0, pos);
    } else {
        throw std::invalid_argument(ELECTRA::Logger::Error("Could not extract tissue diffusivity unit. Expected unit format: area / time"));
    }

    // Apply diffusivity unit.
    std::transform(std::begin(diff_values), std::end(diff_values), std::begin(diff_values),
                   std::bind(std::multiplies<double>(), std::placeholders::_1, units[nominator]/units[denominator]));

}


template<short DIM, short CELL_NODES>
void ConfigElectric<DIM, CELL_NODES>::AssignTransversalRatio(const Parser &parser, const std::unordered_map<std::string, IMP::NodeSet> &model_nodesets,
                                                             int values_num, std::vector<double> &ratio_values) const
{
    if (parser.IsSingleValue("tissue.material.electric.transversal ratio")) {
        ratio_values.clear();
        double value = parser.GetValue<double>("tissue.material.electric.transversal ratio");
        ratio_values.assign(values_num, value);
    } else if (parser.IsArray("tissue.material.electric.transversal ratio")) {
        auto object = parser.GetObject("tissue.material.electric.transversal ratio");
        ratio_values = object.get<std::vector<double>>();

        // Check if the number of the parsed transversal ratio values is the expected.
        if (static_cast<int>(ratio_values.size()) != values_num) {
            std::string error_str = "Could not assign transversal ratio values to electric material. The material has [" + std::to_string(values_num) +
                                    "] and the parsed transversal ratio values are [" + std::to_string(ratio_values.size()) + "]";
            throw std::runtime_error(ELECTRA::Logger::Error(error_str));
        }
    } else if (parser.IsMultiArray("tissue.material.electric.transversal ratio")) {
        this->ObtainValuesFromNodesets(parser, "tissue.material.electric.transversal ratio", model_nodesets, values_num, ratio_values);

    } else {
        std::string error_str = "Could not assign transversal ratio values to electric material. Transversal ratio attribute has unexpected format.";
        throw std::invalid_argument(ELECTRA::Logger::Error(error_str));
    }
}


template<short DIM, short CELL_NODES>
void ConfigElectric<DIM, CELL_NODES>::AssignCapacitance(const Parser &parser, const std::unordered_map<std::string, IMP::NodeSet> &model_nodesets,
                                                        const ELECTRA::MeasureUnits &units, int values_num, std::vector<double> &cap_values) const
{
    if (parser.IsSingleValue("tissue.material.electric.capacitance")) {
        cap_values.clear();
        double value = parser.GetValue<double>("tissue.material.electric.capacitance");
        cap_values.assign(values_num, value);
    } else if (parser.IsArray("tissue.material.electric.capacitance")) {
        auto object = parser.GetObject("tissue.material.electric.capacitance");
        cap_values = object.get<std::vector<double>>();

        // Check if the number of the parsed capacitance values is the expected.
        if (static_cast<int>(cap_values.size()) != values_num) {
            std::string error_str = "Could not assign capacitance to electric material. The material has [" + std::to_string(values_num) +
                                    "] and the parsed capacitance values are [" + std::to_string(cap_values.size()) + "]";
            throw std::runtime_error(ELECTRA::Logger::Error(error_str));
        }
    } else if (parser.IsMultiArray("tissue.material.electric.capacitance")) {
        this->ObtainValuesFromNodesets(parser, "tissue.material.electric.capacitance", model_nodesets, values_num, cap_values);

    } else {
        std::string error_str = "Could not assign capacitance values to electric material. Capacitance attribute has unexpected format.";
        throw std::invalid_argument(ELECTRA::Logger::Error(error_str));
    }

    // Get unit of capacitance.
    std::string cap_unit = parser.GetValue<std::string>("tissue.material.electric.capacitance unit");

    // Apply capacitance unit.
    std::transform(std::begin(cap_values), std::end(cap_values), std::begin(cap_values),
                   std::bind(std::multiplies<double>(), std::placeholders::_1, units[cap_unit]));

}


template<short DIM, short CELL_NODES>
void ConfigElectric<DIM, CELL_NODES>::AssignFibers(const Parser &parser, int values_num, Eigen::MatrixXd &fibers) const
{

    if (parser.IsSingleValue("tissue.material.electric.fibers") && DIM == 1) {
        // Assign single value to all the 1D fibers.
        double dir = parser.GetValue<double>("tissue.material.electric.fibers");
        fibers = dir * Eigen::MatrixXd::Ones(values_num, DIM);

    } else if(parser.IsArray("tissue.material.electric.fibers")) {
        // Assign a single direction vector to ND fibers.
        auto object = parser.GetObject("tissue.material.electric.fibers");
        auto dir = object.get<std::vector<double>>();
        if(dir.size() != DIM) {
            std::string error_str = "Could not extract fibers direction correctly. The provided fiber direction vector has not the same dimensions as the electric material.";
            throw std::runtime_error(ELECTRA::Logger::Error(error_str));
        }

        fibers = Eigen::MatrixXd::Zero(values_num, DIM);
        for (int i = 0; i != values_num; ++i) {
            for (short d = 0; d != DIM; ++d) {
                fibers.coeffRef(i,d) = dir[d];
            }
        }


    } else if(parser.IsMultiArray("tissue.material.electric.fibers")) {
        // Assign a unique direction vector to ND fibers.
        auto object = parser.GetObject("tissue.material.electric.fibers");

        fibers = Eigen::MatrixXd::Zero(values_num, DIM);
        int i = 0; short d = 0;
        for (const auto &dir : object) {
            for (const auto &val : dir) fibers.coeffRef(i,d++) = val;
            i++;
            d = 0;
        }

        if (i != values_num) {
            std::string error_str = "Could not extract fibers direction correctly. The provided fiber vectors are [" + std::to_string(i) +
                                    "] while the model has [" + std::to_string(values_num) + "] nodes.";
            throw std::runtime_error(ELECTRA::Logger::Error(error_str));
        }

    } else {
        std::string error_str = "Could not extract fibers direction correctly. Check the given format of the fibers.";
        throw std::runtime_error(ELECTRA::Logger::Error(error_str));
    }

}


template<short DIM, short CELL_NODES>
void ConfigElectric<DIM, CELL_NODES>::ObtainValuesFromNodesets(const Parser &parser, const std::string &attribute,
                                                               const std::unordered_map<std::string, IMP::NodeSet> &nodesets, int values_num, std::vector<double> &values) const
{
    values.clear();
    values.resize(values_num, 0.);
    int counter = 0;
    double value = 0.;
    auto object = parser.GetObject(attribute);
    if (object[0][0].is_string() && object[0][1].is_string()) {
        // Case where values are given as strings together with nodeset names.
        auto entries = object.get<std::vector<std::pair<std::string,std::string>>>();

        // Assign diffusivity values to the nodes of each nodeset.
        for (const auto &entry : entries) {
            value = std::stod(entry.first);
            for (const auto &id : nodesets.at(entry.second).NodeIds()) {
                values[id] = value;
                counter++;
            }
        }

        // Check that all nodes obtained a diffusivity value.
        if (counter != values_num) {
            std::string error_str = "Could not obtain values from node sets correctly. The total number of nodes is [" + std::to_string(values_num) +
                                    "] while the number of the defined nodes in the node sets is [" + std::to_string(counter) + "]";
            throw std::runtime_error(ELECTRA::Logger::Error(error_str));
        }

    } else if (object[0][0].is_number() && object[0][1].is_string()) {
        // Case where values are given as numbers together with nodeset names.
        auto entries = object.get<std::vector<std::pair<double,std::string>>>();

        // Assign diffusivity values to the nodes of each nodeset.
        for (const auto &entry : entries) {
            value = entry.first;
            for (const auto &id : nodesets.at(entry.second).NodeIds()) {
                values[id] = value;
                counter++;
            }
        }

        // Check that all nodes obtained a diffusivity value.
        if (counter != values_num) {
            std::string error_str = "Could not obtain values from node sets correctly. The total number of nodes is [" + std::to_string(values_num) +
                                    "] while the number of the defined nodes in the node sets is [" + std::to_string(counter) + "]";
            throw std::runtime_error(ELECTRA::Logger::Error(error_str));
        }
    }

}


} // end of namespace APP_ELECTRA

#endif //ELECTRA_APPS_TOOLS_CONFIG_ELECTRIC_TPP_