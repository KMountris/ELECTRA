/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */


#include "ELECTRA/engine/electrophysiology/maccannell.hpp"


namespace ELECTRA {


void Maccannell::SetDataMapping()
{
    using namespace McnVar;
    using namespace McnPrm;
    using namespace McnCur;

    // Set variables mapping.
    this->mapped_data_["v"] = static_cast<std::size_t>(v);
    this->mapped_data_["rf"] = static_cast<std::size_t>(rf);
    this->mapped_data_["sf"] = static_cast<std::size_t>(sf);
    this->mapped_data_["gKv"] = static_cast<std::size_t>(gKv);
    this->mapped_data_["gK1"] = static_cast<std::size_t>(gK1);
    this->mapped_data_["gbna"] = static_cast<std::size_t>(gbna);
    this->mapped_data_["ko"] = static_cast<std::size_t>(ko);
    this->mapped_data_["ki"] = static_cast<std::size_t>(ki);
    this->mapped_data_["nao"] = static_cast<std::size_t>(nao);
    this->mapped_data_["nai"] = static_cast<std::size_t>(nai);

    // Set parameters mapping.
    this->mapped_data_["R"] = static_cast<std::size_t>(R);
    this->mapped_data_["T"] = static_cast<std::size_t>(T);
    this->mapped_data_["Fdy"] = static_cast<std::size_t>(Fdy);
    this->mapped_data_["RTF"] = static_cast<std::size_t>(RTF);
    this->mapped_data_["iRTF"] = static_cast<std::size_t>(iRTF);
    this->mapped_data_["INaKmax"] = static_cast<std::size_t>(INaKmax);
    this->mapped_data_["KmKf"] = static_cast<std::size_t>(KmKf);
    this->mapped_data_["KmNaf"] = static_cast<std::size_t>(KmNaf);
    this->mapped_data_["Bf"] = static_cast<std::size_t>(Bf);
    this->mapped_data_["Vrev"] = static_cast<std::size_t>(Vrev);

    // Set currents mapping.
    this->mapped_data_["INaK"] = static_cast<std::size_t>(INaK);
    this->mapped_data_["IKv"] = static_cast<std::size_t>(IKv);
    this->mapped_data_["IbNa"] = static_cast<std::size_t>(IbNa);
    this->mapped_data_["IK1"] = static_cast<std::size_t>(IK1);

}


Maccannell::Maccannell()
{
    this->model_type_ = EpModelType::MacCannell;
    this->dt_stable_ = 0.02;
    this->var_.resize(11, 0.);
    this->prm_.resize(10, 0.);
    this->cur_.resize(5, 0.);
    this->block_coeff_.resize(4, 0.);

    // Set mapped data.
    this->SetDataMapping();
}


Maccannell::~Maccannell()
{}


void Maccannell::Initialize(CellType cell_type)
{
    using namespace McnVar;
    using namespace McnPrm;

    // Check required cell type.
    if (cell_type != CellType::fibro) {
        std::string error_str = Logger::Error("Could not initialize Maccannel ap model. Expected CellType::fibro.");
        throw std::invalid_argument(error_str);
    }

    //Initialize the model's data.
    this->var_.clear();          this->var_.resize(11, 0.);
    this->prm_.clear();          this->prm_.resize(10, 0.);
    this->cur_.clear();          this->cur_.resize(5, 0.);
    this->block_coeff_.clear();  this->block_coeff_.resize(4, 0.);

    // Set variables.
    this->var_[v]    = -49.6;
    this->var_[dvdt] = 0.;
    this->var_[rf]   = 0.063508346419;
    this->var_[sf]   = 0.976791239;
    this->var_[gKv]  = 0.25; 
    this->var_[gK1]  = 0.4822;
    this->var_[gbna] = 0.0095;
    this->var_[ko]   = 5.4;         // Value from Grandi 2010.
    this->var_[ki]   = 145.;        // Value from Grandi 2010.
    this->var_[nao]  = 140.;        // Value from Grandi 2010.
    this->var_[nai]  = 9.05;        // Value from Grandi 2010.

    // Set parameters.
    this->prm_[R]       = 8.314;
    this->prm_[T]       = 308.0;
    this->prm_[Fdy]     = 96.485;
    this->prm_[RTF]     = this->prm_[R] * this->prm_[T] / this->prm_[Fdy];
    this->prm_[iRTF]    = this->prm_[Fdy] / this->prm_[R] / this->prm_[T];
    this->prm_[INaKmax] = 2.002;
    this->prm_[KmKf]    = 1.;
    this->prm_[KmNaf]   = 11.;
    this->prm_[Bf]      = -200.;
    this->prm_[Vrev]    = -150.;

}


void Maccannell::Compute(double v_new, double dt, double stim_current)
{
    using namespace McnVar;
    using namespace McnPrm;
    using namespace McnCur;

    // Compute Nerst potentials.
    double ENa = this->prm_[RTF] * std::log(this->var_[nao] / this->var_[nai]);
    double EK  = -87.325488706291623;

    // Compute Sodium-potassium pump current.
    this->cur_[INaK] = (1.0-this->block_coeff_[INaK]) * ((this->prm_[INaKmax]*this->var_[ko] / (this->var_[ko] + this->prm_[KmKf])) * (std::pow(this->var_[nai], 1.5) / (std::pow(this->var_[nai], 1.5) + std::pow(this->prm_[KmNaf], 1.5)) ) * ((v_new - this->prm_[Vrev]) / (v_new - this->prm_[Bf])));
                        
    // Compute time and voltage dependent potassium current.
    this->cur_[IKv] = (1.0-this->block_coeff_[IKv]) * (this->var_[gKv]*this->var_[rf]*this->var_[sf] * (v_new - EK));

    // Compute background sodium current.
    this->cur_[IbNa] = (1.0-this->block_coeff_[IbNa]) * (this->var_[gbna] * (v_new - ENa));
    
    // Compute inward rectifying potassium current.
    double aK1 = 0.1 / (1. + std::exp(0.06 * (v_new - EK - 200.)));
    double bK1 = (3. * std::exp(0.0002*(v_new - EK + 100.)) + std::exp(0.1*(v_new - EK - 10.))) / (1. + std::exp(-0.5*(v_new - EK)));
    this->cur_[IK1] = (1.0-this->block_coeff_[IK1]) * ((this->var_[gK1]*aK1*(v_new - EK)) / (aK1 + bK1));

    // Total ionic current.
    this->cur_[McnCur::Iion] = this->cur_[INaK] + this->cur_[IKv] + this->cur_[IbNa] + this->cur_[IK1];

    // Set the new value to the membrante potential time derivative.
    this->var_[dvdt] = - (this->cur_[McnCur::Iion] - stim_current);

    // Update gate r variable.
    double rfinf = 1. / (1. + std::exp(-(v_new + 20.) / 11.));
    double taurf = 20.3 + 138. * std::exp(-std::pow((v_new + 20.) / 25.9, 2.));
    this->var_[rf] = ALGORITHM::RushLarsen(rfinf, this->var_[rf], dt, taurf);

    // Update gate s variable.
    double sfinf = 1. / (1. + std::exp((v_new + 23.) / 7.));
    double tausf = 1574. + 5268. * std::exp(-std::pow((v_new + 23.) / 22.7, 2.));
    this->var_[sf] = ALGORITHM::RushLarsen(sfinf, this->var_[sf], dt, tausf);
    
}


std::string Maccannell::PrintVariables() const
{
    using namespace McnVar;

    // Create output string stream to pass the parameters and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "v = " << this->var_[v]    << "\n";
    oss << "dvdt = " << this->var_[dvdt] << "\n";
    oss << "rf = " << this->var_[rf]   << "\n";
    oss << "sf = " << this->var_[sf]   << "\n";
    oss << "gKv = " << this->var_[gKv]  << "\n";
    oss << "gK1 = " << this->var_[gK1]  << "\n";
    oss << "gbna = " << this->var_[gbna] << "\n";
    oss << "ko = " << this->var_[ko]   << "\n";
    oss << "ki = " << this->var_[ki]   << "\n";
    oss << "nao = " << this->var_[nao]  << "\n";
    oss << "nai = " << this->var_[nai];
    return oss.str();
}


std::string Maccannell::PrintParameters() const
{
    using namespace McnPrm;

    // Create output string stream to pass the parameters and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "R = " << this->prm_[R] << "\n";
    oss << "T = " << this->prm_[T] << "\n";
    oss << "Fdy = " << this->prm_[Fdy] << "\n";
    oss << "RTF = " << this->prm_[RTF] << "\n";
    oss << "iRTF = " << this->prm_[iRTF] << "\n";
    oss << "INaKmax = " << this->prm_[INaKmax] << "\n";
    oss << "KmKf = " << this->prm_[KmKf] << "\n";
    oss << "KmNaf = " << this->prm_[KmNaf] << "\n";
    oss << "Bf = " << this->prm_[Bf] << "\n";
    oss << "Vrev = " << this->prm_[Vrev];
    return oss.str();
}


std::string Maccannell::PrintCurrents() const
{
    using namespace McnCur;

    // Create output string stream to pass the parameters and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "INaK = " << this->cur_[INaK] << "\n";
    oss << "IKv = " << this->cur_[IKv] << "\n";
    oss << "IbNa = " << this->cur_[IbNa] << "\n";
    oss << "IK1 = " << this->cur_[IK1] << "\n";
    oss << "Iion = " << this->cur_[McnCur::Iion];
    return oss.str();
}


std::string Maccannell::PrintBlockCoeffs() const
{
    using namespace McnCur;

    // Create output string stream to pass the parameters and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "INaK = " << this->block_coeff_[INaK] << "\n";
    oss << "IKv = " << this->block_coeff_[IKv] << "\n";
    oss << "IbNa = " << this->block_coeff_[IbNa] << "\n";
    oss << "IK1 = " << this->block_coeff_[IK1];
    return oss.str();
}


} // End of namespace ELECTRA