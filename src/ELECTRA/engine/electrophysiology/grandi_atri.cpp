/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

#include "ELECTRA/engine/electrophysiology/grandi_atri.hpp"


namespace ELECTRA {


void GrandiAtri::SetDataMapping() 
{
    using namespace GrdAtrVar;
    using namespace GrdAtrPrm;
    using namespace GrdAtrCur;

    // Set variables mapping.
    this->mapped_data_["v"] = static_cast<std::size_t>(v);
    this->mapped_data_["mo"] = static_cast<std::size_t>(mo);
    this->mapped_data_["ho"] = static_cast<std::size_t>(ho);
    this->mapped_data_["jo"] = static_cast<std::size_t>(jo);
    this->mapped_data_["d_o"] = static_cast<std::size_t>(d_o);
    this->mapped_data_["fo"] = static_cast<std::size_t>(fo);
    this->mapped_data_["fcaBjo"] = static_cast<std::size_t>(fcaBjo);
    this->mapped_data_["fcaBslo"] = static_cast<std::size_t>(fcaBslo);
    this->mapped_data_["xtoso"] = static_cast<std::size_t>(xtoso);
    this->mapped_data_["ytoso"] = static_cast<std::size_t>(ytoso);
    this->mapped_data_["xtofo"] = static_cast<std::size_t>(xtofo);
    this->mapped_data_["ytofo"] = static_cast<std::size_t>(ytofo);
    this->mapped_data_["xkro"] = static_cast<std::size_t>(xkro);
    this->mapped_data_["xkso"] = static_cast<std::size_t>(xkso);
    this->mapped_data_["RyRro"] = static_cast<std::size_t>(RyRro);
    this->mapped_data_["RyRoo"] = static_cast<std::size_t>(RyRoo);
    this->mapped_data_["RyRio"] = static_cast<std::size_t>(RyRio);
    this->mapped_data_["NaBjo"] = static_cast<std::size_t>(NaBjo);
    this->mapped_data_["NaBslo"] = static_cast<std::size_t>(NaBslo);
    this->mapped_data_["TnCLo"] = static_cast<std::size_t>(TnCLo);
    this->mapped_data_["TnCHco"] = static_cast<std::size_t>(TnCHco);
    this->mapped_data_["TnCHmo"] = static_cast<std::size_t>(TnCHmo);
    this->mapped_data_["CaMo"] = static_cast<std::size_t>(CaMo);
    this->mapped_data_["Myoco"] = static_cast<std::size_t>(Myoco);
    this->mapped_data_["Myomo"] = static_cast<std::size_t>(Myomo);
    this->mapped_data_["SRBo"] = static_cast<std::size_t>(SRBo);
    this->mapped_data_["SLLjo"] = static_cast<std::size_t>(SLLjo);
    this->mapped_data_["SLLslo"] = static_cast<std::size_t>(SLLslo);
    this->mapped_data_["SLHjo"] = static_cast<std::size_t>(SLHjo);
    this->mapped_data_["SLHslo"] = static_cast<std::size_t>(SLHslo);
    this->mapped_data_["Csqnbo"] = static_cast<std::size_t>(Csqnbo);

    // Set parameters mapping.
    this->mapped_data_["Ca_sro"] = static_cast<std::size_t>(Ca_sro);
    this->mapped_data_["Najo"] = static_cast<std::size_t>(Najo);
    this->mapped_data_["Naslo"] = static_cast<std::size_t>(Naslo);
    this->mapped_data_["Naio"] = static_cast<std::size_t>(Naio);
    this->mapped_data_["Kio"] = static_cast<std::size_t>(Kio);
    this->mapped_data_["Cajo"] = static_cast<std::size_t>(Cajo);
    this->mapped_data_["Caslo"] = static_cast<std::size_t>(Caslo);
    this->mapped_data_["Caio"] = static_cast<std::size_t>(Caio);
    this->mapped_data_["rtoso"] = static_cast<std::size_t>(rtoso);
    this->mapped_data_["n41"] = static_cast<std::size_t>(n41);
    this->mapped_data_["C1o"] = static_cast<std::size_t>(C1o);
    this->mapped_data_["C2o"] = static_cast<std::size_t>(C2o);
    this->mapped_data_["C3o"] = static_cast<std::size_t>(C3o);
    this->mapped_data_["C4o"] = static_cast<std::size_t>(C4o);
    this->mapped_data_["C5o"] = static_cast<std::size_t>(C5o);
    this->mapped_data_["C6o"] = static_cast<std::size_t>(C6o);
    this->mapped_data_["C7o"] = static_cast<std::size_t>(C7o);
    this->mapped_data_["C8o"] = static_cast<std::size_t>(C8o);
    this->mapped_data_["C9o"] = static_cast<std::size_t>(C9o);
    this->mapped_data_["C10o"] = static_cast<std::size_t>(C10o);
    this->mapped_data_["C11o"] = static_cast<std::size_t>(C11o);
    this->mapped_data_["C12o"] = static_cast<std::size_t>(C12o);
    this->mapped_data_["C13o"] = static_cast<std::size_t>(C13o);
    this->mapped_data_["C14o"] = static_cast<std::size_t>(C14o);
    this->mapped_data_["C15o"] = static_cast<std::size_t>(C15o);
    this->mapped_data_["O1o"] = static_cast<std::size_t>(O1o);
    this->mapped_data_["rkuro"] = static_cast<std::size_t>(rkuro);
    this->mapped_data_["skuro"] = static_cast<std::size_t>(skuro);
    this->mapped_data_["n60"] = static_cast<std::size_t>(n60);
    this->mapped_data_["n61"] = static_cast<std::size_t>(n61);
    this->mapped_data_["n62"] = static_cast<std::size_t>(n62);
    this->mapped_data_["ikcaoo"] = static_cast<std::size_t>(ikcaoo);
    this->mapped_data_["Ach"] = static_cast<std::size_t>(Ach);
    this->mapped_data_["ISO"] = static_cast<std::size_t>(ISO);
    this->mapped_data_["ISOdlf1muM"] = static_cast<std::size_t>(ISOdlf1muM);
    this->mapped_data_["ISOdlf1nM"] = static_cast<std::size_t>(ISOdlf1nM);
    this->mapped_data_["drugvector"] = static_cast<std::size_t>(drugvector);
    this->mapped_data_["AF"] = static_cast<std::size_t>(AF);
    this->mapped_data_["RA"] = static_cast<std::size_t>(RA);
    this->mapped_data_["R"] = static_cast<std::size_t>(R);
    this->mapped_data_["Fdy"] = static_cast<std::size_t>(Fdy);
    this->mapped_data_["T"] = static_cast<std::size_t>(T);
    this->mapped_data_["l"] = static_cast<std::size_t>(l);
    this->mapped_data_["rad"] = static_cast<std::size_t>(rad);
    this->mapped_data_["Cmem"] = static_cast<std::size_t>(Cmem);

    // Set currents mapping.
    this->mapped_data_["INa"] = static_cast<std::size_t>(INa);
    this->mapped_data_["INaL"] = static_cast<std::size_t>(INaL);
    this->mapped_data_["INab"] = static_cast<std::size_t>(INab);
    this->mapped_data_["INaK"] = static_cast<std::size_t>(INaK);
    this->mapped_data_["IKr"] = static_cast<std::size_t>(IKr);
    this->mapped_data_["IKs"] = static_cast<std::size_t>(IKs);
    this->mapped_data_["IKp"] = static_cast<std::size_t>(IKp);
    this->mapped_data_["IKach"] = static_cast<std::size_t>(IKach);
    this->mapped_data_["IKCa"] = static_cast<std::size_t>(IKCa);
    this->mapped_data_["Ito"] = static_cast<std::size_t>(Ito);
    this->mapped_data_["IKur"] = static_cast<std::size_t>(IKur);
    this->mapped_data_["IKi"] = static_cast<std::size_t>(IKi);
    this->mapped_data_["IClb"] = static_cast<std::size_t>(IClb);
    this->mapped_data_["ICa"] = static_cast<std::size_t>(ICa);
    this->mapped_data_["ICaK"] = static_cast<std::size_t>(ICaK);
    this->mapped_data_["ICaNa"] = static_cast<std::size_t>(ICaNa);
    this->mapped_data_["ICatot"] = static_cast<std::size_t>(ICatot);
    this->mapped_data_["Incx"] = static_cast<std::size_t>(Incx);
    this->mapped_data_["Ipca"] = static_cast<std::size_t>(Ipca);
    this->mapped_data_["ICab"] = static_cast<std::size_t>(ICab);
    
}


GrandiAtri::GrandiAtri()
{
    //Initialize the model's data.
    this->model_type_ = EpModelType::GrandiAtri;
    this->dt_stable_ = 0.005;
    this->var_.resize(32, 0.);
    this->prm_.resize(45, 0.);
    this->cur_.resize(21, 0.);
    this->block_coeff_.resize(20, 0.);

    // Set mapped data.
    this->SetDataMapping();
}


GrandiAtri::~GrandiAtri()
{}


void GrandiAtri::Initialize(CellType cell_type)
{
    using namespace GrdAtrVar;
    using namespace GrdAtrPrm;

    if (cell_type != CellType::atrial) {
        std::string error_str = "Could not initialize the Grandi atrial ap model. Expected cell type: ELECTRA::CellType::atrial";
        throw std::invalid_argument(Logger::Error(error_str));
    }

    //Initialize the model's data.
    this->var_.clear();           this->var_.resize(32, 0.);
    this->prm_.clear();           this->prm_.resize(45, 0.);
    this->cur_.clear();           this->cur_.resize(21, 0.);
    this->block_coeff_.clear();   this->block_coeff_.resize(20, 0.);

    // Set the model variables.
    this->var_[v]       = -73.749725;
    this->var_[dvdt]    = 0.;
    this->var_[mo]      = 0.017812827; 
    this->var_[ho]      = 0.32819992;
    this->var_[jo]      = 0.28812384;
    this->var_[d_o]     = 0.000020563217;   
    this->var_[fo]      = 0.99847596;
    this->var_[fcaBjo]  = 0.043146994;  
    this->var_[fcaBslo] = 0.031431224;
    this->var_[xtoso]   = 0.004051574;
    this->var_[ytoso]   = 0.9945511;
    this->var_[xtofo]   = 0.0013399836;
    this->var_[ytofo]   = 0.94734859;
    this->var_[xkro]    = 0.0012226996;
    this->var_[xkso]    = 0.0073281652;
    this->var_[RyRro]   = 0.80768195;
    this->var_[RyRoo]   = 0.0000018335839;
    this->var_[RyRio]   = 0.00000043659156;
    this->var_[NaBjo]   = 3.612005;
    this->var_[NaBslo]  = 0.7881974;
    this->var_[TnCLo]   = 0.017742203;
    this->var_[TnCHco]  = 0.1271877;
    this->var_[TnCHmo]  = 0.0060177108;
    this->var_[CaMo]    = 0.00067501214;
    this->var_[Myoco]   = 0.0037309625;
    this->var_[Myomo]   = 0.13575808;
    this->var_[SRBo]    = 0.0043194537;
    this->var_[SLLjo]   = 0.013057019;
    this->var_[SLLslo]  = 0.020605241;
    this->var_[SLHjo]   = 0.10171044;
    this->var_[SLHslo]  = 0.18649672;
    this->var_[Csqnbo]  = 1.1204378;

    // Set the model parameters.
    this->prm_[Ca_sro]     = 0.49222997;
    this->prm_[Najo]       = 9.1414252;
    this->prm_[Naslo]      = 9.1419134;
    this->prm_[Naio]       = 9.1420919;
    this->prm_[Kio]        = 120.;
    this->prm_[Cajo]       = 0.000313349;
    this->prm_[Caslo]      = 0.00022415727;
    this->prm_[Caio]       = 0.00020247969;
    this->prm_[rtoso]      = 0.9946;
    this->prm_[n41]        = 1.;
    this->prm_[C1o]        = 0.0015;
    this->prm_[C2o]        = 0.0244;
    this->prm_[C3o]        = 0.1494;
    this->prm_[C4o]        = 0.4071;
    this->prm_[C5o]        = 0.4161;
    this->prm_[C6o]        = 0.;
    this->prm_[C7o]        = 0.0001;
    this->prm_[C8o]        = 0.0006;
    this->prm_[C9o]        = 0.0008;
    this->prm_[C10o]       = 0.;
    this->prm_[C11o]       = 0.;
    this->prm_[C12o]       = 0.;
    this->prm_[C13o]       = 0.;
    this->prm_[C14o]       = 0.;
    this->prm_[C15o]       = 0.;
    this->prm_[O1o]        = 0.;
    this->prm_[rkuro]      = 0.00037902416;
    this->prm_[skuro]      = 0.96508959;
    this->prm_[n60]        = 0.0096934451;
    this->prm_[n61]        = 0.044603201;
    this->prm_[n62]        = 0.;
    this->prm_[ikcaoo]     = 0.0065214488;
    this->prm_[Ach]        = 0.;
    this->prm_[ISO]        = 0.;
    this->prm_[ISOdlf1muM] = 0.;
    this->prm_[ISOdlf1nM]  = 0.;
    this->prm_[drugvector] = 1.;
    this->prm_[AF]         = 0.;
    this->prm_[RA]         = 0.;
    this->prm_[R]          = 8314.;
    this->prm_[Fdy]        = 96485.;
    this->prm_[T]          = 310.;
    this->prm_[l]          = 100.;
    this->prm_[rad]        = 10.25;
    this->prm_[Cmem]       = 1.1e-10;

}


void GrandiAtri::Compute(double v_new, double dt, double stim_current)
{
    using namespace GrdAtrVar;
    using namespace GrdAtrPrm;
    using namespace GrdAtrCur;

    const double pi = 3.141592653589793;

    // Compute constants.
    double FoRT = this->prm_[Fdy] / this->prm_[R] / this->prm_[T];
    double Qpow = (this->prm_[T] - 310.) / 10.;
    double v_cell = pi * this->prm_[rad]*this->prm_[rad] * this->prm_[l] * 1e-15;

    // Compute cell geometry constants.
    double Vmyo = 0.65 * v_cell; 
    double Vsr = 0.035 * v_cell; 
    double Vsl = 0.02 * v_cell; 
    double Vjunc = 0.000539 * v_cell; 
    
    // double SAjunc = 20150. * 3.14 * 2. * junctionLength*junctionRadius;     // [um^2]
    double J_ca_juncsl = 8.24130542278e-13;                                 // [L/msec]
    double J_ca_slmyo = 3.72425607985e-12;                                  // [L/msec]
    double J_na_juncsl = 1.83127823220608e-14;                              // [L/msec]
    double J_na_slmyo = 1.63862792221979e-12;                               // [L/msec]

    // Compute fractional currents in compartments.
    double Fjunc = 0.11;   
    double Fsl = 1. - Fjunc;
    double Fjunc_CaL = 0.9; 
    double Fsl_CaL = 1. - Fjunc_CaL;

    // Compute fixed ion concentrations.
    double Cli = 15.;   // Intracellular Cl  [mM]
    double Clo = 150.;  // Extracellular Cl  [mM]
    double Ko  = 5.4;  // Extracellular K   [mM]
    double Nao = 140.;  // Extracellular Na  [mM]
    double Cao = 1.8;  // Extracellular Ca  [mM]
    double Mgi = 1.;    // Intracellular Mg  [mM]

    // Compute Nernst potentials.
    double ena_junc = (1./FoRT) * std::log(Nao/this->prm_[Najo]);    // [mV]
    double ena_sl = (1./FoRT) * std::log(Nao/this->prm_[Naslo]);     // [mV]
    double ek = (1./FoRT) * std::log(Ko/this->prm_[Kio]);	         // [mV]
    double eca_junc = (1./FoRT/2.) * std::log(Cao/this->prm_[Cajo]); // [mV]
    double eca_sl = (1./FoRT/2.) * std::log(Cao/this->prm_[Caslo]);  // [mV]
    double ecl = (1./FoRT) * std::log(Cli/Clo);                      // [mV]

    // Compute Na transport parameters.
    double GNa = 23. * (1. - 0.1*this->prm_[AF]);       // [mS/uF]
    double GNaB = 0.000597;                             // [mS/uF] 
    double IbarNaK = 1.26;                              // [uA/uF]
    double KmNaip = 11. * (1. - 0.25*this->prm_[ISO]);  // [mM]
    double KmKo =1.5;                                   // [mM]

    // Compute K current parameters.
    double pNaK = 0.01833;      
    double gkp = 0.002;

    // Compute Cl current parameters.
    double GClCa = 0.0548;   // [mS/uF]
    double GClB  = 0.009;    // [mS/uF]
    double KdClCa = 0.1;     // [mM]

    // Compute I_Ca parameters.
    double pNa = (1. + 0.5*this->prm_[ISO]) * (1. + 1.5*this->prm_[ISOdlf1muM]) * (1. + 0.135 * this->prm_[ISOdlf1nM]) * (1. - 0.5*this->prm_[AF]) * 0.75e-8;  // [cm/sec]
    double pCa = (1. + 0.5*this->prm_[ISO]) * (1. + 1.5*this->prm_[ISOdlf1muM]) * (1. + 0.135 * this->prm_[ISOdlf1nM]) * (1. - 0.5*this->prm_[AF]) * 0.00027;  // [cm/sec]
    double pK = (1. + 0.5*this->prm_[ISO]) * (1. + 1.5*this->prm_[ISOdlf1muM]) * (1. + 0.135*this->prm_[ISOdlf1nM]) * (1. - 0.5*this->prm_[AF]) * 1.35e-7;     // [cm/sec]
    double Q10CaL = 1.8;       

    // Compute Ca transport parameters.
    double IbarNCX = (1. + 0.4*this->prm_[AF]) * 3.15;  // [uA/uF]
    double KmCai = 0.00359;                             // [mM]
    double KmCao = 1.3;                                 // [mM]
    double KmNai = 12.29;                               // [mM]
    double KmNao = 87.5;                                // [mM]
    double ksat = 0.27;          
    double nu = 0.35;
    double Kdact = 0.000384;                            // [mM]
    double Q10NCX = 1.57;
    double IbarSLCaP =  0.0471;
    double KmPCa = 0.0005;                              // [mM] 
    double GCaB = 0.00060643;                           // [uA/uF]
    double Q10SLCaP = 2.35;     

    // Compute SR flux parameters
    double Q10SRCaP = 2.6;
    double Vmax_SRCaP = 0.0053114;                          // [mM/msec]
    double Kmf = (2.5-1.25*this->prm_[ISO]) * 0.000246;     // [mM]
    double Kmr = 1.7;                                       // [mM]
    double hillSRCaP = 1.787;                               // [mM]
    double ks = 25.;                                        // [1/ms]      

    double koCa = 10. + 20.*this->prm_[AF] + 10.*this->prm_[ISO]*(1.-this->prm_[AF]);   // [mM^-2 1/ms]
    
    double kom = 0.06;      // [1/ms]     
    double kiCa = 0.5;      // [1/mM/ms]
    double kim = 0.005;     // [1/ms]
    double ec50SR = 0.45;   // [mM]
    
    // Buffering parameters
    double Bmax_Naj = 7.561;        // [mM]
    double Bmax_Nasl = 1.65;        // [mM]
    double koff_na = 0.001;         // [1/ms]
    double kon_na = 0.0001;         // [1/mM/ms]
    double Bmax_TnClow = 0.07;      // [mM]

    double koff_tncl = (1.+ 0.5*this->prm_[ISO])*0.0196;    // [1/ms] 
    
    double kon_tncl = 32.7;             // [1/mM/ms]
    double Bmax_TnChigh = 0.14;         // [mM]
    double koff_tnchca = 0.000032;      // [1/ms] 
    double kon_tnchca = 2.37;           // [1/mM/ms]
    double koff_tnchmg = 0.00333;       // [1/ms] 
    double kon_tnchmg = 0.003;          // [1/mM/ms]
    double Bmax_CaM = 0.024;            // [mM] **? about setting to 0 in c-code**   % CaM buffering
    double koff_cam = 0.238;            // [1/ms] 
    double kon_cam = 34.;               // [1/mM/ms]
    double Bmax_myosin = 0.14;          // Myosin buffering [mM]
    double koff_myoca = 0.00046;        // [1/ms]
    double kon_myoca = 13.8;            // [1/mM/ms]
    double koff_myomg = 0.000057;       // [1/ms]
    double kon_myomg = 0.0157;          // [1/mM/ms]
    double Bmax_SR = 19. * 0.0009;      // [mM] (Bers text says 47e-3) 19e-3
    double koff_sr = 0.060;             // [1/ms]
    double kon_sr = 100.;               // [1/mM/ms]
    
    double Bmax_SLlowsl = 0.0374 * Vmyo / Vsl;          // SL buffering [mM]
    double Bmax_SLlowj = 0.0046  * Vmyo/Vjunc * 0.1;    // low junction reduction factor [mM]
    double koff_sll = 1.3;                              // [1/ms]
    double kon_sll = 100.;                              // [1/mM/ms]
    double Bmax_SLhighsl = 0.0134 * Vmyo / Vsl;         // [mM] 
    double Bmax_SLhighj = 0.00165 * Vmyo / Vjunc * 0.1; // high junction reduction factor [mM]
    double koff_slh = 0.03;                             // [1/ms]
    double kon_slh = 100.;                              // [1/mM/ms]
    double Bmax_Csqn = 0.14 * Vmyo / Vsr;               // Csqn buffering [mM]
    double koff_csqn = 65.;                             // [1/ms] 
    double kon_csqn = 100.;                             // [1/mM/ms] 

    // Solve gate variable mo with Rash Larsen.
    double mss = 1. / std::pow(1. + std::exp(-(56.86 + v_new) / 9.03), 2.);
    double taum = 0.1292 * std::exp(-std::pow((v_new + 45.79)/15.54, 2.)) + 0.06487 * std::exp(-std::pow((v_new - 4.823)/51.12, 2.));
    this->var_[mo] = ALGORITHM::RushLarsen(mss, this->var_[mo], dt, taum);                     

    // Solve gate variable ho with Rash Larsen.
    double ah = 0.;
    double bh = 0.77 / (0.13 * (1.+std::exp(-(v_new + 10.66) / 11.1)));
    if (v_new < -40.) { 
        ah = 0.057 * std::exp(-(v_new + 80.) / 6.8);
        bh = 2.7*std::exp(0.079*v_new) + 310000.*std::exp(0.3485*v_new);
    }       
    double tauh = 1. / (ah + bh); 
    double hss = 1. / std::pow(1. + std::exp((v_new+71.55) / 7.43 ), 2.);
    this->var_[ho] = ALGORITHM::RushLarsen(hss, this->var_[ho], dt, tauh);
    
    // Solve gate variable jo with Rash Larsen.
    double aj = 0.;
    double bj = (0.6*std::exp(0.057*v_new)) / (1. + std::exp(-0.1*(v_new+32.))); 
    if (v_new < - 40.) { 
        aj = ((-25428*std::exp(0.2444*v_new) - 0.000006948*std::exp(-0.04391*v_new)) * (v_new+37.78)) / (1. + std::exp(0.311 * (v_new+79.23))); 
        bj = (0.02424*std::exp(-0.01052*v_new)) / (1.+std::exp(-0.1378 * (v_new+40.14)));
    }    
    double tauj = 1. / (aj + bj);
    double jss = 1. / std::pow(1.+std::exp((v_new+71.55) / 7.43), 2.);
    this->var_[jo] = ALGORITHM::RushLarsen(jss, this->var_[jo], dt, tauj);       

    // Compute sodium current.
    double INa_junc = Fjunc * GNa * this->var_[mo]*this->var_[mo]*this->var_[mo] * this->var_[ho] * this->var_[jo] * (v_new-ena_junc);
    double INa_sl = Fsl * GNa * this->var_[mo]*this->var_[mo]*this->var_[mo] * this->var_[ho] * this->var_[jo] * (v_new-ena_sl);
    this->cur_[INa] = (1.0-this->block_coeff_[INa]) * (INa_junc + INa_sl);

    // Solving gating variable n61 with Rush Larsen.
    double GNaL = 0.0025*this->prm_[AF];
    double aml = 0.32*(v_new+47.13) / (1.-std::exp(-0.1*(v_new+47.13)));
    double bml = 0.08*std::exp(-v_new/11.);
    double hlinf = 1. / (1.+std::exp((v_new+91.)/6.1));
    double tauhl = 600.;
    this->prm_[n61] = ALGORITHM::RushLarsen(hlinf, this->prm_[n61], dt, tauhl);

    double dn60 = aml * (1.0-this->prm_[n60]) - bml*this->prm_[n60];
    this->prm_[n60] = ALGORITHM::ForwardEuler(this->prm_[n60], dt, dn60);

    // Compute late sodium current.
    double INaL_junc = Fjunc * GNaL * this->prm_[n60]*this->prm_[n60]*this->prm_[n60] * this->prm_[n61] * (v_new-ena_junc);
    double INaL_sl = Fsl * GNaL * this->prm_[n60]*this->prm_[n60]*this->prm_[n60] * this->prm_[n61] * (v_new-ena_sl);
    this->cur_[INaL] = (1.0-this->block_coeff_[INaL]) * (INaL_junc + INaL_sl);

    this->prm_[n62] = this->cur_[INaL];

    // Compute sodium background current
    double INab_junc = Fjunc * GNaB * (v_new-ena_junc);
    double INab_sl = Fsl * GNaB * (v_new-ena_sl);
    this->cur_[INab] = (1.0-this->block_coeff_[INab]) * (INab_junc + INab_sl);

    // Compute sodium/potassium pump current.
    double sigma = (std::exp(Nao/67.3)-1.) / 7.;
    double fnak = 1. / (1. + 0.1245*std::exp(-0.1*v_new*FoRT) + 0.0365*sigma*std::exp(-v_new*FoRT));
    double INaK_junc = Fjunc * IbarNaK * fnak * Ko / (1.+std::pow(KmNaip/this->prm_[Najo], 4.)) / (Ko+KmKo);
    double INaK_sl = Fsl * IbarNaK * fnak * Ko / (1.+std::pow(KmNaip/this->prm_[Naslo], 4.)) / (Ko+KmKo);
    this->cur_[INaK] = (1.0-this->block_coeff_[INaK]) * (INaK_junc + INaK_sl);

    // Solving xkro gating variable with Rush Larsen.
    double gkr = 0.035*std::sqrt(Ko/5.4);
    double xrss = 1. / (1.+std::exp(-(v_new+10.) / 5.));
    double tauxr = 550. / (1.+std::exp((-22.-v_new) / 9.)) * 6. / (1.+std::exp((v_new+11.)/9.)) + 230. / (1.+std::exp((v_new+40.)/20.));
    this->var_[xkro] = ALGORITHM::RushLarsen(xrss, this->var_[xkro], dt, tauxr);

    // Compute rapidly activating potassium current.
    double rkr = 1. / (1.+std::exp((v_new+74.) / 24.));
    this->cur_[IKr] = (1.0-this->block_coeff_[IKr]) * (gkr * this->var_[xkro] * rkr * (v_new-ek));

    // Solving xkso gating variable with Rush Larsen.
    double gks_junc = (1. + this->prm_[AF] + 2.*this->prm_[ISO]) * (1. + 0.58*this->prm_[ISOdlf1muM]) * (1. + 0.3*this->prm_[ISOdlf1nM]) * 0.0035;
    double gks_sl = (1. + this->prm_[AF] + 2.*this->prm_[ISO]) * (1. + 0.58*this->prm_[ISOdlf1muM]) * (1. + 0.3*this->prm_[ISOdlf1nM]) * 0.0035;
    double xsss = 1. / (1. + std::exp(-(v_new + 40.*this->prm_[ISO] + 3.8) / 14.25));
    double tauxs = 990.1 / (1. + std::exp(-(v_new + 40.*this->prm_[ISO] + 2.436) / 14.12)); 
    this->var_[xkso] = ALGORITHM::RushLarsen(xsss, this->var_[xkso], dt, tauxs);
    
    // Compute slowly activating potassium current.
    double eks = (1./FoRT) * std::log((Ko + pNaK*Nao) / (this->prm_[Kio] + pNaK*this->prm_[Naio]));
    double IKs_junc = Fjunc * gks_junc * this->var_[xkso]*this->var_[xkso] * (v_new-eks);
    double IKs_sl = Fsl * gks_sl * this->var_[xkso]*this->var_[xkso] * (v_new-eks);                                                                                                                                   
    this->cur_[IKs] = (1.0-this->block_coeff_[IKs]) * (IKs_junc + IKs_sl);

    // Compute plateau potassium current.
    double kp_kp = 1. / (1.+std::exp(7.488 - v_new/5.98));
    double IKp_junc = Fjunc * gkp * kp_kp * (v_new-ek);
    double IKp_sl = Fsl * gkp * kp_kp * (v_new-ek);
    this->cur_[IKp] = (1.0-this->block_coeff_[IKp]) * (IKp_junc + IKp_sl);

    // Compute acetylcholine potassium current.
    this->cur_[IKach] = (1.0-this->block_coeff_[IKach]) * (1. / (1. + std::pow(0.03/this->prm_[Ach], 2.1)) * (v_new-ek) * (0.08 + 0.04 / (1. + std::exp((v_new + 91.) / 12.))));
    
    //Compute potassium-calcium current as in Krizanek, Echebarria & Penaranda.
    double tau = 3.;
    double gKCa = 6.*this->prm_[drugvector];
    double numzss = (this->prm_[Caio]/0.0025)*(this->prm_[Caio]/0.0025); 
    double denomzss = 1. + numzss;
    double yss = numzss / denomzss;
    this->prm_[ikcaoo] = ALGORITHM::RushLarsen(yss, this->prm_[ikcaoo], dt, tau);
    this->cur_[IKCa] = (1.0-this->block_coeff_[IKCa]) * (gKCa * this->prm_[ikcaoo]*this->prm_[ikcaoo] * (v_new-ek));
    
    // Compute transient outward potassium current (fast component) with modification for human myocytes.    
    double GtoFast = (1. - 0.7*this->prm_[AF]) * (1. - 0.38*this->prm_[ISOdlf1muM]) * (1. - 0.106*this->prm_[ISOdlf1nM]) * 0.165; // [nS/pF]

    // Equations for activation. 
    double xtoss = 1. / (1. + std::exp(-(v_new+1.) / 11.));
    double tauxtof = 3.5 * std::exp(-(v_new/30.)*(v_new/30.)) + 1.5;
    this->var_[xtofo] = ALGORITHM::RushLarsen(xtoss, this->var_[xtofo], dt, tauxtof);
    
    // Equations for inactivation.
    double ytoss = 1. / (1. + std::exp((v_new+40.5) / 11.5));
    double tauytof = 25.635 * std::exp(-std::pow((v_new+52.45)/15.8827, 2.)) + 24.14;    
    this->var_[ytofo] = ALGORITHM::RushLarsen(ytoss, this->var_[ytofo], dt, tauytof);

    this->cur_[Ito] = (1.0-this->block_coeff_[Ito]) * (GtoFast * this->var_[xtofo] * this->var_[ytofo] * (v_new-ek));

    // Compute ultra rapid delayed rectifier outward potassium current. Based on Maleckar et al. 2009 - EG
    // Equations for activation.
    double Gkur = (1. - 0.5*this->prm_[AF]) * (1. + 2.*this->prm_[ISO]) * 0.045*(1. + 0.2*this->prm_[RA]);  // [nS/pF]
    double xkurss = 1. / (1. + std::exp((v_new+6.) / -8.6));
    double tauxkur = 9. / (1. + std::exp((v_new+5.)/12.)) + 0.5;
    this->prm_[rkuro] = ALGORITHM::RushLarsen(xkurss, this->prm_[rkuro], dt, tauxkur);

    // Equations for inactivation.
    double ykurss = 1. / (1. + std::exp((v_new+7.5) / 10.));
    double tauykur = 590. / (1. + std::exp((v_new+60.) / 10.)) + 3050.;
    this->prm_[skuro] = ALGORITHM::RushLarsen(ykurss, this->prm_[skuro], dt, tauykur);

    this->cur_[IKur] = (1.0-this->block_coeff_[IKur]) * (Gkur * this->prm_[rkuro] * this->prm_[skuro] * (v_new-ek));

    // Compute time-independent potassium current.
    double aki = 1.02 / (1. + std::exp(0.2385*(v_new-ek-59.215)));
    double bki = (0.49124 * std::exp(0.08032*(v_new+5.476-ek)) + std::exp(0.06175*(v_new-ek-594.31))) / (1. + std::exp(-0.5143*(v_new-ek+4.753)));
    double kiss = aki / (aki+bki);
    
    // Multiplying IK1 by 0.15 to scale it to single cell isolated atrial cell resting potential
    this->cur_[IKi] = (1.0-this->block_coeff_[IKi]) * (1.+this->prm_[AF]) * 0.0525 * std::sqrt(Ko/5.4) * kiss * (v_new-ek);

    // Compute calcium-activated Cl Current and Cl background current.
    double IClCa_junc = Fjunc * GClCa / (1. + KdClCa/this->prm_[Cajo]) * (v_new-ecl);
    double IClCa_sl = Fsl * GClCa / (1. + KdClCa/this->prm_[Caslo]) * (v_new-ecl);
    double IClCa = IClCa_junc + IClCa_sl;
    this->cur_[IClb] = (1.0-this->block_coeff_[IClb]) * (GClB * (v_new-ecl));

    // Compute L-type calcium current.
    double dss = 1. / (1. + std::exp(-(v_new+3.*this->prm_[ISO] + 9.) / 6.)); // %in Maleckar v1/2=-9 S=6 (mV); Courtemanche v1/2=-9 S=5.8 (mV)
    double taud = dss * (1. - std::exp(-(v_new + 3*this->prm_[ISO] + 9.) / 6.)) / (0.035*(v_new + 3.*this->prm_[ISO] + 9.)); 
    this->var_[d_o] = ALGORITHM::RushLarsen(dss, this->var_[d_o], dt, taud);

    double fss = 1. / (1. + std::exp((v_new + 3.*this->prm_[ISO] + 30.) / 7.)) + 0.2 / (1. + std::exp((50.-v_new-3.*this->prm_[ISO]) / 20.)); 
    double tauf = 1. / (0.0197 * std::exp(-std::pow(0.0337*(v_new + 3.*this->prm_[ISO]+25.), 2.)) + 0.02);
    this->var_[fo] = ALGORITHM::RushLarsen(fss, this->var_[fo], dt, tauf);

    double dfcaBjo = 1.7*this->prm_[Cajo]*(1.0-this->var_[fcaBjo]) - 0.0119*this->var_[fcaBjo];
    this->var_[fcaBjo] = ALGORITHM::ForwardEuler(this->var_[fcaBjo], dt, dfcaBjo);

    double dfcaBslo = 1.7*this->prm_[Caslo]*(1.0-this->var_[fcaBslo]) - 0.0119*this->var_[fcaBslo];
    this->var_[fcaBslo] = ALGORITHM::ForwardEuler(this->var_[fcaBslo], dt, dfcaBslo);

    double fcaCaMSL = 0.;
    double fcaCaj = 0.;
    double ibarca_j  = 4.*pCa*(v_new*this->prm_[Fdy]*FoRT) * (0.341*this->prm_[Cajo]*std::exp(2.*v_new*FoRT) - 0.341*Cao) / (std::exp(2.*v_new*FoRT) - 1.);
    double ibarca_sl = 4.*pCa*(v_new*this->prm_[Fdy]*FoRT) * (0.341*this->prm_[Caslo]*std::exp(2.*v_new*FoRT) - 0.341*Cao) / (std::exp(2.*v_new*FoRT) - 1.);
    double ibark = pK*(v_new*this->prm_[Fdy]*FoRT)*(0.75*this->prm_[Kio]*std::exp(v_new*FoRT) - 0.75*Ko) / (std::exp(v_new*FoRT) - 1.);
    double ibarna_j = pNa*(v_new*this->prm_[Fdy]*FoRT) * (0.75*this->prm_[Najo]*std::exp(v_new*FoRT) - 0.75*Nao)  / (std::exp(v_new*FoRT) - 1.);
    double ibarna_sl = pNa*(v_new*this->prm_[Fdy]*FoRT) *(0.75*this->prm_[Naslo]*std::exp(v_new*FoRT) - 0.75*Nao)  / (std::exp(v_new*FoRT) - 1.);

    // Compute calcium and calcium-potassiun currents.
    double ICa_junc = 0.45*(Fjunc_CaL * ibarca_j * this->var_[d_o] * this->var_[fo] * ((1.-this->var_[fcaBjo])+fcaCaj) * std::pow(Q10CaL, Qpow));
    double ICa_sl = 0.45*(Fsl_CaL * ibarca_sl * this->var_[d_o]*this->var_[fo] * ((1.-this->var_[fcaBslo]) + fcaCaMSL) * std::pow(Q10CaL, Qpow));
    this->cur_[ICa] = (1.0-this->block_coeff_[ICa]) * (ICa_junc + ICa_sl);
    this->cur_[ICaK] = (1.0-this->block_coeff_[ICaK]) * (0.45*(ibark * this->var_[d_o] * this->var_[fo] * (Fjunc_CaL * (fcaCaj + (1.-this->var_[fcaBjo])) + Fsl_CaL*(fcaCaMSL + (1.-this->var_[fcaBslo]))) * std::pow(Q10CaL, Qpow)));
    
    // Compute calcium sodium current.
    double ICaNa_junc = 0.45*(Fjunc_CaL * ibarna_j * this->var_[d_o] * this->var_[fo] * ((1.-this->var_[fcaBjo]) + fcaCaj) * std::pow(Q10CaL, Qpow));
    double ICaNa_sl = 0.45*(Fsl_CaL * ibarna_sl * this->var_[d_o] * this->var_[fo] * ((1.-this->var_[fcaBslo]) + fcaCaMSL) * std::pow(Q10CaL, Qpow));
    this->cur_[ICaNa]= (1.0-this->block_coeff_[ICaNa]) * (ICaNa_junc + ICaNa_sl);
    this->cur_[ICatot] = (1.0-this->block_coeff_[ICatot]) * (this->cur_[ICa] + this->cur_[ICaK] + this->cur_[ICaNa]);

    // Compute sodium/calcium exchanger flux current.
    double Ka_junc = 1. / (1. + (Kdact/this->prm_[Cajo])*(Kdact/this->prm_[Cajo]));
    double Ka_sl = 1. / (1. + (Kdact/this->prm_[Caslo])*(Kdact/this->prm_[Caslo]));
    double s1_junc = std::exp(nu*v_new*FoRT) * this->prm_[Najo]*this->prm_[Najo]*this->prm_[Najo] * Cao;
    double s1_sl = std::exp(nu*v_new*FoRT) * this->prm_[Naslo]*this->prm_[Naslo]*this->prm_[Naslo] * Cao;
    double s2_junc = std::exp((nu-1.)*v_new*FoRT) * Nao*Nao*Nao * this->prm_[Cajo];
    double s3_junc = KmCai * Nao*Nao*Nao * (1. + std::pow(this->prm_[Najo]/KmNai, 3.)) + KmNao*KmNao*KmNao * this->prm_[Cajo] * 
                    (1. + this->prm_[Cajo]/KmCai) + KmCao*std::pow(this->prm_[Najo], 3.) + std::pow(this->prm_[Najo], 3.)*Cao + Nao*Nao*Nao*this->prm_[Cajo];
    double s2_sl = std::exp((nu-1)*v_new*FoRT) * Nao*Nao*Nao*this->prm_[Caslo];
    double s3_sl = KmCai*Nao*Nao*Nao*(1. + std::pow(this->prm_[Naslo]/KmNai, 3.)) + KmNao*KmNao*KmNao*this->prm_[Caslo]*(1. + this->prm_[Caslo]/KmCai) + 
                   KmCao*this->prm_[Naslo]*this->prm_[Naslo]*this->prm_[Naslo] + Cao*this->prm_[Naslo]*this->prm_[Naslo]*this->prm_[Naslo] + Nao*Nao*Nao*this->prm_[Caslo];

    double Incx_junc = Fjunc * IbarNCX * std::pow(Q10NCX, Qpow) * Ka_junc*(s1_junc-s2_junc) / s3_junc / (1.+ksat*std::exp((nu-1)*v_new*FoRT));
    double Incx_sl = Fsl * IbarNCX * std::pow(Q10NCX, Qpow) * Ka_sl * (s1_sl-s2_sl) / s3_sl / (1.+ksat*std::exp((nu-1)*v_new*FoRT));
    this->cur_[Incx] = (1.0-this->block_coeff_[Incx]) * (Incx_junc + Incx_sl);

    // Compute the sarcolemmal calcium pump current.
    double Ipca_junc = Fjunc * std::pow(Q10SLCaP, Qpow) * IbarSLCaP * std::pow(this->prm_[Cajo], 1.6) / (std::pow(KmPCa, 1.6) + std::pow(this->prm_[Cajo], 1.6));
    double Ipca_sl = Fsl * std::pow(Q10SLCaP, Qpow) * IbarSLCaP * std::pow(this->prm_[Caslo], 1.6) / (std::pow(KmPCa, 1.6) + std::pow(this->prm_[Caslo], 1.6));
    this->cur_[Ipca] = (1.0-this->block_coeff_[Ipca]) * (Ipca_junc + Ipca_sl);

    // Compute calcium background current.
    double ICab_junc = Fjunc * GCaB * (v_new-eca_junc);
    double ICab_sl = Fsl * GCaB * (v_new-eca_sl);
    this->cur_[ICab] = (1.0-this->block_coeff_[ICab]) * (ICab_junc + ICab_sl);

    // Compute sarcomere fluxes: Calcium Release, SR Ca pump, SR Ca leak.
    double MaxSR = 15; 
    double MinSR = 1;
    double kCaSR = MaxSR - (MaxSR-MinSR) / (1. + std::pow(ec50SR/this->prm_[Ca_sro], 2.5));
    double koSRCa = koCa / kCaSR;
    double kiSRCa = kiCa * kCaSR;
    double RI = 1. - this->var_[RyRro] - this->var_[RyRoo] - this->var_[RyRio];
    double dRyRro = (kim*RI - kiSRCa*this->prm_[Cajo]*this->var_[RyRro]) - (koSRCa*this->prm_[Cajo]*this->prm_[Cajo]*this->var_[RyRro] - kom*this->var_[RyRoo]);
    this->var_[RyRro] = ALGORITHM::ForwardEuler(this->var_[RyRro], dt, dRyRro);

    double dRyRoo = (koSRCa*this->prm_[Cajo]*this->prm_[Cajo]*this->var_[RyRro] - kom*this->var_[RyRoo]) - (kiSRCa*this->prm_[Cajo]*this->var_[RyRoo] - kim*this->var_[RyRio]);
    this->var_[RyRoo] = ALGORITHM::ForwardEuler(this->var_[RyRoo], dt, dRyRoo);

    double dRyRio = (kiSRCa*this->prm_[Cajo]*this->var_[RyRoo] - kim*this->var_[RyRio]) - (kom*this->var_[RyRio] - koSRCa*this->prm_[Cajo]*this->prm_[Cajo]*RI);
    this->var_[RyRio] = ALGORITHM::ForwardEuler(this->var_[RyRio], dt, dRyRio);

    double J_SRCarel = ks * this->var_[RyRoo] * (this->prm_[Ca_sro] - this->prm_[Cajo]);
    double J_serca = 1. * std::pow(Q10SRCaP, Qpow) * Vmax_SRCaP * (std::pow(this->prm_[Caio]/Kmf, hillSRCaP) - std::pow(this->prm_[Ca_sro]/Kmr, hillSRCaP)) /
                    (1. + std::pow(this->prm_[Caio]/Kmf, hillSRCaP) + std::pow(this->prm_[Ca_sro]/Kmr, hillSRCaP));
    double J_SRleak = (1. + 0.25*this->prm_[AF]) * 0.000005348 * (this->prm_[Ca_sro] - this->prm_[Cajo]);   // [mM/ms]

    // Sodium and Calcium Buffering
    double dNaBjo = kon_na * this->prm_[Najo] * (Bmax_Naj-this->var_[NaBjo]) - koff_na*this->var_[NaBjo];
    this->var_[NaBjo] = ALGORITHM::ForwardEuler(this->var_[NaBjo], dt, dNaBjo);

    double dNaBslo = kon_na * this->prm_[Naslo] * (Bmax_Nasl-this->var_[NaBslo]) - koff_na*this->var_[NaBslo];
    this->var_[NaBslo] = ALGORITHM::ForwardEuler(this->var_[NaBslo], dt, dNaBslo);

    // Cytosolic Ca Buffers
    double dTnCLo = kon_tncl * this->prm_[Caio] * (Bmax_TnClow - this->var_[TnCLo]) - koff_tncl*this->var_[TnCLo];
    this->var_[TnCLo] = ALGORITHM::ForwardEuler(this->var_[TnCLo], dt, dTnCLo);

    double dTnCHco = kon_tnchca * this->prm_[Caio] * (Bmax_TnChigh - this->var_[TnCHco] - this->var_[TnCHmo]) - koff_tnchca*this->var_[TnCHco];
    this->var_[TnCHco] = ALGORITHM::ForwardEuler(this->var_[TnCHco], dt, dTnCHco);
    
    double dTnCHmo = kon_tnchmg * Mgi * (Bmax_TnChigh - this->var_[TnCHco] - this->var_[TnCHmo]) - koff_tnchmg*this->var_[TnCHmo];
    this->var_[TnCHmo] = ALGORITHM::ForwardEuler(this->var_[TnCHmo], dt, dTnCHmo);

    double dCaMo = kon_cam * this->prm_[Caio] * (Bmax_CaM - this->var_[CaMo]) - koff_cam*this->var_[CaMo];
    this->var_[CaMo] = ALGORITHM::ForwardEuler(this->var_[CaMo], dt, dCaMo);

    double dMyoco = kon_myoca * this->prm_[Caio] * (Bmax_myosin - this->var_[Myoco] - this->var_[Myomo]) - koff_myoca*this->var_[Myoco];
    this->var_[Myoco] = ALGORITHM::ForwardEuler(this->var_[Myoco], dt, dMyoco);

    double dMyomo = kon_myomg * Mgi * (Bmax_myosin - this->var_[Myoco] - this->var_[Myomo]) - koff_myomg*this->var_[Myomo];
    this->var_[Myomo] = ALGORITHM::ForwardEuler(this->var_[Myomo], dt, dMyomo);

    double dSRBo = kon_sr * this->prm_[Caio] * (Bmax_SR - this->var_[SRBo]) - koff_sr*this->var_[SRBo];
    this->var_[SRBo] = ALGORITHM::ForwardEuler(this->var_[SRBo], dt, dSRBo);
    double J_CaB_cytosol = dTnCLo + dTnCHco + dTnCHmo + dCaMo + dMyoco + dMyomo + dSRBo; 

    // Junctional and SL Ca Buffers.
    double dSLLjo  = kon_sll * this->prm_[Cajo] * (Bmax_SLlowj - this->var_[SLLjo]) - koff_sll*this->var_[SLLjo];
    this->var_[SLLjo] = ALGORITHM::ForwardEuler(this->var_[SLLjo], dt, dSLLjo);

    double dSLLslo = kon_sll * this->prm_[Caslo] * (Bmax_SLlowsl - this->var_[SLLslo]) - koff_sll*this->var_[SLLslo];
    this->var_[SLLslo] = ALGORITHM::ForwardEuler(this->var_[SLLslo], dt, dSLLslo);

    double dSLHjo  = kon_slh * this->prm_[Cajo] * (Bmax_SLhighj - this->var_[SLHjo]) - koff_slh*this->var_[SLHjo];
    this->var_[SLHjo] = ALGORITHM::ForwardEuler(this->var_[SLHjo], dt, dSLHjo);
    
    double dSLHslo = kon_slh * this->prm_[Caslo] * (Bmax_SLhighsl - this->var_[SLHslo]) - koff_slh*this->var_[SLHslo];
    this->var_[SLHslo] = ALGORITHM::ForwardEuler(this->var_[SLHslo], dt, dSLHslo);

    double J_CaB_junction = dSLLjo + dSLHjo;
    double J_CaB_sl = dSLLslo + dSLHslo;

    // Sarcomere calcium concentrations.
    double dCsqnbo = kon_csqn * this->prm_[Ca_sro] * (Bmax_Csqn - this->var_[Csqnbo]) - koff_csqn*this->var_[Csqnbo];
    this->var_[Csqnbo] = ALGORITHM::ForwardEuler(this->var_[Csqnbo], dt, dCsqnbo);

    double dCa_sro = J_serca - (J_SRleak*Vmyo/Vsr + J_SRCarel) - dCsqnbo;
    this->prm_[Ca_sro] = ALGORITHM::ForwardEuler(this->prm_[Ca_sro], dt, dCa_sro);
    
    // Compute sodium concentrations.
    double INa_tot_junc = INa_junc + INab_junc + 3.*Incx_junc + 3.*INaK_junc + ICaNa_junc + INaL_junc;
    double INa_tot_sl = INa_sl + INab_sl + 3.*Incx_sl + 3.*INaK_sl + ICaNa_sl + INaL_sl;

    double dNajo = -INa_tot_junc * this->prm_[Cmem] / (Vjunc*this->prm_[Fdy]) + J_na_juncsl/Vjunc*(this->prm_[Naslo] - this->prm_[Najo]) - dNaBjo;
    this->prm_[Najo] = ALGORITHM::ForwardEuler(this->prm_[Najo], dt, dNajo);

    double dNaslo = -INa_tot_sl * this->prm_[Cmem] / (Vsl*this->prm_[Fdy]) + J_na_juncsl/Vsl*(this->prm_[Najo] - this->prm_[Naslo]) + J_na_slmyo/Vsl*(this->prm_[Naio]-this->prm_[Naslo])-dNaBslo;
    this->prm_[Naslo] = ALGORITHM::ForwardEuler(this->prm_[Naslo], dt, dNaslo);

    double dNaio = J_na_slmyo / Vmyo * (this->prm_[Naslo] - this->prm_[Naio]);
    this->prm_[Naio] = ALGORITHM::ForwardEuler(this->prm_[Naio], dt, dNaio);

    // Compute potassium concentration.
    double IK_tot = this->cur_[Ito] + this->cur_[IKr] + this->cur_[IKs] + this->cur_[IKi] - 2.*this->cur_[INaK] + 
                     this->cur_[ICaK] + this->cur_[IKp] + this->cur_[IKur] + this->cur_[IKach] + this->cur_[IKCa];
    // double dKio = 0.; 
    // this->prm_[Kio] = ALGORITHM::ForwardEuler(this->prm_[Kio], dt, dKio);

    // Compute calcium concentrations.
    double ICa_tot_junc = ICa_junc + ICab_junc + Ipca_junc - 2.*Incx_junc;
    double ICa_tot_sl = ICa_sl + ICab_sl + Ipca_sl - 2.*Incx_sl;
    double dCajo = -ICa_tot_junc*this->prm_[Cmem]/(Vjunc*2.*this->prm_[Fdy]) + J_ca_juncsl/Vjunc*(this->prm_[Caslo]-this->prm_[Cajo]) - 
                     J_CaB_junction + J_SRCarel*Vsr/Vjunc + J_SRleak*Vmyo/Vjunc;
    this->prm_[Cajo] = ALGORITHM::ForwardEuler(this->prm_[Cajo], dt, dCajo);

    double dCaslo = -ICa_tot_sl*this->prm_[Cmem]/(Vsl*2*this->prm_[Fdy]) + J_ca_juncsl/Vsl*(this->prm_[Cajo] - this->prm_[Caslo]) + 
                      J_ca_slmyo/Vsl*(this->prm_[Caio] - this->prm_[Caslo]) - J_CaB_sl;
    this->prm_[Caslo] = ALGORITHM::ForwardEuler(this->prm_[Caslo], dt, dCaslo);

    double dCaio = -J_serca*Vsr/Vmyo - J_CaB_cytosol +J_ca_slmyo/Vmyo*(this->prm_[Caslo] - this->prm_[Caio]);
    this->prm_[Caio] = ALGORITHM::ForwardEuler(this->prm_[Caio], dt, dCaio);

    // Membrane Potential
    double INa_tot = INa_tot_junc + INa_tot_sl;
    double ICl_tot = IClCa + this->cur_[IClb];
    double ICa_tot = ICa_tot_junc + ICa_tot_sl;

    // Add zero stimulus current to avoid compilation warning.
    this->cur_[GrdAtrCur::Iion] = INa_tot + ICl_tot + ICa_tot + IK_tot;

    // Set the new value to the membrante potential time derivative.
    this->var_[dvdt] = - (this->cur_[GrdAtrCur::Iion] - stim_current);

}


std::string GrandiAtri::PrintVariables() const 
{
    using namespace GrdAtrVar;

    // Create output string stream to pass the variables and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "v = " << this->var_[v] << "\n";
    oss << "dvdt = " << this->var_[dvdt] << "\n";
    oss << "mo = " << this->var_[mo] << "\n";
    oss << "ho = " << this->var_[ho] << "\n";
    oss << "jo = " << this->var_[jo] << "\n";
    oss << "d_o = " << this->var_[d_o] << "\n";
    oss << "fo = " << this->var_[fo] << "\n";
    oss << "fcaBjo = " << this->var_[fcaBjo] << "\n";
    oss << "fcaBslo = " << this->var_[fcaBslo] << "\n";
    oss << "xtoso = " << this->var_[xtoso] << "\n";
    oss << "ytoso = " << this->var_[ytoso] << "\n";
    oss << "xtofo = " << this->var_[xtofo] << "\n";
    oss << "ytofo = " << this->var_[ytofo] << "\n";
    oss << "xkro = " << this->var_[xkro] << "\n";
    oss << "xkso = " << this->var_[xkso] << "\n";
    oss << "RyRro = " << this->var_[RyRro] << "\n";
    oss << "RyRoo = " << this->var_[RyRoo] << "\n";
    oss << "RyRio = " << this->var_[RyRio] << "\n";
    oss << "NaBjo = " << this->var_[NaBjo] << "\n";
    oss << "NaBslo = " << this->var_[NaBslo] << "\n";
    oss << "TnCLo = " << this->var_[TnCLo] << "\n";
    oss << "TnCHco = " << this->var_[TnCHco] << "\n";
    oss << "TnCHmo = " << this->var_[TnCHmo] << "\n";
    oss << "CaMo = " << this->var_[CaMo] << "\n";
    oss << "Myoco = " << this->var_[Myoco] << "\n";
    oss << "Myomo = " << this->var_[Myomo] << "\n";
    oss << "SRBo = " << this->var_[SRBo] << "\n";
    oss << "SLLjo = " << this->var_[SLLjo] << "\n";
    oss << "SLLslo = " << this->var_[SLLslo] << "\n";
    oss << "SLHjo = " << this->var_[SLHjo] << "\n";
    oss << "SLHslo = " << this->var_[SLHslo] << "\n";
    oss << "Csqnbo = " << this->var_[Csqnbo];
    return oss.str();

}


std::string GrandiAtri::PrintParameters() const 
{
    using namespace GrdAtrPrm;    

    // Create output string stream to pass the parameters and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "Ca_sro = " << this->prm_[Ca_sro] << "\n";
    oss << "Najo = " << this->prm_[Najo] << "\n";
    oss << "Naslo = " << this->prm_[Naslo] << "\n";
    oss << "Naio = " << this->prm_[Naio] << "\n";
    oss << "Kio = " << this->prm_[Kio] << "\n";
    oss << "Cajo = " << this->prm_[Cajo] << "\n";
    oss << "Caslo = " << this->prm_[Caslo] << "\n";
    oss << "Caio = " << this->prm_[Caio] << "\n";
    oss << "rtoso = " << this->prm_[rtoso] << "\n";
    oss << "n41 = " << this->prm_[n41] << "\n";
    oss << "C1o = " << this->prm_[C1o] << "\n";
    oss << "C2o = " << this->prm_[C2o] << "\n";
    oss << "C3o = " << this->prm_[C3o] << "\n";
    oss << "C4o = " << this->prm_[C4o] << "\n";
    oss << "C5o = " << this->prm_[C5o] << "\n";
    oss << "C6o = " << this->prm_[C6o] << "\n";
    oss << "C7o = " << this->prm_[C7o] << "\n";
    oss << "C8o = " << this->prm_[C8o] << "\n";
    oss << "C9o = " << this->prm_[C9o] << "\n";
    oss << "C10o = " << this->prm_[C10o] << "\n";
    oss << "C11o = " << this->prm_[C11o] << "\n";
    oss << "C12o = " << this->prm_[C12o] << "\n";
    oss << "C13o = " << this->prm_[C13o] << "\n";
    oss << "C14o = " << this->prm_[C14o] << "\n";
    oss << "C15o = " << this->prm_[C15o] << "\n";
    oss << "O1o = " << this->prm_[O1o] << "\n";
    oss << "rkuro = " << this->prm_[rkuro] << "\n";
    oss << "skuro = " << this->prm_[skuro] << "\n";
    oss << "n60 = " << this->prm_[n60] << "\n";
    oss << "n61 = " << this->prm_[n61] << "\n";
    oss << "n62 = " << this->prm_[n62] << "\n";
    oss << "ikcaoo = " << this->prm_[ikcaoo] << "\n";
    oss << "Ach = " << this->prm_[Ach] << "\n";
    oss << "ISO = " << this->prm_[ISO] << "\n";
    oss << "ISOdlf1muM = " << this->prm_[ISOdlf1muM] << "\n";
    oss << "ISOdlf1nM = " << this->prm_[ISOdlf1nM] << "\n";
    oss << "drugvector = " << this->prm_[drugvector] << "\n";
    oss << "AF = " << this->prm_[AF] << "\n";
    oss << "RA = " << this->prm_[RA] << "\n";
    oss << "R = " << this->prm_[R] << "\n";
    oss << "Fdy = " << this->prm_[Fdy] << "\n";
    oss << "T = " << this->prm_[T] << "\n";
    oss << "l = " << this->prm_[l] << "\n";
    oss << "rad = " << this->prm_[rad] << "\n";
    oss << "Cmem = " << this->prm_[Cmem];
    return oss.str();

}


std::string GrandiAtri::PrintCurrents() const
{
    using namespace GrdAtrCur;

    // Create output string stream to pass the currents and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "INa = " << this->cur_[INa] << "\n";
    oss << "INaL = " << this->cur_[INaL] << "\n";
    oss << "INab = " << this->cur_[INab] << "\n";
    oss << "INaK = " << this->cur_[INaK] << "\n";
    oss << "IKr = " << this->cur_[IKr] << "\n";
    oss << "IKs = " << this->cur_[IKs] << "\n";
    oss << "IKp = " << this->cur_[IKp] << "\n";
    oss << "IKach = " << this->cur_[IKach] << "\n";
    oss << "IKCa = " << this->cur_[IKCa] << "\n";
    oss << "Ito = " << this->cur_[Ito] << "\n";
    oss << "IKur = " << this->cur_[IKur] << "\n";
    oss << "IKi = " << this->cur_[IKi] << "\n";
    oss << "IClb = " << this->cur_[IClb] << "\n";
    oss << "ICa = " << this->cur_[ICa] << "\n";
    oss << "ICaK = " << this->cur_[ICaK] << "\n";
    oss << "ICaNa = " << this->cur_[ICaNa] << "\n";
    oss << "ICatot = " << this->cur_[ICatot] << "\n";
    oss << "Incx = " << this->cur_[Incx] << "\n";
    oss << "Ipca = " << this->cur_[Ipca] << "\n";
    oss << "ICab = " << this->cur_[ICab] << "\n";
    oss << "Iion = " << this->cur_[GrdAtrCur::Iion]; 
    return oss.str();

}


std::string GrandiAtri::PrintBlockCoeffs() const
{
    using namespace GrdAtrCur;

    // Create output string stream to pass the currents and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "INa = " << this->block_coeff_[INa] << "\n";
    oss << "INaL = " << this->block_coeff_[INaL] << "\n";
    oss << "INab = " << this->block_coeff_[INab] << "\n";
    oss << "INaK = " << this->block_coeff_[INaK] << "\n";
    oss << "IKr = " << this->block_coeff_[IKr] << "\n";
    oss << "IKs = " << this->block_coeff_[IKs] << "\n";
    oss << "IKp = " << this->block_coeff_[IKp] << "\n";
    oss << "IKach = " << this->block_coeff_[IKach] << "\n";
    oss << "IKCa = " << this->block_coeff_[IKCa] << "\n";
    oss << "Ito = " << this->block_coeff_[Ito] << "\n";
    oss << "IKur = " << this->block_coeff_[IKur] << "\n";
    oss << "IKi = " << this->block_coeff_[IKi] << "\n";
    oss << "IClb = " << this->block_coeff_[IClb] << "\n";
    oss << "ICa = " << this->block_coeff_[ICa] << "\n";
    oss << "ICaK = " << this->block_coeff_[ICaK] << "\n";
    oss << "ICaNa = " << this->block_coeff_[ICaNa] << "\n";
    oss << "ICatot = " << this->block_coeff_[ICatot] << "\n";
    oss << "Incx = " << this->block_coeff_[Incx] << "\n";
    oss << "Ipca = " << this->block_coeff_[Ipca] << "\n";
    oss << "ICab = " << this->block_coeff_[ICab];
    return oss.str();

}



} // End of namespace ELECTRA