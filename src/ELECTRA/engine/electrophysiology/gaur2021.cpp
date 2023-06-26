/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

#include "ELECTRA/engine/electrophysiology/gaur2021.hpp"


namespace ELECTRA {

void Gaur2021::SetDataMapping()
{
    using namespace Gaur21Var;
    using namespace Gaur21Prm;
    using namespace Gaur21Cur;

    // Set variables mapping.
    this->mapped_data_["v"] = static_cast<std::size_t>(v);
    this->mapped_data_["CICR__A"] = static_cast<std::size_t>(CICR__A);
    this->mapped_data_["CICR__Jrel1"] = static_cast<std::size_t>(CICR__Jrel1);
    this->mapped_data_["CICR__Jrel2"] = static_cast<std::size_t>(CICR__Jrel2);
    this->mapped_data_["ionic_concentrations__cacsr"] = static_cast<std::size_t>(ionic_concentrations__cacsr);
    this->mapped_data_["ionic_concentrations__cai"] = static_cast<std::size_t>(ionic_concentrations__cai);
    this->mapped_data_["ionic_concentrations__cai2"] = static_cast<std::size_t>(ionic_concentrations__cai2);
    this->mapped_data_["ionic_concentrations__cajsr"] = static_cast<std::size_t>(ionic_concentrations__cajsr);
    this->mapped_data_["ionic_concentrations__cass"] = static_cast<std::size_t>(ionic_concentrations__cass);
    this->mapped_data_["CICR__tjsrol"] = static_cast<std::size_t>(CICR__tjsrol);
    this->mapped_data_["CaMK__CaMKt"] = static_cast<std::size_t>(CaMK__CaMKt);
    this->mapped_data_["ionic_concentrations__ki"] = static_cast<std::size_t>(ionic_concentrations__ki);
    this->mapped_data_["ionic_concentrations__nai"] = static_cast<std::size_t>(ionic_concentrations__nai);
    this->mapped_data_["ICaL__d"] = static_cast<std::size_t>(ICaL__d);
    this->mapped_data_["ICaL__fca"] = static_cast<std::size_t>(ICaL__fca);
    this->mapped_data_["ICaL__ff"] = static_cast<std::size_t>(ICaL__ff);
    this->mapped_data_["ICaL__fs"] = static_cast<std::size_t>(ICaL__fs);
    this->mapped_data_["ionic_concentrations__kss"] = static_cast<std::size_t>(ionic_concentrations__kss);
    this->mapped_data_["ionic_concentrations__nass"] = static_cast<std::size_t>(ionic_concentrations__nass);
    this->mapped_data_["IKr__xr"] = static_cast<std::size_t>(IKr__xr);
    this->mapped_data_["IKs__xs1"] = static_cast<std::size_t>(IKs__xs1);
    this->mapped_data_["IKs__xs2"] = static_cast<std::size_t>(IKs__xs2);
    this->mapped_data_["INaL__hl"] = static_cast<std::size_t>(INaL__hl);
    this->mapped_data_["INaL__ml"] = static_cast<std::size_t>(INaL__ml);
    this->mapped_data_["ITo__aa"] = static_cast<std::size_t>(ITo__aa);
    this->mapped_data_["I_Na__h"] = static_cast<std::size_t>(I_Na__h);
    this->mapped_data_["I_Na__j"] = static_cast<std::size_t>(I_Na__j);
    this->mapped_data_["I_Na__m"] = static_cast<std::size_t>(I_Na__m);
    this->mapped_data_["ionic_concentrations__cansr"] = static_cast<std::size_t>(ionic_concentrations__cansr);


    // Set parameters mapping.
    this->mapped_data_["cell__F"] = static_cast<std::size_t>(cell__F);
    this->mapped_data_["CaMK__KmCaMK"] = static_cast<std::size_t>(CaMK__KmCaMK);
    this->mapped_data_["CICR__SOICR"] = static_cast<std::size_t>(CICR__SOICR);
    this->mapped_data_["CICR__grelbarjsrol"] = static_cast<std::size_t>(CICR__grelbarjsrol);
    this->mapped_data_["CICR__tau_gap"] = static_cast<std::size_t>(CICR__tau_gap);
    this->mapped_data_["CICR__tauoff"] = static_cast<std::size_t>(CICR__tauoff);
    this->mapped_data_["CICR__tauon"] = static_cast<std::size_t>(CICR__tauon);
    this->mapped_data_["CaMK__CaMKo"] = static_cast<std::size_t>(CaMK__CaMKo);
    this->mapped_data_["CaMK__KmCaM"] = static_cast<std::size_t>(CaMK__KmCaM);
    this->mapped_data_["CaMK__PKNa"] = static_cast<std::size_t>(CaMK__PKNa);
    this->mapped_data_["cell__R"] = static_cast<std::size_t>(cell__R);
    this->mapped_data_["cell__T"] = static_cast<std::size_t>(cell__T);
    this->mapped_data_["CaMK__aCaMK"] = static_cast<std::size_t>(CaMK__aCaMK);
    this->mapped_data_["CaMK__bCaMK"] = static_cast<std::size_t>(CaMK__bCaMK);
    this->mapped_data_["cell__cli"] = static_cast<std::size_t>(cell__cli);
    this->mapped_data_["cell__clo"] = static_cast<std::size_t>(cell__clo);
    this->mapped_data_["cell__ko"] = static_cast<std::size_t>(cell__ko);
    this->mapped_data_["cell__nao"] = static_cast<std::size_t>(cell__nao);
    this->mapped_data_["ICaL__PCa"] = static_cast<std::size_t>(ICaL__PCa);
    this->mapped_data_["ICaL__vhalf_d"] = static_cast<std::size_t>(ICaL__vhalf_d);
    this->mapped_data_["ICaL__zca"] = static_cast<std::size_t>(ICaL__zca);
    this->mapped_data_["IKb__GKb"] = static_cast<std::size_t>(IKb__GKb);
    this->mapped_data_["IKs__GKs"] = static_cast<std::size_t>(IKs__GKs);
    this->mapped_data_["INaCa_i__KmCaAct"] = static_cast<std::size_t>(INaCa_i__KmCaAct);
    this->mapped_data_["INaCa_i__kasymm"] = static_cast<std::size_t>(INaCa_i__kasymm);
    this->mapped_data_["INaCa_i__kcaoff"] = static_cast<std::size_t>(INaCa_i__kcaoff);
    this->mapped_data_["INaCa_i__kcaon"] = static_cast<std::size_t>(INaCa_i__kcaon);
    this->mapped_data_["INaCa_i__kna1"] = static_cast<std::size_t>(INaCa_i__kna1);
    this->mapped_data_["INaCa_i__kna2"] = static_cast<std::size_t>(INaCa_i__kna2);
    this->mapped_data_["INaCa_i__kna3"] = static_cast<std::size_t>(INaCa_i__kna3);
    this->mapped_data_["INaCa_i__qca"] = static_cast<std::size_t>(INaCa_i__qca);
    this->mapped_data_["INaCa_i__qna"] = static_cast<std::size_t>(INaCa_i__qna);
    this->mapped_data_["INaCa_i__wca"] = static_cast<std::size_t>(INaCa_i__wca);
    this->mapped_data_["INaCa_i__wna"] = static_cast<std::size_t>(INaCa_i__wna);
    this->mapped_data_["INaCa_i__wnaca"] = static_cast<std::size_t>(INaCa_i__wnaca);
    this->mapped_data_["INaCa_i__zna"] = static_cast<std::size_t>(INaCa_i__zna);
    this->mapped_data_["INaK__H"] = static_cast<std::size_t>(INaK__H);
    this->mapped_data_["INaK__Khp"] = static_cast<std::size_t>(INaK__Khp);
    this->mapped_data_["INaK__Kki"] = static_cast<std::size_t>(INaK__Kki);
    this->mapped_data_["INaK__Kko"] = static_cast<std::size_t>(INaK__Kko);
    this->mapped_data_["INaK__Kmgatp"] = static_cast<std::size_t>(INaK__Kmgatp);
    this->mapped_data_["INaK__Knai0"] = static_cast<std::size_t>(INaK__Knai0);
    this->mapped_data_["INaK__Knao0"] = static_cast<std::size_t>(INaK__Knao0);
    this->mapped_data_["INaK__Knap"] = static_cast<std::size_t>(INaK__Knap);
    this->mapped_data_["INaK__Kxkur"] = static_cast<std::size_t>(INaK__Kxkur);
    this->mapped_data_["INaK__MgADP"] = static_cast<std::size_t>(INaK__MgADP);
    this->mapped_data_["INaK__MgATP"] = static_cast<std::size_t>(INaK__MgATP);
    this->mapped_data_["INaK__eP"] = static_cast<std::size_t>(INaK__eP);
    this->mapped_data_["INaK__k1m"] = static_cast<std::size_t>(INaK__k1m);
    this->mapped_data_["INaK__k1p"] = static_cast<std::size_t>(INaK__k1p);
    this->mapped_data_["INaK__k2m"] = static_cast<std::size_t>(INaK__k2m);
    this->mapped_data_["INaK__k2p"] = static_cast<std::size_t>(INaK__k2p);
    this->mapped_data_["INaK__k3m"] = static_cast<std::size_t>(INaK__k3m);
    this->mapped_data_["INaK__k3p2"] = static_cast<std::size_t>(INaK__k3p2);
    this->mapped_data_["INaK__k4m"] = static_cast<std::size_t>(INaK__k4m);
    this->mapped_data_["INaK__k4p2"] = static_cast<std::size_t>(INaK__k4p2);
    this->mapped_data_["INaK__zk"] = static_cast<std::size_t>(INaK__zk);
    this->mapped_data_["INaL__GNaL"] = static_cast<std::size_t>(INaL__GNaL);
    this->mapped_data_["INaL__tau_hl"] = static_cast<std::size_t>(INaL__tau_hl);
    this->mapped_data_["INab__PNab"] = static_cast<std::size_t>(INab__PNab);
    this->mapped_data_["I_Na__GNa"] = static_cast<std::size_t>(I_Na__GNa);
    this->mapped_data_["IpCa__GpCa"] = static_cast<std::size_t>(IpCa__GpCa);
    this->mapped_data_["SR_uptake__BSLmax"] = static_cast<std::size_t>(SR_uptake__BSLmax);
    this->mapped_data_["SR_uptake__BSRmax"] = static_cast<std::size_t>(SR_uptake__BSRmax);
    this->mapped_data_["SR_uptake__KmBSL"] = static_cast<std::size_t>(SR_uptake__KmBSL);
    this->mapped_data_["SR_uptake__KmBSR"] = static_cast<std::size_t>(SR_uptake__KmBSR);
    this->mapped_data_["SR_uptake__cmdnmax"] = static_cast<std::size_t>(SR_uptake__cmdnmax);
    this->mapped_data_["SR_uptake__csqnmax"] = static_cast<std::size_t>(SR_uptake__csqnmax);
    this->mapped_data_["SR_uptake__kmcmdn"] = static_cast<std::size_t>(SR_uptake__kmcmdn);
    this->mapped_data_["SR_uptake__kmcsqn"] = static_cast<std::size_t>(SR_uptake__kmcsqn);
    this->mapped_data_["SR_uptake__kmtrpn"] = static_cast<std::size_t>(SR_uptake__kmtrpn);
    this->mapped_data_["SR_uptake__trpnmax"] = static_cast<std::size_t>(SR_uptake__trpnmax);
    this->mapped_data_["cell__L"] = static_cast<std::size_t>(cell__L);
    this->mapped_data_["cell__pi"] = static_cast<std::size_t>(cell__pi);
    this->mapped_data_["cell__rad"] = static_cast<std::size_t>(cell__rad);
    this->mapped_data_["cell__vmyo1frac"] = static_cast<std::size_t>(cell__vmyo1frac);
    this->mapped_data_["stimulus__duration"] = static_cast<std::size_t>(stimulus__duration);
    this->mapped_data_["stimulus__offset"] = static_cast<std::size_t>(stimulus__offset);
    this->mapped_data_["stimulus__period"] = static_cast<std::size_t>(stimulus__period);
    this->mapped_data_["CaMK__ECl"] = static_cast<std::size_t>(CaMK__ECl);
    this->mapped_data_["ICaL__PCaK"] = static_cast<std::size_t>(ICaL__PCaK);
    this->mapped_data_["ICaL__PCaNa"] = static_cast<std::size_t>(ICaL__PCaNa);
    this->mapped_data_["ICaL__vhalff"] = static_cast<std::size_t>(ICaL__vhalff);
    this->mapped_data_["ICab__PCab"] = static_cast<std::size_t>(ICab__PCab);
    this->mapped_data_["IK1__GK1"] = static_cast<std::size_t>(IK1__GK1);
    this->mapped_data_["IKr__GKr"] = static_cast<std::size_t>(IKr__GKr);
    this->mapped_data_["INaK__delta"] = static_cast<std::size_t>(INaK__delta);
    this->mapped_data_["INaCa_i__Gncx"] = static_cast<std::size_t>(INaCa_i__Gncx);
    this->mapped_data_["ITo__Gto"] = static_cast<std::size_t>(ITo__Gto);
    this->mapped_data_["stimulus__stimulus_amplitude"] = static_cast<std::size_t>(stimulus__stimulus_amplitude);
    this->mapped_data_["INaCa_i__h10"] = static_cast<std::size_t>(INaCa_i__h10);
    this->mapped_data_["INaCa_i__h11"] = static_cast<std::size_t>(INaCa_i__h11);
    this->mapped_data_["INaCa_i__h12"] = static_cast<std::size_t>(INaCa_i__h12);
    this->mapped_data_["INaCa_i__k2"] = static_cast<std::size_t>(INaCa_i__k2);
    this->mapped_data_["INaCa_i__k5"] = static_cast<std::size_t>(INaCa_i__k5);
    this->mapped_data_["INaCa_ss__h101"] = static_cast<std::size_t>(INaCa_ss__h101);
    this->mapped_data_["INaCa_ss__h1111"] = static_cast<std::size_t>(INaCa_ss__h1111);
    this->mapped_data_["INaCa_ss__h121"] = static_cast<std::size_t>(INaCa_ss__h121);
    this->mapped_data_["INaCa_ss__k21"] = static_cast<std::size_t>(INaCa_ss__k21);
    this->mapped_data_["INaCa_ss__k51"] = static_cast<std::size_t>(INaCa_ss__k51);
    this->mapped_data_["INaK__Pnak"] = static_cast<std::size_t>(INaK__Pnak);
    this->mapped_data_["INaK__a2"] = static_cast<std::size_t>(INaK__a2);
    this->mapped_data_["INaK__a4"] = static_cast<std::size_t>(INaK__a4);
    this->mapped_data_["INaK__b1"] = static_cast<std::size_t>(INaK__b1);
    this->mapped_data_["cell__Ageo"] = static_cast<std::size_t>(cell__Ageo);
    this->mapped_data_["cell__Acap"] = static_cast<std::size_t>(cell__Acap);
    this->mapped_data_["cell__cao"] = static_cast<std::size_t>(cell__cao);
    this->mapped_data_["INaCa_i__k1"] = static_cast<std::size_t>(INaCa_i__k1);
    this->mapped_data_["INaCa_ss__k11"] = static_cast<std::size_t>(INaCa_ss__k11);
    this->mapped_data_["cell__vcell"] = static_cast<std::size_t>(cell__vcell);
    this->mapped_data_["cell__vjsr"] = static_cast<std::size_t>(cell__vjsr);
    this->mapped_data_["cell__vcsr"] = static_cast<std::size_t>(cell__vcsr);
    this->mapped_data_["cell__vmyo"] = static_cast<std::size_t>(cell__vmyo);
    this->mapped_data_["cell__vmyo1"] = static_cast<std::size_t>(cell__vmyo1);
    this->mapped_data_["cell__vmyo2"] = static_cast<std::size_t>(cell__vmyo2);
    this->mapped_data_["cell__vnsr"] = static_cast<std::size_t>(cell__vnsr);
    this->mapped_data_["cell__vnsr1"] = static_cast<std::size_t>(cell__vnsr1);
    this->mapped_data_["cell__vnsr2"] = static_cast<std::size_t>(cell__vnsr2);
    this->mapped_data_["cell__vss"] = static_cast<std::size_t>(cell__vss);


    // Set currents mapping.
    this->mapped_data_["INa"] = static_cast<std::size_t>(INa);
    this->mapped_data_["INaL"] = static_cast<std::size_t>(INaL);
    this->mapped_data_["ICaL"] = static_cast<std::size_t>(ICaL);
    this->mapped_data_["ICaNa"] = static_cast<std::size_t>(ICaNa);
    this->mapped_data_["ICaK"] = static_cast<std::size_t>(ICaK);
    this->mapped_data_["IKr"] = static_cast<std::size_t>(IKr);
    this->mapped_data_["IKs"] = static_cast<std::size_t>(IKs);
    this->mapped_data_["IK1"] = static_cast<std::size_t>(IK1);
    this->mapped_data_["ITo"] = static_cast<std::size_t>(ITo);
    this->mapped_data_["INaCa_i"] = static_cast<std::size_t>(INaCa_i);
    this->mapped_data_["INaCa_ss"] = static_cast<std::size_t>(INaCa_ss);
    this->mapped_data_["INaK"] = static_cast<std::size_t>(INaK);
    this->mapped_data_["INab"] = static_cast<std::size_t>(INab);
    this->mapped_data_["IKb"] = static_cast<std::size_t>(IKb);
    this->mapped_data_["IpCa"] = static_cast<std::size_t>(IpCa);
    this->mapped_data_["ICab"] = static_cast<std::size_t>(ICab);

}


Gaur2021::Gaur2021()
{
    //Initialize the model's data.
    this->model_type_ = EpModelType::Gaur2021;
    this->dt_stable_ = 0.01;
    this->var_.resize(30, 0.);
    this->prm_.resize(119, 0.);
    this->cur_.resize(17, 0.);
    this->block_coeff_.resize(16, 0.);

    // Set mapped data.
    this->SetDataMapping();
}


Gaur2021::~Gaur2021()
{}


void Gaur2021::Initialize(CellType cell_type)
{
    using namespace Gaur21Var;
    using namespace Gaur21Prm;

    //Initialize the model's data.
    this->var_.clear();           this->var_.resize(30, 0.);
    this->prm_.clear();           this->prm_.resize(119, 0.);
    this->cur_.clear();           this->cur_.resize(17, 0.);
    this->block_coeff_.clear();   this->block_coeff_.resize(16, 0.);

    // Set data according to the cell type.
    switch (cell_type)
    {
        case CellType::ventricular :
            // Set the model variables for ventricular cell type.
            this->var_[v] = -87.3;
            this->var_[dvdt]   = 0.;
            this->var_[CICR__A] = 100.0;
            this->var_[CICR__Jrel1] = 9.71e-22;
            this->var_[CICR__Jrel2] = 8.59e-10;
            this->var_[ionic_concentrations__cacsr] = 1.3;
            this->var_[ionic_concentrations__cai] = 6.75e-05;
            this->var_[ionic_concentrations__cai2] = 6.81e-05;
            this->var_[ionic_concentrations__cajsr] = 1.31;
            this->var_[ionic_concentrations__cass] = 6.75e-05;
            this->var_[CICR__tjsrol] = 100.0;
            this->var_[CaMK__CaMKt] = 8.23910999999999914e-03;
            this->var_[ionic_concentrations__ki] = 140.76;
            this->var_[ionic_concentrations__nai] = 6.43;
            this->var_[ICaL__d] = 4.42e-07;
            this->var_[ICaL__fca] = 0.988;
            this->var_[ICaL__ff] = 0.0;
            this->var_[ICaL__fs] = 0.0;
            this->var_[ionic_concentrations__kss] = 140.76;
            this->var_[ionic_concentrations__nass] = 6.43;
            this->var_[IKr__xr] = 0.153;
            this->var_[IKs__xs1] = 0.054;
            this->var_[IKs__xs2] = 0.046;
            this->var_[INaL__hl] = 0.301;
            this->var_[INaL__ml] = 0.001;
            this->var_[ITo__aa] = 0.995;
            this->var_[I_Na__h] = 0.631;
            this->var_[I_Na__j] = 0.631;
            this->var_[I_Na__m] = 0.0022;
            this->var_[ionic_concentrations__cansr] = 1.37;   

            // Set the model parameters for endocardium cell type.
            this->prm_[cell__F] = 96485.0;    
            this->prm_[CaMK__KmCaMK] = 0.065;
            this->prm_[CICR__SOICR] = 10.0;   
            this->prm_[CICR__grelbarjsrol] = 2.0;
            this->prm_[CICR__tau_gap] = 3.0;
            this->prm_[CICR__tauoff] = 5.0;
            this->prm_[CICR__tauon] = 0.5;
            this->prm_[CaMK__CaMKo] = 0.05;
            this->prm_[CaMK__KmCaM] = 0.0015;
            this->prm_[CaMK__PKNa] = 0.01833;
            this->prm_[cell__R] = 8314.0;
            this->prm_[cell__T] = 310.0;
            this->prm_[CaMK__aCaMK] = 0.05;
            this->prm_[CaMK__bCaMK] = 0.00068;
            this->prm_[cell__cli] = 19.53;
            this->prm_[cell__clo] = 100.0;
            this->prm_[cell__ko] = 5.4;
            this->prm_[cell__nao] = 140.0;
            this->prm_[ICaL__PCa] = 0.0002;
            this->prm_[ICaL__vhalf_d] = 3.4;
            this->prm_[ICaL__zca] = 2.0;
            this->prm_[IKb__GKb] = 0.003;
            this->prm_[IKs__GKs] = 0.021;
            this->prm_[INaCa_i__KmCaAct] = 0.00015;
            this->prm_[INaCa_i__kasymm] = 12.5;
            this->prm_[INaCa_i__kcaoff] = 5000.0;
            this->prm_[INaCa_i__kcaon] = 1500000.0;
            this->prm_[INaCa_i__kna1] = 15.0;
            this->prm_[INaCa_i__kna2] = 5.0;
            this->prm_[INaCa_i__kna3] = 88.12;
            this->prm_[INaCa_i__qca] = 0.167;
            this->prm_[INaCa_i__qna] = 0.5224;
            this->prm_[INaCa_i__wca] = 60000.0;
            this->prm_[INaCa_i__wna] = 60000.0;
            this->prm_[INaCa_i__wnaca] = 5000.0;
            this->prm_[INaCa_i__zna] = 1.0;
            this->prm_[INaK__H] = 1e-07;
            this->prm_[INaK__Khp] = 1.698e-07;
            this->prm_[INaK__Kki] = 0.5;
            this->prm_[INaK__Kko] = 0.3582;
            this->prm_[INaK__Kmgatp] = 1.698e-07;
            this->prm_[INaK__Knai0] = 9.073;
            this->prm_[INaK__Knao0] = 27.78;
            this->prm_[INaK__Knap] = 224.0;
            this->prm_[INaK__Kxkur] = 292.0;
            this->prm_[INaK__MgADP] = 0.05;
            this->prm_[INaK__MgATP] = 9.8;
            this->prm_[INaK__eP] = 4.2;
            this->prm_[INaK__k1m] = 182.4;
            this->prm_[INaK__k1p] = 949.5;
            this->prm_[INaK__k2m] = 39.4;
            this->prm_[INaK__k2p] = 687.2;
            this->prm_[INaK__k3m] = 79300.0;
            this->prm_[INaK__k3p2] = 1899.0;
            this->prm_[INaK__k4m] = 40.0;
            this->prm_[INaK__k4p2] = 639.0;
            this->prm_[INaK__zk] = 1.0;
            this->prm_[INaL__GNaL] = 0.0075;
            this->prm_[INaL__tau_hl] = 600.0;
            this->prm_[INab__PNab] = 3.75e-10;
            this->prm_[I_Na__GNa] = 14.838;
            this->prm_[IpCa__GpCa] = 0.0575;
            this->prm_[SR_uptake__BSLmax] = 1.124;
            this->prm_[SR_uptake__BSRmax] = 0.047;
            this->prm_[SR_uptake__KmBSL] = 0.0087;
            this->prm_[SR_uptake__KmBSR] = 0.00087;
            this->prm_[SR_uptake__cmdnmax] = 0.05;
            this->prm_[SR_uptake__csqnmax] = 10.0;
            this->prm_[SR_uptake__kmcmdn] = 0.00238;
            this->prm_[SR_uptake__kmcsqn] = 0.8;
            this->prm_[SR_uptake__kmtrpn] = 0.0005;
            this->prm_[SR_uptake__trpnmax] = 0.07;
            this->prm_[cell__L] = 0.017;
            this->prm_[cell__pi] = 3.14;
            this->prm_[cell__rad] = 0.0011;
            this->prm_[cell__vmyo1frac] = 0.4;
            this->prm_[stimulus__duration] = 0.5;
            this->prm_[stimulus__offset] = 20.0;
            this->prm_[stimulus__period] = 500.0;
            this->prm_[CaMK__ECl] =  (( this->prm_[cell__R]*this->prm_[cell__T])/this->prm_[cell__F])*std::log(this->prm_[cell__cli]/this->prm_[cell__clo]);
            this->prm_[ICaL__PCaK] =  0.000357400*this->prm_[ICaL__PCa];
            this->prm_[ICaL__PCaNa] =  0.00125000*this->prm_[ICaL__PCa];
            this->prm_[ICaL__vhalff] =  - 22.9000;
            this->prm_[ICab__PCab] =  2.50000e-08;
            this->prm_[IK1__GK1] =  0.75*0.200000;
            this->prm_[IKr__GKr] =  0.5*0.0150000;
            this->prm_[INaK__delta] =  - 0.155000;
            this->prm_[INaCa_i__Gncx] =  0.000800000;
            this->prm_[ITo__Gto] =  0.4*0.5;
            this->prm_[stimulus__stimulus_amplitude] =  - 80.0000;
            this->prm_[INaCa_i__h10] = (this->prm_[INaCa_i__kasymm]+1.00000)+ (this->prm_[cell__nao]/this->prm_[INaCa_i__kna1])*(1.00000+this->prm_[cell__nao]/this->prm_[INaCa_i__kna2]);
            this->prm_[INaCa_i__h11] = ( this->prm_[cell__nao]*this->prm_[cell__nao])/( ( this->prm_[INaCa_i__h10]*this->prm_[INaCa_i__kna1])*this->prm_[INaCa_i__kna2]);
            this->prm_[INaCa_i__h12] = 1/this->prm_[INaCa_i__h10];
            this->prm_[INaCa_i__k2] = this->prm_[INaCa_i__kcaoff];
            this->prm_[INaCa_i__k5] = this->prm_[INaCa_i__kcaoff];
            this->prm_[INaCa_ss__h101] = (this->prm_[INaCa_i__kasymm]+1.00000)+ (this->prm_[cell__nao]/this->prm_[INaCa_i__kna1])*(1.00000+this->prm_[cell__nao]/this->prm_[INaCa_i__kna2]);
            this->prm_[INaCa_ss__h1111] = ( this->prm_[cell__nao]*this->prm_[cell__nao])/( ( this->prm_[INaCa_ss__h101]*this->prm_[INaCa_i__kna1])*this->prm_[INaCa_i__kna2]);
            this->prm_[INaCa_ss__h121] = 1./this->prm_[INaCa_ss__h101];
            this->prm_[INaCa_ss__k21] = this->prm_[INaCa_i__kcaoff];
            this->prm_[INaCa_ss__k51] = this->prm_[INaCa_i__kcaoff];
            this->prm_[INaK__Pnak] =  30.0000;
            this->prm_[INaK__a2] = this->prm_[INaK__k2p];
            this->prm_[INaK__a4] = (( this->prm_[INaK__k4p2]*this->prm_[INaK__MgATP])/this->prm_[INaK__Kmgatp])/(1.00000+this->prm_[INaK__MgATP]/this->prm_[INaK__Kmgatp]);
            this->prm_[INaK__b1] =  this->prm_[INaK__k1m]*this->prm_[INaK__MgADP];
            this->prm_[cell__Ageo] =  ( ( 2*this->prm_[cell__pi])*this->prm_[cell__rad])*this->prm_[cell__rad] + ( ( 2*this->prm_[cell__pi])*this->prm_[cell__rad])*this->prm_[cell__L];
            this->prm_[cell__Acap] =  2*this->prm_[cell__Ageo];
            this->prm_[cell__cao] =  1.80000;
            this->prm_[INaCa_i__k1] =  ( this->prm_[INaCa_i__h12]*this->prm_[cell__cao])*this->prm_[INaCa_i__kcaon];
            this->prm_[INaCa_ss__k11] =  ( this->prm_[INaCa_ss__h121]*this->prm_[cell__cao])*this->prm_[INaCa_i__kcaon];
            this->prm_[cell__vcell] =  ( ( ( 1000*this->prm_[cell__pi])*this->prm_[cell__rad])*this->prm_[cell__rad])*this->prm_[cell__L];
            this->prm_[cell__vjsr] = ( 0.0048*this->prm_[cell__vcell])/2.00000;
            this->prm_[cell__vcsr] = this->prm_[cell__vjsr];
            this->prm_[cell__vmyo] =  0.68*this->prm_[cell__vcell];
            this->prm_[cell__vmyo1] =  this->prm_[cell__vmyo]*this->prm_[cell__vmyo1frac];
            this->prm_[cell__vmyo2] = this->prm_[cell__vmyo] - this->prm_[cell__vmyo1];
            this->prm_[cell__vnsr] =  0.0552*this->prm_[cell__vcell];
            this->prm_[cell__vnsr1] = this->prm_[cell__vnsr]/2.00000;
            this->prm_[cell__vnsr2] = this->prm_[cell__vnsr]/2.00000;
            this->prm_[cell__vss] = ( 0.02*this->prm_[cell__vcell])/2.00000;

            break;
    
        default:
            throw std::invalid_argument(Logger::Error("Could not initialize Gaur2021 ap model. Expected: CellType::ventricular"));
            break;
    }

}


void Gaur2021::Compute(double v_new, double dt, double stim_current)
{
    using namespace Gaur21Var;
    using namespace Gaur21Prm;
    using namespace Gaur21Cur;

    double CICR__diff_A = 1.;
    if (this->var_[CICR__tjsrol] <5.00000) {
      CICR__diff_A = - this->var_[CICR__A]/0.100000;
    }

    double CICR__diff_tjsrol = 1.;
    if (this->var_[ionic_concentrations__cajsr] > this->prm_[CICR__SOICR] && this->var_[CICR__A]>45.0000) {
        CICR__diff_tjsrol = - this->var_[CICR__tjsrol]/0.00100000;
    }
    
    double INaL__hl_inf = 1./(1.00000+std::exp((v_new+91.0000)/6.10000));

    double ICaL__d_inf = 1.00000/(1.00000+std::exp( - (v_new - this->prm_[ICaL__vhalf_d])/6.20000));
    double ICaL__tau_d = 0.600000+1.00000/(std::exp(  - 0.0500000*(v_new+6.00000))+std::exp( 0.0900000*(v_new+14.0000)));

    double IKr__tau_xr = 12.9800+1.00000/( 0.365200*std::exp((v_new - 31.6600)/3.86900)+ 4.12300e-05*std::exp( - (v_new - 47.7800)/20.3800));
    double IKr__xr_inf = 1.00000/(1.00000+std::exp( - (v_new+56.8000)/17.8000));

    double ITo__alpha_aa = 0.0250000/(1.00000+std::exp((v_new+58.0000)/5.00000));
    double ITo__beta_aa = 1.00000/( 5.00000*(1.00000+std::exp((v_new+19.0000)/ - 9.00000)));
    double ITo__aa_inf = ITo__alpha_aa / (ITo__alpha_aa + ITo__beta_aa);
    double ITo__tau_aa = 1.00000 / (ITo__alpha_aa + ITo__beta_aa);
    
    double ICaL__f_inf = 1.00000/(1.00000+std::exp((v_new - this->prm_[ICaL__vhalff])/4.90000))+0.350000/(1.00000+std::exp((45.0000 - v_new)/20.0000));
    double ICaL__ff_inf = ICaL__f_inf;
    double ICaL__tau_ff = 7.00000+1.00000/( 0.00450000*std::exp( - (v_new+20.0000)/10.0000)+ 0.00450000*std::exp((v_new+20.0000)/10.0000));

    double ICaL__fs_inf = ICaL__f_inf;
    double ICaL__tau_fs = 70.0000+1.00000/( 3.50000e-05*std::exp( - (v_new+5.00000)/4.00000)+ 3.50000e-05*std::exp((v_new+5.00000)/6.00000));
    
    double IKs__tau_xs1 = 300.000+1.00000/( 1.00000e-06*std::exp((v_new+50.0000)/20.0000)+ 0.0400000*std::exp( - (v_new+50.0000)/20.0000));
    double IKs__xs1_inf = 1.00000/(1.00000+std::exp( - (v_new - 25.1000)/37.1000));
    
    double IKs__tau_xs2 = 1.00000/( 0.0100000*std::exp((v_new - 50.0000)/20.0000)+ 0.0193000*std::exp( - (v_new+66.5400)/31.0000));
    double IKs__xs2_inf = IKs__xs1_inf;
    
    double INaL__aml =  0.320000 * (- (v_new+47.1300)/(std::exp(  - 0.100000*(v_new+47.1300)) - 1.00000));
    if (v_new == - 47.1300) {
        INaL__aml = 3.2;
    }

    double INaL__bml =  0.0800000*std::exp( - v_new/11.0000);
    double INaL__ml_inf = INaL__aml/(INaL__aml+INaL__bml);
    double INaL__tau_ml = 1.00000/(INaL__aml+INaL__bml);
    
    double I_Na__m_inf = 1.00000/( (1.00000+std::exp(( - 56.8600 - v_new)/9.03000))*(1.00000+std::exp(( - 56.8600 - v_new)/9.03000)));
    double I_Na__aa_m = 1.00000/(1.00000+std::exp(( - 60.0000 - v_new)/5.00000));
    double I_Na__bb_m = 0.100000/(1.00000+std::exp((v_new+35.0000)/5.00000))+0.1/(1.00000+std::exp((v_new - 50.0000)/200.000));
    double I_Na__tau_m =  I_Na__aa_m*I_Na__bb_m;
    
    double CaMK__CaMKb = ( this->prm_[CaMK__CaMKo]*(1.00000 - this->var_[CaMK__CaMKt]))/(1.00000+this->prm_[CaMK__KmCaM]/this->var_[ionic_concentrations__cass]);
    double CaMK__diff_CaMKt =  ( this->prm_[CaMK__aCaMK]*CaMK__CaMKb)*(CaMK__CaMKb+this->var_[CaMK__CaMKt]) -  this->prm_[CaMK__bCaMK]*this->var_[CaMK__CaMKt];
    
    double I_Na__h_inf = 1.00000/( (1.00000+std::exp((v_new+71.5500)/7.43000))*(1.00000+std::exp((v_new+71.5500)/7.43000)));
    double I_Na__aa_h = 0.0570000*std::exp( - (v_new+80.0000)/6.80000);
    if (v_new >= - 40.0000) {
        I_Na__aa_h = 0.00000;
    }
    
    double I_Na__bb_h = 2.70000*std::exp( 0.0790000*v_new)+ 310000*std::exp( 0.348500*v_new);
    if (v_new>= - 40.0000) {
        I_Na__bb_h = 0.770000/( 0.130000*(1.00000+std::exp( - (v_new+10.6600)/11.1000)));
    }

    double I_Na__tau_h = 1.00000/(I_Na__aa_h+I_Na__bb_h);    
    double I_Na__j_inf = I_Na__h_inf;

    double I_Na__aa_j = ((- 25428*std::exp( 0.2444*v_new) -  6.94800e-06*std::exp(- 0.0439100*v_new))*(v_new+37.7800))/(1.00000+std::exp(0.311*(v_new+79.2300)));
    if (v_new >= - 40.0000) {
        I_Na__aa_j = 0.;
    }
    
    double I_Na__bb_j = ( 0.0242400*std::exp(- 0.01052*v_new))/(1.00000+std::exp(- 0.137800*(v_new+40.14)));
    if (v_new >= - 40.0000) {
        I_Na__bb_j = (0.6*std::exp( 0.057*v_new))/(1.00000+std::exp(- 0.1*(v_new+32.0000)));
    } 
    
    double I_Na__tau_j = 1.00000/(I_Na__aa_j+I_Na__bb_j);
    
    double CaMK__vfrt = ( v_new*this->prm_[cell__F])/( this->prm_[cell__R]*this->prm_[cell__T]);

    double ICaL__PhiCaL_temp = ( 2.00000*CaMK__vfrt)/(std::exp( 2.00000*CaMK__vfrt) - 1.00000);
    if (v_new == 0.00000) {
        ICaL__PhiCaL_temp = 1.00000; 
    }
    double ICaL__PhiCaL =  ( ( 2.00000*this->prm_[cell__F])*( this->var_[ionic_concentrations__cass]*std::exp( 2.00000*(( v_new*this->prm_[cell__F])/( this->prm_[cell__R]*this->prm_[cell__T]))) -  0.341000*this->prm_[cell__cao]))*ICaL__PhiCaL_temp;
    double ICaL__f =  this->var_[ICaL__ff]*this->var_[ICaL__fs];
    this->cur_[ICaL] =  ( ( ( this->prm_[ICaL__PCa]*ICaL__PhiCaL)*this->var_[ICaL__d])*ICaL__f)*this->var_[ICaL__fca];
    this->cur_[ICaL] *= (1.0-this->block_coeff_[ICaL]); 
    
    double ICaL__fca_inf = (0.300000/(1.00000 - this->cur_[ICaL]/0.0500000) +0.550000/(1.00000+this->var_[ionic_concentrations__cass]/0.00300000))+0.150000;
    if (this->cur_[ICaL] > 0.) {
        ICaL__fca_inf = (0.300000/(1.00000)+0.550000/(1.00000+this->var_[ionic_concentrations__cass]/0.00300000))+0.150000;
    }
    double CaMK__CaMKa = CaMK__CaMKb+this->var_[CaMK__CaMKt];
    double CaMK__CaMK_f = 1.00000/(1.00000+this->prm_[CaMK__KmCaMK]/CaMK__CaMKa);
    double ICaL__tau_fca = (10.0000*CaMK__CaMK_f+0.500000)+1.00000/(1.00000+this->var_[ionic_concentrations__cass]/0.00300000);
    
    double CaMK__EK =  (( this->prm_[cell__R]*this->prm_[cell__T])/this->prm_[cell__F])*std::log(this->prm_[cell__ko]/this->var_[ionic_concentrations__ki]);
    double IK1__rk1 = 1.00000/(1.00000+std::exp(((v_new+79.3000) -  2.60000*this->prm_[cell__ko])/19.6000));
    this->cur_[IK1] =  ( ( this->prm_[IK1__GK1]*std::pow((this->prm_[cell__ko]/5.40000), 0.5))*IK1__rk1)*(v_new - CaMK__EK);
    this->cur_[IK1] *= (1.0-this->block_coeff_[IK1]); 

    double IKb__xkb = 1.00000/(1.00000+std::exp( - (v_new - 14.4800)/18.3400));
    this->cur_[IKb] =  ( this->prm_[IKb__GKb]*IKb__xkb)*(v_new - CaMK__EK);
    this->cur_[IKb] *= (1.0-this->block_coeff_[IKb]); 

    double IKr__rkr = 1.00000/(1.00000+std::exp((v_new+22.0000)/15.0000));
    this->cur_[IKr] = (((this->prm_[IKr__GKr]*std::pow((this->prm_[cell__ko]/5.40000), 0.5))*this->var_[IKr__xr])*IKr__rkr)*(v_new - CaMK__EK);
    this->cur_[IKr] *= (1.0-this->block_coeff_[IKr]);
    
    double CaMK__EKs =  (( this->prm_[cell__R]*this->prm_[cell__T])/this->prm_[cell__F])*std::log((this->prm_[cell__ko]+ this->prm_[CaMK__PKNa]*this->prm_[cell__nao])/(this->var_[ionic_concentrations__ki]+ this->prm_[CaMK__PKNa]*this->var_[ionic_concentrations__nai]));
    double IKs__KsCa = 1.00000+0.600000/(1.00000+std::pow(3.80000e-05/this->var_[ionic_concentrations__cai], 1.40000));
    this->cur_[IKs] =  ( ( (this->prm_[IKs__GKs]*IKs__KsCa)*this->var_[IKs__xs1])*this->var_[IKs__xs2])*(v_new - CaMK__EKs);
    this->cur_[IKs] *= (1.0-this->block_coeff_[IKs]);

    double INaK__Knai =  this->prm_[INaK__Knai0]*std::exp((( this->prm_[INaK__delta]*v_new)*this->prm_[cell__F])/((3*this->prm_[cell__R])*this->prm_[cell__T]));
    double INaK__a1 = ( this->prm_[INaK__k1p]*(( (this->var_[ionic_concentrations__nai]/INaK__Knai)*(this->var_[ionic_concentrations__nai]/INaK__Knai))*(this->var_[ionic_concentrations__nai]/INaK__Knai)))/(( ( (1.00000+this->var_[ionic_concentrations__nai]/INaK__Knai)*(1.00000+this->var_[ionic_concentrations__nai]/INaK__Knai))*(1.00000+this->var_[ionic_concentrations__nai]/INaK__Knai)+ (1.00000+this->var_[ionic_concentrations__ki]/this->prm_[INaK__Kki])*(1.00000+this->var_[ionic_concentrations__ki]/this->prm_[INaK__Kki])) - 1.00000);
    double INaK__Knao =  this->prm_[INaK__Knao0]*std::exp(( ( (1.00000 - this->prm_[INaK__delta])*v_new)*this->prm_[cell__F])/( (3.00000*this->prm_[cell__R])*this->prm_[cell__T]));
    double INaK__b2 = ( this->prm_[INaK__k2m]*( ( (this->prm_[cell__nao]/INaK__Knao)*(this->prm_[cell__nao]/INaK__Knao))*(this->prm_[cell__nao]/INaK__Knao)))/(( ( (1.00000+this->prm_[cell__nao]/INaK__Knao)*(1.00000+this->prm_[cell__nao]/INaK__Knao))*(1.00000+this->prm_[cell__nao]/INaK__Knao)+ (1.00000+this->prm_[cell__ko]/this->prm_[INaK__Kko])*(1.00000+this->prm_[cell__ko]/this->prm_[INaK__Kko])) - 1.00000);
    double INaK__P = this->prm_[INaK__eP]/(((1.00000+this->prm_[INaK__H]/this->prm_[INaK__Khp])+this->var_[ionic_concentrations__nai]/this->prm_[INaK__Knap])+this->var_[ionic_concentrations__ki]/this->prm_[INaK__Kxkur]);
    double INaK__b3 = ( ( this->prm_[INaK__k3m]*INaK__P)*this->prm_[INaK__H])/(1.00000+this->prm_[INaK__MgATP]/this->prm_[INaK__Kmgatp]);
    double INaK__b4 = ( this->prm_[INaK__k4m]*( (this->var_[ionic_concentrations__ki]/this->prm_[INaK__Kki])*(this->var_[ionic_concentrations__ki]/this->prm_[INaK__Kki])))/(( ( (1.00000+this->var_[ionic_concentrations__nai]/INaK__Knai)*(1.00000+this->var_[ionic_concentrations__nai]/INaK__Knai))*(1.00000+this->var_[ionic_concentrations__nai]/INaK__Knai)+ (1.00000+this->var_[ionic_concentrations__ki]/this->prm_[INaK__Kki])*(1.00000+this->var_[ionic_concentrations__ki]/this->prm_[INaK__Kki])) - 1.00000);
    double INaK__x12 = (( ( this->prm_[INaK__a4]*INaK__a1)*this->prm_[INaK__a2] + ( INaK__b2*INaK__b4)*INaK__b3)+ ( this->prm_[INaK__a2]*INaK__b4)*INaK__b3)+ ( INaK__b3*INaK__a1)*this->prm_[INaK__a2];
    double INaK__a3 = ( this->prm_[INaK__k3p2]*( (this->prm_[cell__ko]/this->prm_[INaK__Kko])*(this->prm_[cell__ko]/this->prm_[INaK__Kko])))/(( ( (1.00000+this->prm_[cell__nao]/INaK__Knao)*(1.00000+this->prm_[cell__nao]/INaK__Knao))*(1.00000+this->prm_[cell__nao]/INaK__Knao)+ (1.00000+this->prm_[cell__ko]/this->prm_[INaK__Kko])*(1.00000+this->prm_[cell__ko]/this->prm_[INaK__Kko])) - 1.00000);
    double INaK__x22 = (( ( INaK__b2*this->prm_[INaK__b1])*INaK__b4+ ( INaK__a1*this->prm_[INaK__a2])*INaK__a3)+ ( INaK__a3*this->prm_[INaK__b1])*INaK__b4)+ ( this->prm_[INaK__a2]*INaK__a3)*INaK__b4;
    double INaK__x32 = (( ( this->prm_[INaK__a2]*INaK__a3)*this->prm_[INaK__a4]+ ( INaK__b3*INaK__b2)*this->prm_[INaK__b1])+ ( INaK__b2*this->prm_[INaK__b1])*this->prm_[INaK__a4])+ ( INaK__a3*this->prm_[INaK__a4])*this->prm_[INaK__b1];
    double INaK__x42 = (( ( INaK__b4*INaK__b3)*INaK__b2+ ( INaK__a3*this->prm_[INaK__a4])*INaK__a1)+ ( INaK__b2*this->prm_[INaK__a4])*INaK__a1)+ ( INaK__b3*INaK__b2)*INaK__a1;
    double INaK__E32 = INaK__x32/(((INaK__x12+INaK__x22)+INaK__x32)+INaK__x42);
    double INaK__E42 = INaK__x42/(((INaK__x12+INaK__x22)+INaK__x32)+INaK__x42);
    double INaK__JnakK =  2.00000*( INaK__E42*this->prm_[INaK__b1] -  INaK__E32*INaK__a1);
    double INaK__E12 = INaK__x12/(((INaK__x12+INaK__x22)+INaK__x32)+INaK__x42);
    double INaK__E22 = INaK__x22/(((INaK__x12+INaK__x22)+INaK__x32)+INaK__x42);
    double INaK__JnakNa =  3.00000*( INaK__E12*INaK__a3 -  INaK__E22*INaK__b3);
    this->cur_[INaK] =  this->prm_[INaK__Pnak]*( this->prm_[INaCa_i__zna]*INaK__JnakNa+ this->prm_[INaK__zk]*INaK__JnakK);
    this->cur_[INaK] *= (1.0-this->block_coeff_[INaK]);

    double diffusion__JdiffK = (this->var_[ionic_concentrations__kss] - this->var_[ionic_concentrations__ki])/2.00000;
    double ionic_concentrations__diff_ki = (  - ((((this->cur_[IKr]+this->cur_[IKs])+this->cur_[IK1])+this->cur_[IKb]) -  2.00000*this->cur_[INaK])*this->prm_[cell__Acap])/( ( 2.00000*this->prm_[cell__F])*this->prm_[cell__vmyo])+( diffusion__JdiffK*this->prm_[cell__vss])/this->prm_[cell__vmyo];
    
    double ICaL__PhiCaK_temp = CaMK__vfrt/(exp(CaMK__vfrt) - 1.00000);
    if (v_new==0.00000) {
        ICaL__PhiCaK_temp = 1.00000;
    }
    double ICaL__PhiCaK =  ( this->prm_[cell__F]*( ( 0.750000*this->var_[ionic_concentrations__kss])*std::exp(CaMK__vfrt) -  0.750000*this->prm_[cell__ko]))*ICaL__PhiCaK_temp;
    this->cur_[ICaK] =  ( ( ( this->prm_[ICaL__PCaK]*ICaL__PhiCaK)*this->var_[ICaL__d])*ICaL__f)*this->var_[ICaL__fca];
    
    double ionic_concentrations__diff_kss = (-this->cur_[ICaK]*this->prm_[cell__Acap])/( ( 2.00000*this->prm_[cell__F])*this->prm_[cell__vss]) - diffusion__JdiffK;
    
    double INaCa_i__hna = std::exp(( ( this->prm_[INaCa_i__qna]*v_new)*this->prm_[cell__F])/( this->prm_[cell__R]*this->prm_[cell__T]));
    double INaCa_i__h7 = 1.00000+ (this->prm_[cell__nao]/this->prm_[INaCa_i__kna3])*(1.00000+1.00000/INaCa_i__hna);
    double INaCa_i__h9 = 1.00000/INaCa_i__h7;
    double INaCa_i__k3p =  INaCa_i__h9*this->prm_[INaCa_i__wca];
    double INaCa_i__h8 = this->prm_[cell__nao]/( ( this->prm_[INaCa_i__kna3]*INaCa_i__hna)*INaCa_i__h7);
    double INaCa_i__k3pp =  INaCa_i__h8*this->prm_[INaCa_i__wnaca];
    double INaCa_i__k3 = INaCa_i__k3p+INaCa_i__k3pp;
    double INaCa_i__h1 = 1.00000+ (this->var_[ionic_concentrations__nai]/this->prm_[INaCa_i__kna3])*(1.00000+INaCa_i__hna);
    double INaCa_i__h3 = 1.00000/INaCa_i__h1;
    double INaCa_i__hca = std::exp(( ( this->prm_[INaCa_i__qca]*v_new)*this->prm_[cell__F])/( this->prm_[cell__R]*this->prm_[cell__T]));
    double INaCa_i__k4p = ( INaCa_i__h3*this->prm_[INaCa_i__wca])/INaCa_i__hca;
    double INaCa_i__h2 = ( this->var_[ionic_concentrations__nai]*INaCa_i__hna)/( this->prm_[INaCa_i__kna3]*INaCa_i__h1);
    double INaCa_i__k4pp =  INaCa_i__h2*this->prm_[INaCa_i__wnaca];
    double INaCa_i__k4 = INaCa_i__k4p+INaCa_i__k4pp;
    double INaCa_i__h4 = 1.00000+ (this->var_[ionic_concentrations__nai]/this->prm_[INaCa_i__kna1])*(1.00000+this->var_[ionic_concentrations__nai]/this->prm_[INaCa_i__kna2]);
    double INaCa_i__h6 = 1.00000/INaCa_i__h4;
    double INaCa_i__k6 =  ( INaCa_i__h6*this->var_[ionic_concentrations__cai])*this->prm_[INaCa_i__kcaon];
    double INaCa_i__h5 = ( this->var_[ionic_concentrations__nai]*this->var_[ionic_concentrations__nai])/( ( INaCa_i__h4*this->prm_[INaCa_i__kna1])*this->prm_[INaCa_i__kna2]);
    double INaCa_i__k7 =  ( INaCa_i__h5*INaCa_i__h2)*this->prm_[INaCa_i__wna];
    double INaCa_i__x1 =  ( this->prm_[INaCa_i__k2]*INaCa_i__k4)*(INaCa_i__k7+INaCa_i__k6)+ ( this->prm_[INaCa_i__k5]*INaCa_i__k7)*(this->prm_[INaCa_i__k2]+INaCa_i__k3);
    double INaCa_i__k8 =  ( INaCa_i__h8*this->prm_[INaCa_i__h11])*this->prm_[INaCa_i__wna];
    double INaCa_i__x2 =  ( this->prm_[INaCa_i__k1]*INaCa_i__k7)*(INaCa_i__k4+this->prm_[INaCa_i__k5])+ ( INaCa_i__k4*INaCa_i__k6)*(this->prm_[INaCa_i__k1]+INaCa_i__k8);
    double INaCa_i__x3 =  ( this->prm_[INaCa_i__k1]*INaCa_i__k3)*(INaCa_i__k7+INaCa_i__k6)+ ( INaCa_i__k8*INaCa_i__k6)*(this->prm_[INaCa_i__k2]+INaCa_i__k3);
    double INaCa_i__x4 =  ( this->prm_[INaCa_i__k2]*INaCa_i__k8)*(INaCa_i__k4+this->prm_[INaCa_i__k5])+ ( INaCa_i__k3*this->prm_[INaCa_i__k5])*(this->prm_[INaCa_i__k1]+INaCa_i__k8);
    double INaCa_i__E1 = INaCa_i__x1/(((INaCa_i__x1+INaCa_i__x2)+INaCa_i__x3)+INaCa_i__x4);
    double INaCa_i__E2 = INaCa_i__x2/(((INaCa_i__x1+INaCa_i__x2)+INaCa_i__x3)+INaCa_i__x4);
    double INaCa_i__JncxCa =  INaCa_i__E2*this->prm_[INaCa_i__k2] -  INaCa_i__E1*this->prm_[INaCa_i__k1];
    double INaCa_i__E3 = INaCa_i__x3/(((INaCa_i__x1+INaCa_i__x2)+INaCa_i__x3)+INaCa_i__x4);
    double INaCa_i__E4 = INaCa_i__x4/(((INaCa_i__x1+INaCa_i__x2)+INaCa_i__x3)+INaCa_i__x4);
    double INaCa_i__JncxNa = ( 3.00000*( INaCa_i__E4*INaCa_i__k7 -  INaCa_i__E1*INaCa_i__k8)+ INaCa_i__E3*INaCa_i__k4pp) -  INaCa_i__E2*INaCa_i__k3pp;
    double INaCa_i__allo = 1.00000/(1.00000+ (this->prm_[INaCa_i__KmCaAct]/this->var_[ionic_concentrations__cai])*(this->prm_[INaCa_i__KmCaAct]/this->var_[ionic_concentrations__cai]));
    this->cur_[INaCa_i] =  ( ( 0.800000*this->prm_[INaCa_i__Gncx])*INaCa_i__allo)*( this->prm_[INaCa_i__zna]*INaCa_i__JncxNa+ this->prm_[ICaL__zca]*INaCa_i__JncxCa);
    
    double CaMK__ENa =  (( this->prm_[cell__R]*this->prm_[cell__T])/this->prm_[cell__F])*std::log(this->prm_[cell__nao]/this->var_[ionic_concentrations__nai]);
    this->cur_[INaL] =  ( ( ( ( this->prm_[INaL__GNaL]*this->var_[INaL__ml])*this->var_[INaL__ml])*this->var_[INaL__ml])*this->var_[INaL__hl])*(v_new - CaMK__ENa);
    this->cur_[INaL] *= (1.0-this->block_coeff_[INaL]);

    double INab_temp = CaMK__vfrt/(std::exp(CaMK__vfrt) - 1.00000);
    if (v_new == 0.00000) {
        INab_temp = 1.00000;
    } 
    this->cur_[INab] =  ( ( this->prm_[INab__PNab]*this->prm_[cell__F])*( this->var_[ionic_concentrations__nai]*std::exp(CaMK__vfrt) - this->prm_[cell__nao])) * INab_temp;
    this->cur_[INab] *= (1.0-this->block_coeff_[INab]);

    this->cur_[INa] =  ( ( ( ( (this->prm_[I_Na__GNa]*this->var_[I_Na__m])*this->var_[I_Na__m])*this->var_[I_Na__m])*this->var_[I_Na__h])*this->var_[I_Na__j])*(v_new - CaMK__ENa);
    this->cur_[INa] *= (1.0-this->block_coeff_[INa]);

    double diffusion__JdiffNa = (this->var_[ionic_concentrations__nass] - this->var_[ionic_concentrations__nai])/2.00000;
    double ionic_concentrations__diff_nai = (  - ((((this->cur_[INa]+this->cur_[INaL])+ 3.00000*this->cur_[INaCa_i])+ 3.00000*this->cur_[INaK])+this->cur_[INab])*this->prm_[cell__Acap])/( ( 2.00000*this->prm_[cell__F])*this->prm_[cell__vmyo])+( diffusion__JdiffNa*this->prm_[cell__vss])/this->prm_[cell__vmyo];
    
    double INaCa_ss__h71 = 1.00000+ (this->prm_[cell__nao]/this->prm_[INaCa_i__kna3])*(1.00000+1.00000/INaCa_i__hna);
    double INaCa_ss__h91 = 1.00000/INaCa_ss__h71;
    double INaCa_ss__k3p1 =  INaCa_ss__h91*this->prm_[INaCa_i__wca];
    double INaCa_ss__h81 = this->prm_[cell__nao]/( ( this->prm_[INaCa_i__kna3]*INaCa_i__hna)*INaCa_ss__h71);
    double INaCa_ss__k3pp1 =  INaCa_ss__h81*this->prm_[INaCa_i__wnaca];
    double INaCa_ss__k31 = INaCa_ss__k3p1+INaCa_ss__k3pp1;
    double INaCa_ss__h111 = 1.00000+ (this->var_[ionic_concentrations__nass]/this->prm_[INaCa_i__kna3])*(1.00000+INaCa_i__hna);
    double INaCa_ss__h31 = 1.00000/INaCa_ss__h111;
    double INaCa_ss__k4p1 = ( INaCa_ss__h31*this->prm_[INaCa_i__wca])/INaCa_i__hca;
    double INaCa_ss__h21 = ( this->var_[ionic_concentrations__nass]*INaCa_i__hna)/( this->prm_[INaCa_i__kna3]*INaCa_ss__h111);
    double INaCa_ss__k4pp1 =  INaCa_ss__h21*this->prm_[INaCa_i__wnaca];
    double INaCa_ss__k41 = INaCa_ss__k4p1+INaCa_ss__k4pp1;
    double INaCa_ss__h41 = 1.00000+ (this->var_[ionic_concentrations__nass]/this->prm_[INaCa_i__kna1])*(1.00000+this->var_[ionic_concentrations__nass]/this->prm_[INaCa_i__kna2]);
    double INaCa_ss__h61 = 1.00000/INaCa_ss__h41;
    double INaCa_ss__k61 =  ( INaCa_ss__h61*this->var_[ionic_concentrations__cass])*this->prm_[INaCa_i__kcaon];
    double INaCa_ss__h51 = ( this->var_[ionic_concentrations__nass]*this->var_[ionic_concentrations__nass])/( ( INaCa_ss__h41*this->prm_[INaCa_i__kna1])*this->prm_[INaCa_i__kna2]);
    double INaCa_ss__k71 =  ( INaCa_ss__h51*INaCa_ss__h21)*this->prm_[INaCa_i__wna];
    double INaCa_ss__x11 =  ( this->prm_[INaCa_ss__k21]*INaCa_ss__k41)*(INaCa_ss__k71+INaCa_ss__k61)+ ( this->prm_[INaCa_ss__k51]*INaCa_ss__k71)*(this->prm_[INaCa_ss__k21]+INaCa_ss__k31);
    double INaCa_ss__k81 =  ( INaCa_ss__h81*this->prm_[INaCa_ss__h1111])*this->prm_[INaCa_i__wna];
    double INaCa_ss__x21 =  ( this->prm_[INaCa_ss__k11]*INaCa_ss__k71)*(INaCa_ss__k41+this->prm_[INaCa_ss__k51])+ ( INaCa_ss__k41*INaCa_ss__k61)*(this->prm_[INaCa_ss__k11]+INaCa_ss__k81);
    double INaCa_ss__x31 =  ( this->prm_[INaCa_ss__k11]*INaCa_ss__k31)*(INaCa_ss__k71+INaCa_ss__k61)+ ( INaCa_ss__k81*INaCa_ss__k61)*(this->prm_[INaCa_ss__k21]+INaCa_ss__k31);
    double INaCa_ss__x41 =  ( this->prm_[INaCa_ss__k21]*INaCa_ss__k81)*(INaCa_ss__k41+this->prm_[INaCa_ss__k51])+ ( INaCa_ss__k31*this->prm_[INaCa_ss__k51])*(this->prm_[INaCa_ss__k11]+INaCa_ss__k81);
    double INaCa_ss__E11 = INaCa_ss__x11/(((INaCa_ss__x11+INaCa_ss__x21)+INaCa_ss__x31)+INaCa_ss__x41);
    double INaCa_ss__E21 = INaCa_ss__x21/(((INaCa_ss__x11+INaCa_ss__x21)+INaCa_ss__x31)+INaCa_ss__x41);
    double INaCa_ss__JncxCa1 =  INaCa_ss__E21*this->prm_[INaCa_ss__k21] -  INaCa_ss__E11*this->prm_[INaCa_ss__k11];
    double INaCa_ss__E31 = INaCa_ss__x31/(((INaCa_ss__x11+INaCa_ss__x21)+INaCa_ss__x31)+INaCa_ss__x41);
    double INaCa_ss__E41 = INaCa_ss__x41/(((INaCa_ss__x11+INaCa_ss__x21)+INaCa_ss__x31)+INaCa_ss__x41);
    double INaCa_ss__JncxNa1 = ( 3.00000*( INaCa_ss__E41*INaCa_ss__k71 -  INaCa_ss__E11*INaCa_ss__k81)+ INaCa_ss__E31*INaCa_ss__k4pp1) -  INaCa_ss__E21*INaCa_ss__k3pp1;
    double INaCa_ss__allo1 = 1.00000/(1.00000+ (this->prm_[INaCa_i__KmCaAct]/this->var_[ionic_concentrations__cass])*(this->prm_[INaCa_i__KmCaAct]/this->var_[ionic_concentrations__cass]));
    this->cur_[INaCa_ss] =  ( ( 0.200000*this->prm_[INaCa_i__Gncx])*INaCa_ss__allo1)*( this->prm_[INaCa_i__zna]*INaCa_ss__JncxNa1+ this->prm_[ICaL__zca]*INaCa_ss__JncxCa1);
    this->cur_[INaCa_ss] *= (1.0-this->block_coeff_[INaCa_ss]);

    double ICaL__PhiCaNa_temp = CaMK__vfrt/(std::exp(CaMK__vfrt) - 1.00000);
    if (v_new == 0.00000) {
        ICaL__PhiCaNa_temp = 1.0;
    }
    double ICaL__PhiCaNa =  ( this->prm_[cell__F]*( ( 0.750000*this->var_[ionic_concentrations__nass])*std::exp(CaMK__vfrt) -  0.750000*this->prm_[cell__nao]))*ICaL__PhiCaNa_temp;
    this->cur_[ICaNa] =  ( ( ( this->prm_[ICaL__PCaNa]*ICaL__PhiCaNa)*this->var_[ICaL__d])*ICaL__f)*this->var_[ICaL__fca];
    this->cur_[ICaNa] *= (1.0-this->block_coeff_[ICaNa]);

    double ionic_concentrations__diff_nass = (  - (this->cur_[ICaNa]+ 3.00000*this->cur_[INaCa_ss])*this->prm_[cell__Acap])/( ( 2.00000*this->prm_[cell__F])*this->prm_[cell__vss]) - diffusion__JdiffNa;
    double ICab_temp = (2.00000*CaMK__vfrt) / (std::exp(2.00000*CaMK__vfrt) - 1.00000);
    if (v_new == 0.00000) {
        ICab_temp = 1.00000;
    }
    this->cur_[ICab] =  ( ( ( this->prm_[ICab__PCab]*2.00000)*this->prm_[cell__F])*( this->var_[ionic_concentrations__cai]*std::exp( 2.00000*CaMK__vfrt) -  0.341000*this->prm_[cell__cao]))*ICab_temp;
    this->cur_[ICab] *= (1.0-this->block_coeff_[ICab]);

    double CICR__greljsrol =  ( this->prm_[CICR__grelbarjsrol]*(1.00000 - std::exp( - this->var_[CICR__tjsrol]/this->prm_[CICR__tauon])))*std::exp( - this->var_[CICR__tjsrol]/this->prm_[CICR__tauoff]);
    double CICR__Jrelol =  CICR__greljsrol*(this->var_[ionic_concentrations__cajsr] - this->var_[ionic_concentrations__cass]);
    double CICR__Jrel = this->var_[CICR__Jrel1]+CICR__Jrelol;
    double ITo__kito2 = 1.00000 - 1.00000/(1.00000+ (CICR__Jrel/0.400000)*(CICR__Jrel/0.400000));
    double ITo__rito2 = 1.00000/(1.00000+std::exp( - (v_new+10.0000)/5.00000));
    this->cur_[ITo] =  ( ( ( this->prm_[ITo__Gto]*this->var_[ITo__aa])*ITo__rito2)*ITo__kito2)*(v_new - this->prm_[CaMK__ECl]);
    this->cur_[INa] *= (1.0-this->block_coeff_[INa]);

    this->cur_[IpCa] = ( this->prm_[IpCa__GpCa]*this->var_[ionic_concentrations__cai])/(0.000500000+this->var_[ionic_concentrations__cai]);
    this->cur_[IpCa] *= (1.0-this->block_coeff_[IpCa]);

    // Compute the total Iion current.    
    this->cur_[Gaur21Cur::Iion] = ((((((((((((((this->cur_[INa]+this->cur_[INaL])+this->cur_[ICaL])+this->cur_[ICaNa])+this->cur_[ICaK])+this->cur_[IKr])+this->cur_[IKs])+this->cur_[IK1])+this->cur_[ITo])+this->cur_[INaCa_i])+this->cur_[INaCa_ss])+this->cur_[INaK])+this->cur_[INab])+this->cur_[IKb])+this->cur_[IpCa])+this->cur_[ICab];
    
    // Set the new value to the membrante potential time derivative.
    this->var_[dvdt] = - (this->cur_[Gaur21Cur::Iion] - stim_current);


    double SR_uptake__Jtr = (this->var_[ionic_concentrations__cansr] - this->var_[ionic_concentrations__cajsr])/100.000;
    double ionic_concentrations__Bcajsr = 1.00000/(1.00000+( this->prm_[SR_uptake__csqnmax]*this->prm_[SR_uptake__kmcsqn])/( (this->prm_[SR_uptake__kmcsqn]+this->var_[ionic_concentrations__cajsr])*(this->prm_[SR_uptake__kmcsqn]+this->var_[ionic_concentrations__cajsr])));
    double ionic_concentrations__diff_cajsr =  ionic_concentrations__Bcajsr*(SR_uptake__Jtr - CICR__Jrel);
    
    double SR_uptake__Jtr2 = (this->var_[ionic_concentrations__cansr] - this->var_[ionic_concentrations__cacsr])/100.000;
    double ionic_concentrations__Bcacsr = 1.00000/(1.00000+( this->prm_[SR_uptake__csqnmax]*this->prm_[SR_uptake__kmcsqn])/( (this->prm_[SR_uptake__kmcsqn]+this->var_[ionic_concentrations__cacsr])*(this->prm_[SR_uptake__kmcsqn]+this->var_[ionic_concentrations__cacsr])));
    double ionic_concentrations__diff_cacsr =  ionic_concentrations__Bcacsr*(SR_uptake__Jtr2 - this->var_[CICR__Jrel2]);
    
    double SR_uptake__Jleak = ( 0.00500000*this->var_[ionic_concentrations__cansr])/15.0000;
    double SR_uptake__Jupnp2 = ( ( 1.00000*0.00500000)*this->var_[ionic_concentrations__cai2])/(this->var_[ionic_concentrations__cai2]+0.00100000);
    double SR_uptake__Jupp2 = ( ( ( 1.00000*3.00000)*0.00500000)*this->var_[ionic_concentrations__cai2])/((this->var_[ionic_concentrations__cai2]+0.00100000) - 0.000200000);
    double SR_uptake__fJupp = 1.00000/(1.00000+this->prm_[CaMK__KmCaMK]/CaMK__CaMKa);
    double SR_uptake__Jup2 = ( (1.00000 - SR_uptake__fJupp)*SR_uptake__Jupnp2+ SR_uptake__fJupp*SR_uptake__Jupp2) - SR_uptake__Jleak;
    double SR_uptake__Jupnp = ( ( 1.00000*0.00500000)*this->var_[ionic_concentrations__cai])/(this->var_[ionic_concentrations__cai]+0.00100000);
    double SR_uptake__Jupp = ( ( ( 1.00000*3.00000)*0.00500000)*this->var_[ionic_concentrations__cai])/((this->var_[ionic_concentrations__cai]+0.00100000) - 0.000200000);
    double SR_uptake__Jup = ( (1.00000 - SR_uptake__fJupp)*SR_uptake__Jupnp+ SR_uptake__fJupp*SR_uptake__Jupp) - SR_uptake__Jleak;
    double ionic_concentrations__diff_cansr = (( SR_uptake__Jup*(this->prm_[cell__vnsr1]/this->prm_[cell__vnsr])+ SR_uptake__Jup2*(this->prm_[cell__vnsr2]/this->prm_[cell__vnsr])) - ( SR_uptake__Jtr*this->prm_[cell__vjsr])/this->prm_[cell__vnsr]) - ( SR_uptake__Jtr2*this->prm_[cell__vcsr])/this->prm_[cell__vnsr];
    
    double diffusion__Jdiff = (this->var_[ionic_concentrations__cass] - this->var_[ionic_concentrations__cai])/0.200000;
    double CICR__Rel1 = ( ( - this->cur_[ICaL]+ 2.00000*this->cur_[INaCa_ss])*(this->prm_[cell__Acap]/( ( 2.00000*this->prm_[cell__vss])*this->prm_[cell__F]))+ CICR__Jrel*(this->prm_[cell__vjsr]/this->prm_[cell__vss])) - diffusion__Jdiff;
    
    double CICR__Jrel1_inf = 0.00000;
    if (CICR__Rel1 > 0.00000) { 
        CICR__Jrel1_inf = ((( 1.00000*60.0000)*CICR__Rel1)*(1.00000+1.00000/(1.00000+std::pow(this->prm_[CaMK__KmCaMK]/CaMK__CaMKa, 8.00000)))) / (1.00000+std::pow(0.750000/this->var_[ionic_concentrations__cajsr], 8.00000));
    }
    double CICR__trel1factor = ( ( 1.00000*20.0000)*(1.00000+1.00000/(1.00000+std::pow(this->prm_[CaMK__KmCaMK]/CaMK__CaMKa, 8.00000)))) / (1.00000+std::pow(0.500000/this->var_[ionic_concentrations__cajsr], 8.00000));
    double CICR__tau_Jrel1 = CICR__trel1factor;
    if (0.00100000 > CICR__trel1factor) {
        CICR__tau_Jrel1 = 0.00100000;
    }

    double CICR__Jgap = (this->var_[ionic_concentrations__cai] - this->var_[ionic_concentrations__cai2])/this->prm_[CICR__tau_gap];
    double CICR__Rel2 = (CICR__Jgap+ this->var_[CICR__Jrel2]*(this->prm_[cell__vcsr]/this->prm_[cell__vmyo2])) -  SR_uptake__Jup2*(this->prm_[cell__vnsr2]/this->prm_[cell__vmyo2]);
    double CICR__Jrel2_inf = 0.00000;
    if (CICR__Rel2 > 0.00000) {
        CICR__Jrel2_inf = ( 1.00000*( ( 250.000*CICR__Rel2)*(1.00000+1.00000/(1.00000+std::pow(this->prm_[CaMK__KmCaMK]/CaMK__CaMKa, 8.00000))))) / (1.00000+std::pow(0.750000/this->var_[ionic_concentrations__cacsr], 8.00000));
    }
 
    double CICR__trel2factor = ( 50.0000*(1.00000+1.00000/(1.00000+std::pow(this->prm_[CaMK__KmCaMK]/CaMK__CaMKa, 8.00000)))) / (1.00000+std::pow(0.500000/this->var_[ionic_concentrations__cacsr], 8.00000));
    double CICR__tau_Jrel2 = CICR__trel2factor;
    if (0.00100000 > CICR__trel2factor) {
        CICR__tau_Jrel2 = 0.00100000;
    }
    
    double ionic_concentrations__Bcai = 1.00000/((1.00000+( this->prm_[SR_uptake__cmdnmax]*this->prm_[SR_uptake__kmcmdn]) / ( (this->prm_[SR_uptake__kmcmdn]+this->var_[ionic_concentrations__cai])*(this->prm_[SR_uptake__kmcmdn]+this->var_[ionic_concentrations__cai])))+( this->prm_[SR_uptake__trpnmax]*this->prm_[SR_uptake__kmtrpn])/( (this->prm_[SR_uptake__kmtrpn]+this->var_[ionic_concentrations__cai])*(this->prm_[SR_uptake__kmtrpn]+this->var_[ionic_concentrations__cai])));
    double ionic_concentrations__diff_cai =  ionic_concentrations__Bcai*((((  - ((this->cur_[IpCa]+this->cur_[ICab]) -  2.00000*this->cur_[INaCa_i])*this->prm_[cell__Acap])/( ( ( 2.00000*2.00000)*this->prm_[cell__F])*this->prm_[cell__vmyo1]) - ( SR_uptake__Jup*this->prm_[cell__vnsr1])/this->prm_[cell__vmyo1])+( diffusion__Jdiff*this->prm_[cell__vss])/this->prm_[cell__vmyo1]) - CICR__Jgap);
    
    double ionic_concentrations__Bcai2 = 1.00000/((1.00000+( this->prm_[SR_uptake__cmdnmax]*this->prm_[SR_uptake__kmcmdn])/( (this->prm_[SR_uptake__kmcmdn]+this->var_[ionic_concentrations__cai2])*(this->prm_[SR_uptake__kmcmdn]+this->var_[ionic_concentrations__cai2])))+( this->prm_[SR_uptake__trpnmax]*this->prm_[SR_uptake__kmtrpn])/( (this->prm_[SR_uptake__kmtrpn]+this->var_[ionic_concentrations__cai2])*(this->prm_[SR_uptake__kmtrpn]+this->var_[ionic_concentrations__cai2])));
    double ionic_concentrations__diff_cai2 =  ionic_concentrations__Bcai2*(( this->var_[CICR__Jrel2]*(this->prm_[cell__vcsr]/this->prm_[cell__vmyo2])+( CICR__Jgap*this->prm_[cell__vmyo2])/this->prm_[cell__vmyo1]) - ( SR_uptake__Jup2*this->prm_[cell__vnsr2])/this->prm_[cell__vmyo2]);
    
    double ionic_concentrations__Bcass = 1.00000/((1.00000+( this->prm_[SR_uptake__BSRmax]*this->prm_[SR_uptake__KmBSR])/( (this->prm_[SR_uptake__KmBSR]+this->var_[ionic_concentrations__cass])*(this->prm_[SR_uptake__KmBSR]+this->var_[ionic_concentrations__cass])))+( this->prm_[SR_uptake__BSLmax]*this->prm_[SR_uptake__KmBSL])/( (this->prm_[SR_uptake__KmBSL]+this->var_[ionic_concentrations__cass])*(this->prm_[SR_uptake__KmBSL]+this->var_[ionic_concentrations__cass])));
    double ionic_concentrations__diff_cass =  ionic_concentrations__Bcass*(((  - (this->cur_[ICaL] -  2.00000*this->cur_[INaCa_ss])*this->prm_[cell__Acap])/( ( ( 2.00000*2.00000)*this->prm_[cell__F])*this->prm_[cell__vss])+( CICR__Jrel*this->prm_[cell__vjsr])/this->prm_[cell__vss]) - diffusion__Jdiff);

    // Update with Fordward Euler for non-gating variables
    this->var_[CICR__A] = ALGORITHM::ForwardEuler(this->var_[CICR__A], dt, CICR__diff_A);
    this->var_[CICR__tjsrol] = ALGORITHM::ForwardEuler(this->var_[CICR__tjsrol], dt, CICR__diff_tjsrol);
    this->var_[CaMK__CaMKt] = ALGORITHM::ForwardEuler(this->var_[CaMK__CaMKt], dt, CaMK__diff_CaMKt);
    this->var_[ionic_concentrations__ki] = ALGORITHM::ForwardEuler(this->var_[ionic_concentrations__ki], dt, ionic_concentrations__diff_ki);
    this->var_[ionic_concentrations__kss] = ALGORITHM::ForwardEuler(this->var_[ionic_concentrations__kss], dt, ionic_concentrations__diff_kss); 
    this->var_[ionic_concentrations__nai] = ALGORITHM::ForwardEuler(this->var_[ionic_concentrations__nai], dt, ionic_concentrations__diff_nai);
    this->var_[ionic_concentrations__nass] = ALGORITHM::ForwardEuler(this->var_[ionic_concentrations__nass], dt, ionic_concentrations__diff_nass);
    this->var_[ionic_concentrations__cajsr] = ALGORITHM::ForwardEuler(this->var_[ionic_concentrations__cajsr], dt, ionic_concentrations__diff_cajsr);
    this->var_[ionic_concentrations__cacsr] = ALGORITHM::ForwardEuler(this->var_[ionic_concentrations__cacsr], dt, ionic_concentrations__diff_cacsr);
    this->var_[ionic_concentrations__cansr] = ALGORITHM::ForwardEuler(this->var_[ionic_concentrations__cansr], dt, ionic_concentrations__diff_cansr);
    this->var_[ionic_concentrations__cai] = ALGORITHM::ForwardEuler(this->var_[ionic_concentrations__cai], dt, ionic_concentrations__diff_cai);
    this->var_[ionic_concentrations__cai2] = ALGORITHM::ForwardEuler(this->var_[ionic_concentrations__cai2], dt, ionic_concentrations__diff_cai2);
    this->var_[ionic_concentrations__cass] = ALGORITHM::ForwardEuler(this->var_[ionic_concentrations__cass], dt, ionic_concentrations__diff_cass);
    
    // Rush Larsen for gating variables
    this->var_[INaL__hl] = ALGORITHM::RushLarsen(INaL__hl_inf, this->var_[INaL__hl], dt, this->prm_[INaL__tau_hl]);
    this->var_[INaL__ml] = ALGORITHM::RushLarsen(INaL__ml_inf, this->var_[INaL__ml], dt, INaL__tau_ml);
    this->var_[I_Na__m] = ALGORITHM::RushLarsen(I_Na__m_inf, this->var_[I_Na__m], dt, I_Na__tau_m);
    this->var_[I_Na__h] = ALGORITHM::RushLarsen(I_Na__h_inf, this->var_[I_Na__h], dt, I_Na__tau_h);
    this->var_[I_Na__j] = ALGORITHM::RushLarsen(I_Na__j_inf, this->var_[I_Na__j], dt, I_Na__tau_j);
    this->var_[ICaL__d] = ALGORITHM::RushLarsen(ICaL__d_inf, this->var_[ICaL__d], dt, ICaL__tau_d);
    this->var_[ICaL__fca] = ALGORITHM::RushLarsen(ICaL__fca_inf, this->var_[ICaL__fca], dt, ICaL__tau_fca);
    this->var_[IKr__xr] = ALGORITHM::RushLarsen(IKr__xr_inf, this->var_[IKr__xr], dt, IKr__tau_xr);
    this->var_[ITo__aa] = ALGORITHM::RushLarsen(ITo__aa_inf, this->var_[ITo__aa], dt, ITo__tau_aa);
    this->var_[ICaL__ff] = ALGORITHM::RushLarsen(ICaL__ff_inf, this->var_[ICaL__ff], dt, ICaL__tau_ff);
    this->var_[ICaL__fs] = ALGORITHM::RushLarsen(ICaL__fs_inf, this->var_[ICaL__fs], dt, ICaL__tau_fs);
    this->var_[IKs__xs1] = ALGORITHM::RushLarsen(IKs__xs1_inf, this->var_[IKs__xs1], dt, IKs__tau_xs1);
    this->var_[IKs__xs2] = ALGORITHM::RushLarsen(IKs__xs2_inf, this->var_[IKs__xs2], dt, IKs__tau_xs2);
    this->var_[CICR__Jrel1] = ALGORITHM::RushLarsen(CICR__Jrel1_inf, this->var_[CICR__Jrel1], dt, CICR__tau_Jrel1);
    this->var_[CICR__Jrel2] = ALGORITHM::RushLarsen(CICR__Jrel2_inf, this->var_[CICR__Jrel2], dt, CICR__tau_Jrel2);
    
}


std::string Gaur2021::PrintVariables() const 
{
    using namespace Gaur21Var;

    // Create output string stream to pass the variables and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "v = " << this->var_[v] << "\n";
    oss << "dvdt = " << this->var_[dvdt] << "\n";
    oss << "CICR__A = " << this->var_[CICR__A] << "\n";
    oss << "CICR__Jrel1 = " << this->var_[CICR__Jrel1] << "\n";
    oss << "CICR__Jrel2 = " << this->var_[CICR__Jrel2] << "\n";
    oss << "ionic_concentrations__cacsr = " << this->var_[ionic_concentrations__cacsr] << "\n";
    oss << "ionic_concentrations__cai = " << this->var_[ionic_concentrations__cai] << "\n";
    oss << "ionic_concentrations__cai2 = " << this->var_[ionic_concentrations__cai2] << "\n";
    oss << "ionic_concentrations__cajsr = " << this->var_[ionic_concentrations__cajsr] << "\n";
    oss << "ionic_concentrations__cass = " << this->var_[ionic_concentrations__cass] << "\n";
    oss << "CICR__tjsrol = " << this->var_[CICR__tjsrol] << "\n";
    oss << "CaMK__CaMKt = " << this->var_[CaMK__CaMKt] << "\n";
    oss << "ionic_concentrations__ki = " << this->var_[ionic_concentrations__ki] << "\n";
    oss << "ionic_concentrations__nai = " << this->var_[ionic_concentrations__nai] << "\n";
    oss << "ICaL__d = " << this->var_[ICaL__d] << "\n";
    oss << "ICaL__fca = " << this->var_[ICaL__fca] << "\n";
    oss << "ICaL__ff = " << this->var_[ICaL__ff] << "\n";
    oss << "ICaL__fs = " << this->var_[ICaL__fs] << "\n";
    oss << "ionic_concentrations__kss = " << this->var_[ionic_concentrations__kss] << "\n";
    oss << "ionic_concentrations__nass = " << this->var_[ionic_concentrations__nass] << "\n";
    oss << "IKr__xr = " << this->var_[IKr__xr] << "\n";
    oss << "IKs__xs1 = " << this->var_[IKs__xs1] << "\n";
    oss << "IKs__xs2 = " << this->var_[IKs__xs2] << "\n";
    oss << "INaL__hl = " << this->var_[INaL__hl] << "\n";
    oss << "INaL__ml = " << this->var_[INaL__ml] << "\n";
    oss << "ITo__aa = " << this->var_[ITo__aa] << "\n";
    oss << "I_Na__h = " << this->var_[I_Na__h] << "\n";
    oss << "I_Na__j = " << this->var_[I_Na__j] << "\n";
    oss << "I_Na__m = " << this->var_[I_Na__m] << "\n";
    oss << "ionic_concentrations__cansr = " << this->var_[ionic_concentrations__cansr];
    return oss.str();

}


std::string Gaur2021::PrintParameters() const 
{
    using namespace Gaur21Prm;

    // Create output string stream to pass the parameters and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "cell__F = " << this->prm_[cell__F] << "\n";
    oss << "CaMK__KmCaMK = " << this->prm_[CaMK__KmCaMK] << "\n";
    oss << "CICR__SOICR = " << this->prm_[CICR__SOICR] << "\n";
    oss << "CICR__grelbarjsrol = " << this->prm_[CICR__grelbarjsrol] << "\n";
    oss << "CICR__tau_gap = " << this->prm_[CICR__tau_gap] << "\n";
    oss << "CICR__tauoff = " << this->prm_[CICR__tauoff] << "\n";
    oss << "CICR__tauon = " << this->prm_[CICR__tauon] << "\n";
    oss << "CaMK__CaMKo = " << this->prm_[CaMK__CaMKo] << "\n";
    oss << "CaMK__KmCaM = " << this->prm_[CaMK__KmCaM] << "\n";
    oss << "CaMK__PKNa = " << this->prm_[CaMK__PKNa] << "\n";
    oss << "cell__R = " << this->prm_[cell__R] << "\n";
    oss << "cell__T = " << this->prm_[cell__T] << "\n";
    oss << "CaMK__aCaMK = " << this->prm_[CaMK__aCaMK] << "\n";
    oss << "CaMK__bCaMK = " << this->prm_[CaMK__bCaMK] << "\n";
    oss << "cell__cli = " << this->prm_[cell__cli] << "\n";
    oss << "cell__clo = " << this->prm_[cell__clo] << "\n";
    oss << "cell__ko = " << this->prm_[cell__ko] << "\n";
    oss << "cell__nao = " << this->prm_[cell__nao] << "\n";
    oss << "ICaL__PCa = " << this->prm_[ICaL__PCa] << "\n";
    oss << "ICaL__vhalf_d = " << this->prm_[ICaL__vhalf_d] << "\n";
    oss << "ICaL__zca = " << this->prm_[ICaL__zca] << "\n";
    oss << "IKb__GKb = " << this->prm_[IKb__GKb] << "\n";
    oss << "IKs__GKs = " << this->prm_[IKs__GKs] << "\n";
    oss << "INaCa_i__KmCaAct = " << this->prm_[INaCa_i__KmCaAct] << "\n";
    oss << "INaCa_i__kasymm = " << this->prm_[INaCa_i__kasymm] << "\n";
    oss << "INaCa_i__kcaoff = " << this->prm_[INaCa_i__kcaoff] << "\n";
    oss << "INaCa_i__kcaon = " << this->prm_[INaCa_i__kcaon] << "\n";
    oss << "INaCa_i__kna1 = " << this->prm_[INaCa_i__kna1] << "\n";
    oss << "INaCa_i__kna2 = " << this->prm_[INaCa_i__kna2] << "\n";
    oss << "INaCa_i__kna3 = " << this->prm_[INaCa_i__kna3] << "\n";
    oss << "INaCa_i__qca = " << this->prm_[INaCa_i__qca] << "\n";
    oss << "INaCa_i__qna = " << this->prm_[INaCa_i__qna] << "\n";
    oss << "INaCa_i__wca = " << this->prm_[INaCa_i__wca] << "\n";
    oss << "INaCa_i__wna = " << this->prm_[INaCa_i__wna] << "\n";
    oss << "INaCa_i__wnaca = " << this->prm_[INaCa_i__wnaca] << "\n";
    oss << "INaCa_i__zna = " << this->prm_[INaCa_i__zna] << "\n";
    oss << "INaK__H = " << this->prm_[INaK__H] << "\n";
    oss << "INaK__Khp = " << this->prm_[INaK__Khp] << "\n";
    oss << "INaK__Kki = " << this->prm_[INaK__Kki] << "\n";
    oss << "INaK__Kko = " << this->prm_[INaK__Kko] << "\n";
    oss << "INaK__Kmgatp = " << this->prm_[INaK__Kmgatp] << "\n";
    oss << "INaK__Knai0 = " << this->prm_[INaK__Knai0] << "\n";
    oss << "INaK__Knao0 = " << this->prm_[INaK__Knao0] << "\n";
    oss << "INaK__Knap = " << this->prm_[INaK__Knap] << "\n";
    oss << "INaK__Kxkur = " << this->prm_[INaK__Kxkur] << "\n";
    oss << "INaK__MgADP = " << this->prm_[INaK__MgADP] << "\n";
    oss << "INaK__MgATP = " << this->prm_[INaK__MgATP] << "\n";
    oss << "INaK__eP = " << this->prm_[INaK__eP] << "\n";
    oss << "INaK__k1m = " << this->prm_[INaK__k1m] << "\n";
    oss << "INaK__k1p = " << this->prm_[INaK__k1p] << "\n";
    oss << "INaK__k2m = " << this->prm_[INaK__k2m] << "\n";
    oss << "INaK__k2p = " << this->prm_[INaK__k2p] << "\n";
    oss << "INaK__k3m = " << this->prm_[INaK__k3m] << "\n";
    oss << "INaK__k3p2 = " << this->prm_[INaK__k3p2] << "\n";
    oss << "INaK__k4m = " << this->prm_[INaK__k4m] << "\n";
    oss << "INaK__k4p2 = " << this->prm_[INaK__k4p2] << "\n";
    oss << "INaK__zk = " << this->prm_[INaK__zk] << "\n";
    oss << "INaL__GNaL = " << this->prm_[INaL__GNaL] << "\n";
    oss << "INaL__tau_hl = " << this->prm_[INaL__tau_hl] << "\n";
    oss << "INab__PNab = " << this->prm_[INab__PNab] << "\n";
    oss << "I_Na__GNa = " << this->prm_[I_Na__GNa] << "\n";
    oss << "IpCa__GpCa = " << this->prm_[IpCa__GpCa] << "\n";
    oss << "SR_uptake__BSLmax = " << this->prm_[SR_uptake__BSLmax] << "\n";
    oss << "SR_uptake__BSRmax = " << this->prm_[SR_uptake__BSRmax] << "\n";
    oss << "SR_uptake__KmBSL = " << this->prm_[SR_uptake__KmBSL] << "\n";
    oss << "SR_uptake__KmBSR = " << this->prm_[SR_uptake__KmBSR] << "\n";
    oss << "SR_uptake__cmdnmax = " << this->prm_[SR_uptake__cmdnmax] << "\n";
    oss << "SR_uptake__csqnmax = " << this->prm_[SR_uptake__csqnmax] << "\n";
    oss << "SR_uptake__kmcmdn = " << this->prm_[SR_uptake__kmcmdn] << "\n";
    oss << "SR_uptake__kmcsqn = " << this->prm_[SR_uptake__kmcsqn] << "\n";
    oss << "SR_uptake__kmtrpn = " << this->prm_[SR_uptake__kmtrpn] << "\n";
    oss << "SR_uptake__trpnmax = " << this->prm_[SR_uptake__trpnmax] << "\n";
    oss << "cell__L = " << this->prm_[cell__L] << "\n";
    oss << "cell__pi = " << this->prm_[cell__pi] << "\n";
    oss << "cell__rad = " << this->prm_[cell__rad] << "\n";
    oss << "cell__vmyo1frac = " << this->prm_[cell__vmyo1frac] << "\n";
    oss << "stimulus__duration = " << this->prm_[stimulus__duration] << "\n";
    oss << "stimulus__offset = " << this->prm_[stimulus__offset] << "\n";
    oss << "stimulus__period = " << this->prm_[stimulus__period] << "\n";
    oss << "CaMK__ECl = " << this->prm_[CaMK__ECl] << "\n";
    oss << "ICaL__PCaK = " << this->prm_[ICaL__PCaK] << "\n";
    oss << "ICaL__PCaNa = " << this->prm_[ICaL__PCaNa] << "\n";
    oss << "ICaL__vhalff = " << this->prm_[ICaL__vhalff] << "\n";
    oss << "ICab__PCab = " << this->prm_[ICab__PCab] << "\n";
    oss << "IK1__GK1 = " << this->prm_[IK1__GK1] << "\n";
    oss << "IKr__GKr = " << this->prm_[IKr__GKr] << "\n";
    oss << "INaK__delta = " << this->prm_[INaK__delta] << "\n";
    oss << "INaCa_i__Gncx = " << this->prm_[INaCa_i__Gncx] << "\n";
    oss << "ITo__Gto = " << this->prm_[ITo__Gto] << "\n";
    oss << "stimulus__stimulus_amplitude = " << this->prm_[stimulus__stimulus_amplitude] << "\n";
    oss << "INaCa_i__h10 = " << this->prm_[INaCa_i__h10] << "\n";
    oss << "INaCa_i__h11 = " << this->prm_[INaCa_i__h11] << "\n";
    oss << "INaCa_i__h12 = " << this->prm_[INaCa_i__h12] << "\n";
    oss << "INaCa_i__k2 = " << this->prm_[INaCa_i__k2] << "\n";
    oss << "INaCa_i__k5 = " << this->prm_[INaCa_i__k5] << "\n";
    oss << "INaCa_ss__h101 = " << this->prm_[INaCa_ss__h101] << "\n";
    oss << "INaCa_ss__h1111 = " << this->prm_[INaCa_ss__h1111] << "\n";
    oss << "INaCa_ss__h121 = " << this->prm_[INaCa_ss__h121] << "\n";
    oss << "INaCa_ss__k21 = " << this->prm_[INaCa_ss__k21] << "\n";
    oss << "INaCa_ss__k51 = " << this->prm_[INaCa_ss__k51] << "\n";
    oss << "INaK__Pnak = " << this->prm_[INaK__Pnak] << "\n";
    oss << "INaK__a2 = " << this->prm_[INaK__a2] << "\n";
    oss << "INaK__a4 = " << this->prm_[INaK__a4] << "\n";
    oss << "INaK__b1 = " << this->prm_[INaK__b1] << "\n";
    oss << "cell__Ageo = " << this->prm_[cell__Ageo] << "\n";
    oss << "cell__Acap = " << this->prm_[cell__Acap] << "\n";
    oss << "cell__cao = " << this->prm_[cell__cao] << "\n";
    oss << "INaCa_i__k1 = " << this->prm_[INaCa_i__k1] << "\n";
    oss << "INaCa_ss__k11 = " << this->prm_[INaCa_ss__k11] << "\n";
    oss << "cell__vcell = " << this->prm_[cell__vcell] << "\n";
    oss << "cell__vjsr = " << this->prm_[cell__vjsr] << "\n";
    oss << "cell__vcsr = " << this->prm_[cell__vcsr] << "\n";
    oss << "cell__vmyo = " << this->prm_[cell__vmyo] << "\n";
    oss << "cell__vmyo1 = " << this->prm_[cell__vmyo1] << "\n";
    oss << "cell__vmyo2 = " << this->prm_[cell__vmyo2] << "\n";
    oss << "cell__vnsr = " << this->prm_[cell__vnsr] << "\n";
    oss << "cell__vnsr1 = " << this->prm_[cell__vnsr1] << "\n";
    oss << "cell__vnsr2 = " << this->prm_[cell__vnsr2] << "\n";
    oss << "cell__vss = " << this->prm_[cell__vss];
    return oss.str();

}


std::string Gaur2021::PrintCurrents() const
{
    using namespace Gaur21Cur;

    // Create output string stream to pass the currents and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "INa = " << this->cur_[INa] << "\n";
    oss << "INaL = " << this->cur_[INaL] << "\n";
    oss << "ICaL = " << this->cur_[ICaL] << "\n";
    oss << "ICaNa = " << this->cur_[ICaNa] << "\n";
    oss << "ICaK = " << this->cur_[ICaK] << "\n";
    oss << "IKr = " << this->cur_[IKr] << "\n";
    oss << "IKs = " << this->cur_[IKs] << "\n";
    oss << "IK1 = " << this->cur_[IK1] << "\n";
    oss << "ITo = " << this->cur_[ITo] << "\n";
    oss << "INaCa_i = " << this->cur_[INaCa_i] << "\n";
    oss << "INaCa_ss = " << this->cur_[INaCa_ss] << "\n";
    oss << "INaK = " << this->cur_[INaK] << "\n";
    oss << "INab = " << this->cur_[INab] << "\n";
    oss << "IKb = " << this->cur_[IKb] << "\n";
    oss << "IpCa = " << this->cur_[IpCa] << "\n";
    oss << "ICab = " << this->cur_[ICab] << "\n";
    oss << "Iion = " << this->cur_[Gaur21Cur::Iion];
    return oss.str();

}


std::string Gaur2021::PrintBlockCoeffs() const
{
    using namespace Gaur21Cur;

    // Create output string stream to pass the currents and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "INa = " << this->block_coeff_[INa] << "\n";
    oss << "INaL = " << this->block_coeff_[INaL] << "\n";
    oss << "ICaL = " << this->block_coeff_[ICaL] << "\n";
    oss << "ICaNa = " << this->block_coeff_[ICaNa] << "\n";
    oss << "ICaK = " << this->block_coeff_[ICaK] << "\n";
    oss << "IKr = " << this->block_coeff_[IKr] << "\n";
    oss << "IKs = " << this->block_coeff_[IKs] << "\n";
    oss << "IK1 = " << this->block_coeff_[IK1] << "\n";
    oss << "ITo = " << this->block_coeff_[ITo] << "\n";
    oss << "INaCa_i = " << this->block_coeff_[INaCa_i] << "\n";
    oss << "INaCa_ss = " << this->block_coeff_[INaCa_ss] << "\n";
    oss << "INaK = " << this->block_coeff_[INaK] << "\n";
    oss << "INab = " << this->block_coeff_[INab] << "\n";
    oss << "IKb = " << this->block_coeff_[IKb] << "\n";
    oss << "IpCa = " << this->block_coeff_[IpCa] << "\n";
    oss << "ICab = " << this->block_coeff_[ICab] << "\n";
    return oss.str();

}

} // End of namespace ELECTRA
