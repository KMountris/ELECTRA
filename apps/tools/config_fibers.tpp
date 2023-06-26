/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

#ifndef ELECTRA_APPS_TOOLS_CONFIG_FIBERS_TPP_
#define ELECTRA_APPS_TOOLS_CONFIG_FIBERS_TPP_

#include "config_fibers.hpp"

namespace APP_ELECTRA
{

template<short DIM, short CELL_NODES>
ConfigFibers<DIM,CELL_NODES>::ConfigFibers()
{}


template<short DIM, short CELL_NODES>
ConfigFibers<DIM,CELL_NODES>::~ConfigFibers()
{}


template<short DIM, short CELL_NODES>
void ConfigFibers<DIM,CELL_NODES>::SetVentriFiberRules(const Parser &parser, VentriFiberRules &ventri_fiber_rules, std::ostream &stream) const
{
    stream << ELECTRA::Logger::Message("Assigned ventricular fiber rules:\n");

    // Obtain alpha angle from venticular rules if available.
    if (parser.HasAttribute("fibers.ventricular rules.alpha angle")) {
        stream << ELECTRA::Logger::Message("Alpha angle rules:\n");
    } else {
        stream << ELECTRA::Logger::Message("The default alpha angle rules are used.\n");
    }

    // Get alpha angle rules.
    if (parser.HasAttribute("fibers.ventricular rules.alpha angle.left endo wall")) {
        double alpha_lv_endo = parser.GetValue<double>("fibers.ventricular rules.alpha angle.left endo wall");
        ventri_fiber_rules.SetAlphaLvEndo(alpha_lv_endo);
        stream << ELECTRA::Logger::Message("        - left endo wall: ") << alpha_lv_endo << " deg\n";
    }
    if (parser.HasAttribute("fibers.ventricular rules.alpha angle.left epi wall")) {
        double alpha_lv_epi = parser.GetValue<double>("fibers.ventricular rules.alpha angle.left epi wall");
        ventri_fiber_rules.SetAlphaLvEpi(alpha_lv_epi);
        stream << ELECTRA::Logger::Message("        - left epi wall: ") << alpha_lv_epi << " deg\n";
    }
    if (parser.HasAttribute("fibers.ventricular rules.alpha angle.right endo wall")) {
        double alpha_rv_endo = parser.GetValue<double>("fibers.ventricular rules.alpha angle.right endo wall");
        ventri_fiber_rules.SetAlphaRvEndo(alpha_rv_endo);
        stream << ELECTRA::Logger::Message("        - right endo wall: ") << alpha_rv_endo << " deg\n";
    }
    if (parser.HasAttribute("fibers.ventricular rules.alpha angle.right epi wall")) {
        double alpha_rv_epi = parser.GetValue<double>("fibers.ventricular rules.alpha angle.right epi wall");
        ventri_fiber_rules.SetAlphaRvEpi(alpha_rv_epi);
        stream << ELECTRA::Logger::Message("        - right epi wall: ") << alpha_rv_epi << " deg\n";
    }
    if (parser.HasAttribute("fibers.ventricular rules.alpha angle.left endo outflow tract")) {
        double alpha_lv_ot_endo = parser.GetValue<double>("fibers.ventricular rules.alpha angle.left endo outflow tract");
        ventri_fiber_rules.SetAlphaOtLvEndo(alpha_lv_ot_endo);
        stream << ELECTRA::Logger::Message("        - left endo outflow tract: ") << alpha_lv_ot_endo << " deg\n";
    }
    if (parser.HasAttribute("fibers.ventricular rules.alpha angle.left epi outflow tract")) {
        double alpha_lv_ot_epi = parser.GetValue<double>("fibers.ventricular rules.alpha angle.left epi outflow tract");
        ventri_fiber_rules.SetAlphaOtLvEpi(alpha_lv_ot_epi);
        stream << ELECTRA::Logger::Message("        - left epi outflow tract: ") << alpha_lv_ot_epi << " deg\n";
    }
    if (parser.HasAttribute("fibers.ventricular rules.alpha angle.right endo outflow tract")) {
        double alpha_rv_ot_endo = parser.GetValue<double>("fibers.ventricular rules.alpha angle.right endo outflow tract");
        ventri_fiber_rules.SetAlphaOtRvEndo(alpha_rv_ot_endo);
        stream << ELECTRA::Logger::Message("        - right endo outflow tract: ") << alpha_rv_ot_endo << " deg\n";
    }
    if (parser.HasAttribute("fibers.ventricular rules.alpha angle.right epi outflow tract")) {
        double alpha_rv_ot_epi = parser.GetValue<double>("fibers.ventricular rules.alpha angle.right epi outflow tract");
        ventri_fiber_rules.SetAlphaOtRvEpi(alpha_rv_ot_epi);
        stream << ELECTRA::Logger::Message("        - right epi outflow tract: ") << alpha_rv_ot_epi << " deg\n";
    }
    if (parser.HasAttribute("fibers.ventricular rules.alpha angle.septum wall")) {
        double alpha_septum = parser.GetValue<double>("fibers.ventricular rules.alpha angle.septum wall");
        ventri_fiber_rules.SetAlphaSeptum(alpha_septum);
        stream << ELECTRA::Logger::Message("        - septum wall: ") << alpha_septum << " deg\n";
    }

    // Obtain beta angle from venticular rules if available.
    if (parser.HasAttribute("fibers.ventricular rules.beta angle")) {
        stream << ELECTRA::Logger::Message("Beta angle rules:\n");
    } else {
        stream << ELECTRA::Logger::Message("The default beta angle rules are used.\n");
    }

    // Get beta angle rules.
    if (parser.HasAttribute("fibers.ventricular rules.beta angle.left endo wall")) {
        double beta_lv_endo = parser.GetValue<double>("fibers.ventricular rules.beta angle.left endo wall");
        ventri_fiber_rules.SetBetaLvEndo(beta_lv_endo);
        stream << ELECTRA::Logger::Message("        - left endo wall: ") << beta_lv_endo << " deg\n";
    }
    if (parser.HasAttribute("fibers.ventricular rules.beta angle.left epi wall")) {
        double beta_lv_epi = parser.GetValue<double>("fibers.ventricular rules.beta angle.left epi wall");
        ventri_fiber_rules.SetBetaLvEpi(beta_lv_epi);
        stream << ELECTRA::Logger::Message("        - left epi wall: ") << beta_lv_epi << " deg\n";
    }
    if (parser.HasAttribute("fibers.ventricular rules.beta angle.right endo wall")) {
        double beta_rv_endo = parser.GetValue<double>("fibers.ventricular rules.beta angle.right endo wall");
        ventri_fiber_rules.SetBetaRvEndo(beta_rv_endo);
        stream << ELECTRA::Logger::Message("        - right endo wall: ") << beta_rv_endo << " deg\n";
    }
    if (parser.HasAttribute("fibers.ventricular rules.beta angle.right epi wall")) {
        double beta_rv_epi = parser.GetValue<double>("fibers.ventricular rules.beta angle.right epi wall");
        ventri_fiber_rules.SetBetaRvEpi(beta_rv_epi);
        stream << ELECTRA::Logger::Message("        - right epi wall: ") << beta_rv_epi << " deg\n";
    }
    if (parser.HasAttribute("fibers.ventricular rules.beta angle.left endo outflow tract")) {
        double beta_lv_ot_endo = parser.GetValue<double>("fibers.ventricular rules.beta angle.left endo outflow tract");
        ventri_fiber_rules.SetBetaOtLvEndo(beta_lv_ot_endo);
        stream << ELECTRA::Logger::Message("        - left endo outflow tract: ") << beta_lv_ot_endo << " deg\n";
    }
    if (parser.HasAttribute("fibers.ventricular rules.beta angle.left epi outflow tract")) {
        double beta_lv_ot_epi = parser.GetValue<double>("fibers.ventricular rules.beta angle.left epi outflow tract");
        ventri_fiber_rules.SetBetaOtLvEpi(beta_lv_ot_epi);
        stream << ELECTRA::Logger::Message("        - left epi outflow tract: ") << beta_lv_ot_epi << " deg\n";
    }
    if (parser.HasAttribute("fibers.ventricular rules.beta angle.right endo outflow tract")) {
        double beta_rv_ot_endo = parser.GetValue<double>("fibers.ventricular rules.beta angle.right endo outflow tract");
        ventri_fiber_rules.SetBetaOtRvEndo(beta_rv_ot_endo);
        stream << ELECTRA::Logger::Message("        - right endo outflow tract: ") << beta_rv_ot_endo << " deg\n";
    }
    if (parser.HasAttribute("fibers.ventricular rules.beta angle.right epi outflow tract")) {
        double beta_rv_ot_epi = parser.GetValue<double>("fibers.ventricular rules.beta angle.right epi outflow tract");
        ventri_fiber_rules.SetBetaOtRvEpi(beta_rv_ot_epi);
        stream << ELECTRA::Logger::Message("        - right epi outflow tract: ") << beta_rv_ot_epi << " deg\n";
    }

}


template<short DIM, short CELL_NODES>
void ConfigFibers<DIM,CELL_NODES>::SetVentriTags(const Parser &parser, VentriTags &ventri_tags, std::ostream &stream) const
{
    // Parse ventricular tags.
    std::string lv_endo = parser.GetValue<std::string>("fibers.ventricular tags.left endo wall");
    std::string rv_endo = parser.GetValue<std::string>("fibers.ventricular tags.right endo wall");
    std::string lv_apex = parser.GetValue<std::string>("fibers.ventricular tags.left apex");
    std::string rv_apex = parser.GetValue<std::string>("fibers.ventricular tags.right apex");
    std::string epi = parser.GetValue<std::string>("fibers.ventricular tags.epi wall");
    std::string mitral_ring = parser.GetValue<std::string>("fibers.ventricular tags.mitral valve");
    std::string tricuspid_ring = parser.GetValue<std::string>("fibers.ventricular tags.tricuspid valve");

    // Assign ventricular tags.
    ventri_tags.SetLvEndoWallTag(lv_endo);
    ventri_tags.SetRvEndoWallTag(rv_endo);
    ventri_tags.SetLvApexTag(lv_apex);
    ventri_tags.SetRvApexTag(rv_apex);
    ventri_tags.SetEpiWallTag(epi);
    ventri_tags.SetMitralValveTag(mitral_ring);
    ventri_tags.SetTricuspidValveTag(tricuspid_ring);

    // Set optional tags for biventricular geometries.
    std::string aortic_ring = "";
    if (parser.HasAttribute("fibers.ventricular tags.aortic valve")) {
        aortic_ring = parser.GetValue<std::string>("fibers.ventricular tags.aortic valve");
        ventri_tags.SetAorticValveTag(aortic_ring);
    }
    std::string pulmonary_ring = "";
    if (parser.HasAttribute("fibers.ventricular tags.pulmonary valve")) {
        pulmonary_ring = parser.GetValue<std::string>("fibers.ventricular tags.pulmonary valve");
        ventri_tags.SetPulmonaryValveTag(pulmonary_ring);
    }

    // Print ventricular tags information.
    stream << ELECTRA::Logger::Message("Assigned ventricular tags:\n");
    stream << ELECTRA::Logger::Message("        - ") << lv_endo << "\n";
    stream << ELECTRA::Logger::Message("        - ") << rv_endo << "\n";
    stream << ELECTRA::Logger::Message("        - ") << lv_apex << "\n";
    stream << ELECTRA::Logger::Message("        - ") << rv_apex << "\n";
    stream << ELECTRA::Logger::Message("        - ") << epi << "\n";
    stream << ELECTRA::Logger::Message("        - ") << mitral_ring << "\n";
    if (parser.HasAttribute("fibers.ventricular tags.aortic valve"))
        stream << ELECTRA::Logger::Message("        - ") << aortic_ring << "\n";
    stream << ELECTRA::Logger::Message("        - ") << tricuspid_ring << "\n";
    if (parser.HasAttribute("fibers.ventricular tags.pulmonary valve"))
        stream << ELECTRA::Logger::Message("        - ") << pulmonary_ring << "\n";

}


template<short DIM, short CELL_NODES>
void ConfigFibers<DIM,CELL_NODES>::SetAtriFiberRules(const Parser &parser, AtriFiberRules &atri_fiber_rules, std::ostream &stream) const
{
    stream << ELECTRA::Logger::Message("Assigned atrial fiber rules:\n");

    // Get angle rules.
    if (parser.HasAttribute("fibers.atrial rules.tricuspid valve")) {
        double angle = parser.GetValue<double>("fibers.atrial rules.tricuspid valve");
        atri_fiber_rules.SetTricuspidValveRule(angle);
        stream << ELECTRA::Logger::Message("        - tricuspid valve: ") << angle << " deg\n";
    }
    if (parser.HasAttribute("fibers.atrial rules.inferior caval vein")) {
        double angle = parser.GetValue<double>("fibers.atrial rules.inferior caval vein");
        atri_fiber_rules.SetInfCavalVeinRule(angle);
        stream << ELECTRA::Logger::Message("        - inferior caval vein: ") << angle << " deg\n";
    }
    if (parser.HasAttribute("fibers.atrial rules.superior caval vein")) {
        double angle = parser.GetValue<double>("fibers.atrial rules.superior caval vein");
        atri_fiber_rules.SetSupCavalVeinRule(angle);
        stream << ELECTRA::Logger::Message("        - superior caval vein: ") << angle << " deg\n";
    }
    if (parser.HasAttribute("fibers.atrial rules.crista terminal plus")) {
        double angle = parser.GetValue<double>("fibers.atrial rules.crista terminal plus");
        atri_fiber_rules.SetCristaTerminalPlusRule(angle);
        stream << ELECTRA::Logger::Message("        - crista terminal plus: ") << angle << " deg\n";
    }
    if (parser.HasAttribute("fibers.atrial rules.crista terminal minus")) {
        double angle = parser.GetValue<double>("fibers.atrial rules.crista terminal minus");
        atri_fiber_rules.SetCristaTerminalMinusRule(angle);
        stream << ELECTRA::Logger::Message("        - crista terminal minus: ") << angle << " deg\n";
    }
    if (parser.HasAttribute("fibers.atrial rules.inter caval bundle")) {
        double angle = parser.GetValue<double>("fibers.atrial rules.inter caval bundle");
        atri_fiber_rules.SetInterCavalBundleRule(angle);
        stream << ELECTRA::Logger::Message("        - inter caval bundle: ") << angle << " deg\n";
    }
    if (parser.HasAttribute("fibers.atrial rules.right septum wall")) {
        double angle = parser.GetValue<double>("fibers.atrial rules.right septum wall");
        atri_fiber_rules.SetRightSeptumWallRule(angle);
        stream << ELECTRA::Logger::Message("        - right septum wall: ") << angle << " deg\n";
    }
    if (parser.HasAttribute("fibers.atrial rules.right lateral wall")) {
        double angle = parser.GetValue<double>("fibers.atrial rules.right lateral wall");
        atri_fiber_rules.SetRightLateralWallRule(angle);
        stream << ELECTRA::Logger::Message("        - right lateral wall: ") << angle << " deg\n";
    }
    if (parser.HasAttribute("fibers.atrial rules.mitral valve")) {
        double angle = parser.GetValue<double>("fibers.atrial rules.mitral valve");
        atri_fiber_rules.SetMitralValveRule(angle);
        stream << ELECTRA::Logger::Message("        - mitral valve: ") << angle << " deg\n";
    }
    if (parser.HasAttribute("fibers.atrial rules.left pulmonary vein")) {
        double angle = parser.GetValue<double>("fibers.atrial rules.left pulmonary vein");
        atri_fiber_rules.SetLeftPulmonVeinRule(angle);
        stream << ELECTRA::Logger::Message("        - left pulmonary vein: ") << angle << " deg\n";
    }
    if (parser.HasAttribute("fibers.atrial rules.right pulmonary vein")) {
        double angle = parser.GetValue<double>("fibers.atrial rules.right pulmonary vein");
        atri_fiber_rules.SetRightPulmonVeinRule(angle);
        stream << ELECTRA::Logger::Message("        - right pulmonary vein: ") << angle << " deg\n";
    }

}


template<short DIM, short CELL_NODES>
void ConfigFibers<DIM,CELL_NODES>::SetAtriTags(const Parser &parser, AtriTags &atri_tags, std::ostream &stream) const
{
    // Parse atrial tags.
    std::string la_endo = parser.GetValue<std::string>("fibers.atrial tags.left endo wall");
    std::string ra_endo = parser.GetValue<std::string>("fibers.atrial tags.right endo wall");
    std::string la_appendage = parser.GetValue<std::string>("fibers.atrial tags.left appendage");
    std::string ra_appendage = parser.GetValue<std::string>("fibers.atrial tags.right appendage");
    std::string ra_top_endo_ring = parser.GetValue<std::string>("fibers.atrial tags.right top endo");
    std::string ra_top_epi_ring = parser.GetValue<std::string>("fibers.atrial tags.right top epi");
    std::string epi = parser.GetValue<std::string>("fibers.atrial tags.epi wall");
    std::string mitral_ring = parser.GetValue<std::string>("fibers.atrial tags.mitral valve");
    std::string left_pulmon_vein_ring = parser.GetValue<std::string>("fibers.atrial tags.left pulmonary vein");
    std::string right_pulmon_vein_ring = parser.GetValue<std::string>("fibers.atrial tags.right pulmonary vein");
    std::string sup_caval_vein_ring = parser.GetValue<std::string>("fibers.atrial tags.superior caval vein");
    std::string inf_caval_vein_ring = parser.GetValue<std::string>("fibers.atrial tags.inferior caval vein");
    std::string tricuspid_septum_ring = parser.GetValue<std::string>("fibers.atrial tags.tricuspid valve septum");
    std::string tricuspid_free_ring = parser.GetValue<std::string>("fibers.atrial tags.tricuspid valve free");

    // Assign atrial tags.
    atri_tags.SetLaEndoWallTag("la_endo");
    atri_tags.SetRaEndoWallTag("ra_endo");
    atri_tags.SetLaAppendageTag("la_appendage");
    atri_tags.SetRaAppendageTag("ra_appendage");
    atri_tags.SetRaTopEndoTag("ra_top_endo_ring");
    atri_tags.SetRaTopEpiTag("ra_top_epi_ring");
    atri_tags.SetEpiWallTag("epi");
    atri_tags.SetMitralValveTag("mitral_ring");
    atri_tags.SetLeftPulmonVeinTag("left_pulmon_vein_ring");
    atri_tags.SetRightPulmonVeinTag("right_pulmon_vein_ring");
    atri_tags.SetSuperiorCavalVeinTag("sup_caval_vein_ring");
    atri_tags.SetInferiorCavalVeinTag("inf_caval_vein_ring");
    atri_tags.SetTricuspidValveSeptumTag("tricuspid_septum_ring");
    atri_tags.SetTricuspidValveFreeTag("tricuspid_free_ring");

    // Print atrial tags information.
    stream << ELECTRA::Logger::Message("Assigned atrial tags:\n");
    stream << ELECTRA::Logger::Message("        - ") << la_endo << "\n";
    stream << ELECTRA::Logger::Message("        - ") << ra_endo << "\n";
    stream << ELECTRA::Logger::Message("        - ") << la_appendage << "\n";
    stream << ELECTRA::Logger::Message("        - ") << ra_appendage << "\n";
    stream << ELECTRA::Logger::Message("        - ") << ra_top_endo_ring << "\n";
    stream << ELECTRA::Logger::Message("        - ") << ra_top_epi_ring << "\n";
    stream << ELECTRA::Logger::Message("        - ") << epi << "\n";
    stream << ELECTRA::Logger::Message("        - ") << mitral_ring << "\n";
    stream << ELECTRA::Logger::Message("        - ") << left_pulmon_vein_ring << "\n";
    stream << ELECTRA::Logger::Message("        - ") << right_pulmon_vein_ring << "\n";
    stream << ELECTRA::Logger::Message("        - ") << sup_caval_vein_ring << "\n";
    stream << ELECTRA::Logger::Message("        - ") << inf_caval_vein_ring << "\n";
    stream << ELECTRA::Logger::Message("        - ") << tricuspid_septum_ring << "\n";
    stream << ELECTRA::Logger::Message("        - ") << tricuspid_free_ring << "\n";
}


template<short DIM, short CELL_NODES>
void ConfigFibers<DIM,CELL_NODES>::ComputeVentriFibers(const Parser &parser, const IMP::Mesh<DIM,CELL_NODES> &mesh, const IMP::Voronoi<DIM> &voro,
                                                       const CLOUDEA::Fpm<DIM> &fpm, const VentriTags &tags, const VentriFiberRules &rules,
                                                       Ldrbm<DIM,CELL_NODES> &ldrbm, std::ostream &stream) const
{
    // Get numerical approximation method.
    std::string method = parser.GetValue<std::string>("numerical approximation.method");
    std::transform(std::begin(method), std::end(method), std::begin(method), ::tolower);

    // Assemble matrices and compute transmurality.
    Timer timer;
    if (method == "fpm") {
        timer.Reset();
        stream << Logger::Message("Assembling Laplacian matrices... ");
        ldrbm.AssembleLaplacian(voro, fpm, fpm.Penalty());
        stream << timer.PrintElapsedTime() << "\n";

        timer.Reset();
        stream << Logger::Message("Computing transmural distance and gradient... ");
        ldrbm.ComputeVentriTransmural(voro, fpm, tags);
        stream << timer.PrintElapsedTime() << "\n";
    } else {
        timer.Reset();
        stream << Logger::Message("Assembling Laplacian matrices... ");
        ldrbm.AssembleLaplacian(mesh);
        stream << timer.PrintElapsedTime() << "\n";

        timer.Reset();
        stream << Logger::Message("Computing transmural distance and gradient... ");
        ldrbm.ComputeVentriTransmural(mesh, tags);
        stream << timer.PrintElapsedTime() << "\n";
    }

    // Partition the volume of geometry.
    stream << Logger::Message("Partioning geometry...\n");
    ldrbm.ExtractVentricleNodalPartitions(mesh);

    // Get septum threshold value.
    double thres = parser.GetValue<double>("fibers.septum threshold");

    // Compute additional distance fields.
    if (method == "fpm") {
        timer.Reset();
        stream << Logger::Message("Computing septal distance and gradient... ");
        ldrbm.ComputeVentriSeptal(voro, tags, thres);
        stream << timer.PrintElapsedTime() << "\n";

        timer.Reset();
        stream << Logger::Message("Computing intraventricular function... ");
        ldrbm.ComputeVentriIntraventricular(voro, tags);
        stream << timer.PrintElapsedTime() << "\n";

        timer.Reset();
        stream << Logger::Message("Computing apicobasal distance and gradient... ");
        ldrbm.ComputeVentriApicoBasal(voro, fpm, tags);
        stream << timer.PrintElapsedTime() << "\n";

    } else {
        timer.Reset();
        stream << Logger::Message("Computing septal distance and gradient... ");
        ldrbm.ComputeVentriSeptal(mesh, tags, thres);
        stream << timer.PrintElapsedTime() << "\n";

        timer.Reset();
        stream << Logger::Message("Computing intraventricular function... ");
        ldrbm.ComputeVentriIntraventricular(mesh, tags);
        stream << timer.PrintElapsedTime() << "\n";

        timer.Reset();
        stream << Logger::Message("Computing apicobasal distance and gradient... ");
        ldrbm.ComputeVentriApicoBasal(mesh, tags);
        stream << timer.PrintElapsedTime() << "\n";
    }

    // Compute the fibers direction.
    timer.Reset();
    stream << Logger::Message("Computing ventricular fibers... ");
    ldrbm.ComputeVentriFibers(rules);
    stream << timer.PrintElapsedTime() << "\n";

}


template<short DIM, short CELL_NODES>
void ConfigFibers<DIM,CELL_NODES>::ComputeAtriFibers(const Parser &parser, const IMP::Mesh<DIM,CELL_NODES> &mesh, const IMP::Voronoi<DIM> &voro,
                                                     const CLOUDEA::Fpm<DIM> &fpm, const AtriTags &tags, const AtriFiberRules &rules,
                                                     Ldrbm<DIM,CELL_NODES> &ldrbm, std::ostream &stream) const
{
    // Get numerical approximation method.
    std::string method = parser.GetValue<std::string>("numerical approximation.method");
    std::transform(std::begin(method), std::end(method), std::begin(method), ::tolower);

    Timer timer;
    // Assemble matrices and compute transmurality.
    if (method == "fpm") {
        timer.Reset();
        stream << Logger::Message("Assembling Laplacian matrices... ");
        ldrbm.AssembleLaplacian(voro, fpm, fpm.Penalty());
        stream << timer.PrintElapsedTime() << "\n";

        timer.Reset();
        stream << Logger::Message("Computing transmural distance and gradient... ");
        ldrbm.ComputeAtriTransmural(voro, fpm, tags);
        stream << timer.PrintElapsedTime() << "\n";

    } else {
        timer.Reset();
        stream << Logger::Message("Assembling Laplacian matrices... ");
        ldrbm.AssembleLaplacian(mesh);
        stream << timer.PrintElapsedTime() << "\n";

        timer.Reset();
        stream << Logger::Message("Computing transmural distance and gradient... ");
        ldrbm.ComputeAtriTransmural(mesh, tags);
        stream << timer.PrintElapsedTime() << "\n";
    }

    // Partition the volume of geometry.
    stream << Logger::Message("Partioning geometry...\n");
    ldrbm.ExtractAtriumNodalPartitions(mesh);
    stream << timer.PrintElapsedTime() << "\n";

    // Compute additional distance fields.
    if (method == "fpm") {
        timer.Reset();
        stream << Logger::Message("Computing appendage to veins distance and gradient... ");
        ldrbm.ComputeAtriAppendageVeins(voro, fpm, tags);
        stream << timer.PrintElapsedTime() << "\n";

        timer.Reset();
        stream << Logger::Message("Computing interveins distance and gradient... ");
        ldrbm.ComputeAtriInterVeins(voro, fpm, tags);
        stream << timer.PrintElapsedTime() << "\n";

        timer.Reset();
        stream << Logger::Message("Computing valve to veins distance and gradient... ");
        ldrbm.ComputeAtriValveVeins(voro, fpm, tags);
        stream << timer.PrintElapsedTime() << "\n";

        timer.Reset();
        stream << Logger::Message("Computing atria tricuspid distance and gradient... ");
        ldrbm.ComputeAtriTricuspid(voro, fpm, tags);
        stream << timer.PrintElapsedTime() << "\n";

    } else {
        timer.Reset();
        stream << Logger::Message("Computing appendage to veins distance and gradient... ");
        ldrbm.ComputeAtriAppendageVeins(mesh, tags);
        stream << timer.PrintElapsedTime() << "\n";

        timer.Reset();
        stream << Logger::Message("Computing interveins distance and gradient... ");
        ldrbm.ComputeAtriInterVeins(mesh, tags);
        stream << timer.PrintElapsedTime() << "\n";

        timer.Reset();
        stream << Logger::Message("Computing valve to veins distance and gradient... ");
        ldrbm.ComputeAtriValveVeins(mesh, tags);
        stream << timer.PrintElapsedTime() << "\n";

        timer.Reset();
        stream << Logger::Message("Computing atria tricuspid distance and gradient... ");
        ldrbm.ComputeAtriTricuspid(mesh, tags);
        stream << timer.PrintElapsedTime() << "\n";
    }

    // Compute the fibers direction.
    timer.Reset();
    stream << Logger::Message("Computing atrial fibers... ");
    ldrbm.ComputeAtriFibers(rules);
    stream << timer.PrintElapsedTime() << "\n";
}


} // end of namespace APP_ELECTRA

#endif //ELECTRA_APPS_TOOLS_CONFIG_FIBERS_TPP_