#include "ProteinConfig.hpp"

#include <fstream>

#include <yaml-cpp/yaml.h>

ProteinConfig::ProteinConfig(std::string proteinConfigFile) {
    YAML::Node node = YAML::LoadFile(proteinConfigFile);

    const int nTypes = node["proteins"].size();
    types.resize(nTypes);
    freeNumber.resize(nTypes);
    fixedLocations.resize(nTypes);

    for (int i = 0; i < nTypes; i++) {
        const auto &p = node["proteins"][i];
        auto &type = types[i];
        // int
        type.tag = p["tag"].as<int>();
        tagLookUp[type.tag] = i;
        // bool
        type.walkOff = p["walkOff"].as<bool>();
        type.bindAntiParallel = p["bindAntiParallel"].as<bool>();
        type.fixedEnd0 = p["fixedEnd0"].as<bool>();
        // double
        type.freeLength = p["freeLength"].as<double>();
        type.kappa = p["kappa"].as<double>();
        type.fstall = p["fstall"].as<double>();
        type.lambda = p["lambda"].as<double>();
        type.diffUnbound = p["diffUnbound"].as<double>();
        type.rc = p["rc"].as<double>();
        // array[2]
        type.vmax[0] = p["vmax"][0].as<double>();
        type.vmax[1] = p["vmax"][1].as<double>();
        type.eps = p["eps"].as<double>();
        type.Ka[0] = p["Ka"][0].as<double>() / 602;
        type.Ka[1] = p["Ka"][1].as<double>() / 602;
        type.Ke[0] = p["Ke"][0].as<double>();
        type.Ke[1] = p["Ke"][1].as<double>();
        type.ko_s[0] = p["ko_s"][0].as<double>();
        type.ko_s[1] = p["ko_s"][1].as<double>();
        type.ko_d[0] = p["ko_d"][0].as<double>();
        type.ko_d[1] = p["ko_d"][1].as<double>();

        freeNumber[i] = p["freeNumber"].as<int>(); // int
        // vector of double
        const int nFixedPerMT = p["fixedLocationPerMT"].size();
        fixedLocations[i].resize(nFixedPerMT);
        for (int j = 0; j < nFixedPerMT; j++) {
            fixedLocations[i][j] = p["fixedLocationPerMT"][j].as<double>();
        }
    }

    // sanity check
    for (int iType = 0; iType < nTypes; iType++) {
        if (types[iType].fixedEnd0 && freeNumber[iType] != 0) {
            printf("Protein Type %d setting error\n", iType);
            printf("Proteins with fixedEnd0 must set freeNumber = 0\n");
            printf("Setting freeNumber = 0 for this type\n");
            freeNumber[iType] = 0;
        }
    }
}

ProteinConfig::~ProteinConfig() {
    freeNumber.clear();
    fixedLocations.clear();
    for (auto &t : types) {
        if (t.LUTablePtr)
            delete t.LUTablePtr;
    }
    types.clear();
}

void ProteinConfig::echo() const {
    const int nTypes = types.size();
    printf("%d Types of proteins\n", nTypes);

    for (int i = 0; i < nTypes; i++) {
        types[i].echo();
        printf("freeNumber: %d\n", freeNumber[i]);
        printf("fixedLocationPerMT: ");
        for (auto &l : fixedLocations[i]) {
            printf("%g,", l);
        }
        printf("\n");
    }
    printf("----------------------------------\n");
    printf("<<< Warning >>> \n");
    printf("bindAntiParalell not implemented yet \n");
    printf("fixedEnd0        not implemented yet \n");
    printf("----------------------------------\n");
}
