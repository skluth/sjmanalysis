#ifndef GENEVENTDATA_HH
#define GENEVENTDATA_HH

// This file is adapted from headers in include/HepMC/Data/ in HepMC3
// and can be used to create the dictionary to read hepmc root files
// Copyright (C) 2014 The HepMC collaboration (see AUTHORS for details)

#include <vector>
#include <string>

namespace HepMC3 {

  namespace Units {
    enum MomentumUnit { MEV, GEV };
    enum LengthUnit { MM, CM };
  };

  struct FourVector {
    double m_v1;
    double m_v2;
    double m_v3;
    double m_v4;
  };
  
  struct GenRunInfoData {
    std::vector<std::string> weight_names;    
    std::vector<std::string> tool_name;       
    std::vector<std::string> tool_version;    
    std::vector<std::string> tool_description;
    std::vector<std::string> attribute_name;  
    std::vector<std::string> attribute_string;
  };
  
  struct GenParticleData {
    int pid;
    int status;
    bool is_mass_set;
    double mass;
    FourVector momentum;
  };

  struct GenVertexData {
    int status;
    FourVector position;
  };
  
  struct GenEventData {
    int event_number;
    Units::MomentumUnit momentum_unit;
    Units::LengthUnit length_unit;  
    std::vector<GenParticleData> particles;
    std::vector<GenVertexData> vertices;
    std::vector<double> weights;
    FourVector event_pos;
    std::vector<int> links1;
    std::vector<int> links2;
    std::vector<int> attribute_id;
    std::vector<std::string> attribute_name;
    std::vector<std::string> attribute_string;
  };
  
} // namespace HepMC3

#endif
