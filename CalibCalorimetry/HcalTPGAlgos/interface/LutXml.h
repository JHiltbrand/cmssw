#ifndef CaloOnlineTools_HcalOnlineDb_LutXml_h
#define CaloOnlineTools_HcalOnlineDb_LutXml_h
// -*- C++ -*-
//
// Package:     CaloOnlineTools/HcalOnlineDb
// Class  :     LutXml
//
/**\class LutXml LutXml.h CalibCalorimetry/HcalTPGAlgos/interface/LutXml.h

 Description: <one line class summary>

 Usage:
    <usage>

*/
//
// Original Author:  Gena Kukartsev, kukarzev@fnal.gov
//         Created:  Tue Mar 18 14:30:33 CDT 2008
//

#include "CalibCalorimetry/HcalTPGAlgos/interface/XMLDOMBlock.h"
#include "CondFormats/HcalObjects/interface/HcalElectronicsMap.h"
#include "DataFormats/DetId/interface/DetId.h"

#include <cstdint>
#include <map>
#include <vector>

class LutXml : public XMLDOMBlock {
public:
  typedef struct _Config {
    _Config();
    std::string infotype;
    int ieta, iphi, depth, crate, slot, topbottom, fiber, fiberchan, lut_type;
    std::string creationtag;
    std::string creationstamp;
    std::string formatrevision;
    std::string targetfirmware;
    int generalizedindex;
    int weight;
    int codedvetothreshold;
    std::vector<unsigned int> lut;
    std::vector<uint64_t> mask;
  } Config;

  LutXml();
  LutXml(XERCES_CPP_NAMESPACE::InputSource& _source);
  LutXml(std::string filename);
  ~LutXml() override;

  void init(void);
  void addLut(Config& _config, XMLDOMBlock* checksums_xml = nullptr);
  std::string& getCurrentBrick(void);

  DetId detid_from_crate(int crate, int slot, int fiber, int fiberch, bool isTrigger, const HcalElectronicsMap* emap);
  int a_to_i(char* inbuf);
  int create_lut_map(const HcalElectronicsMap* emap);

  static std::string get_checksum(std::vector<unsigned int>& lut);

  typedef std::map<uint32_t, std::vector<unsigned int> >::const_iterator const_iterator;
  const_iterator begin() const;
  const_iterator end() const;
  const_iterator find(uint32_t) const;

protected:
  XMLCh* root;
  XMLCh* brick;
  XERCES_CPP_NAMESPACE::DOMElement* addParameter(std::string _name, std::string _type, std::string _value);
  XERCES_CPP_NAMESPACE::DOMElement* addParameter(std::string _name, std::string _type, int _value);

  template <typename T>
  XERCES_CPP_NAMESPACE::DOMElement* addData(std::string _elements, std::string _encoding, const T& _lut);

  XERCES_CPP_NAMESPACE::DOMElement* add_checksum(XERCES_CPP_NAMESPACE::DOMDocument* parent, Config& config);

  XERCES_CPP_NAMESPACE::DOMElement* brickElem;

  std::map<uint32_t, std::vector<unsigned int> > lut_map;
};

#endif
