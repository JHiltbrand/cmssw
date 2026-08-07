#include <cstdlib>
#include <fstream>
#include <memory>
#include <sstream>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>

#include <cstdlib>  // For srand() and rand()

#ifdef HAVE_XDAQ
#include <toolbox/string.h>
#else
#include "CaloOnlineTools/HcalOnlineDb/interface/xdaq_compat.h"  // Replaces toolbox::toString
#endif

#include "OnlineDB/Oracle/interface/Oracle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CaloOnlineTools/HcalOnlineDb/interface/HcalLutManager.h"
#include "CaloOnlineTools/HcalOnlineDb/interface/ZdcLut.h"
#include "CalibCalorimetry/HcalTPGAlgos/interface/XMLProcessor.h"
#include "CalibCalorimetry/HcalTPGAlgos/interface/XMLDOMBlock.h"
#include "CaloOnlineTools/HcalOnlineDb/interface/HcalQIEManager.h"
#include "CaloOnlineTools/HcalOnlineDb/interface/LMap.h"
#include "CaloOnlineTools/HcalOnlineDb/interface/XMLLUTLoader.h"
#include "CaloOnlineTools/HcalOnlineDb/interface/RooGKCounter.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"

#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
XERCES_CPP_NAMESPACE_USE
using namespace std;
using namespace oracle::occi;
using namespace hcal;

/**

   \class HcalLutManager
   \brief Various manipulations with trigger Lookup Tables
   \author Gena Kukartsev, Brown University, March 14, 2008

*/

HcalLutManager::HcalLutManager(void) { init(); }

HcalLutManager::HcalLutManager(std::vector<HcalGenericDetId>& map) {
  init();
  _iter.init(map);
}

HcalLutManager::HcalLutManager(const HcalElectronicsMap* _emap,
                               const HcalChannelQuality* _cq,
                               uint32_t _status_word_to_mask) {
  init();
  emap = _emap;
  cq = _cq;
  status_word_to_mask = _status_word_to_mask;
}

HcalLutManager::HcalLutManager(const HcalDbService* _conditions,
                               const HcalChannelQuality* _cq,
                               uint32_t _status_word_to_mask) {
  init();
  conditions = _conditions;
  emap = conditions->getHcalMapping();
  cq = _cq;
  status_word_to_mask = _status_word_to_mask;
}

void HcalLutManager::init(void) {
  lut_xml = nullptr;
  lut_checksums_xml = nullptr;
  db = nullptr;
  lmap = nullptr;
  emap = nullptr;
  cq = nullptr;
  conditions = nullptr;
  status_word_to_mask = 0x0000;
}

HcalLutManager::~HcalLutManager(void) {
  delete lut_xml;
  delete lut_checksums_xml;
  delete db;
  delete lmap;
}

int HcalLutManager::initChannelIterator(std::vector<HcalGenericDetId>& map) {
  _iter.init(map);
  return _iter.size();
}

std::string& HcalLutManager::getLutXml(std::vector<unsigned int>& _lut) {
  if (lut_xml)
    delete lut_xml;

  lut_xml = new LutXml();

  LutXml::Config _config;
  _config.lut = _lut;
  lut_xml->addLut(_config);
  lut_xml->addLut(_config);
  lut_xml->addLut(_config);

  //return lut_xml->getString();
  return lut_xml->getCurrentBrick();
}

int HcalLutManager::getInt(std::string number) {
  int result;
  sscanf(number.c_str(), "%d", &result);
  return result;
}

HcalSubdetector HcalLutManager::get_subdetector(std::string _det) {
  HcalSubdetector result;
  if (_det.find("HB") != std::string::npos)
    result = HcalBarrel;
  else if (_det.find("HE") != std::string::npos)
    result = HcalEndcap;
  else if (_det.find("HF") != std::string::npos)
    result = HcalForward;
  else if (_det.find("HO") != std::string::npos)
    result = HcalOuter;
  else
    result = HcalOther;

  return result;
}

HcalLutSet HcalLutManager::getLutSetFromFile(std::string _filename, int _type) {
  HcalLutSet _lutset;

  ifstream infile(_filename.c_str());
  std::string buf;

  if (infile.is_open()) {
    edm::LogInfo("HcalLutManager") << "File " << _filename << " is open..." << std::endl
                                   << "Reading LUTs and their eta/phi/depth/subdet ranges...";

    // get label
    getline(infile, _lutset.label);

    if (_type == 1) {  // for linearization LUTs get subdetectors (default)
      //get subdetectors
      getline(infile, buf);
      _lutset.subdet = HcalQIEManager::splitString(buf);
    }

    //get min etas
    std::vector<std::string> buf_vec;
    getline(infile, buf);
    buf_vec = HcalQIEManager::splitString(buf);
    for (std::vector<std::string>::const_iterator iter = buf_vec.begin(); iter != buf_vec.end(); iter++) {
      _lutset.eta_min.push_back(HcalLutManager::getInt(*iter));
    }

    //get max etas
    getline(infile, buf);
    buf_vec = HcalQIEManager::splitString(buf);
    for (std::vector<std::string>::const_iterator iter = buf_vec.begin(); iter != buf_vec.end(); iter++) {
      _lutset.eta_max.push_back(HcalLutManager::getInt(*iter));
    }

    //get min phis
    getline(infile, buf);
    buf_vec = HcalQIEManager::splitString(buf);
    for (std::vector<std::string>::const_iterator iter = buf_vec.begin(); iter != buf_vec.end(); iter++) {
      _lutset.phi_min.push_back(HcalLutManager::getInt(*iter));
    }

    //get max phis
    getline(infile, buf);
    buf_vec = HcalQIEManager::splitString(buf);
    for (std::vector<std::string>::const_iterator iter = buf_vec.begin(); iter != buf_vec.end(); iter++) {
      _lutset.phi_max.push_back(HcalLutManager::getInt(*iter));
    }

    if (_type == 1) {  // for linearization LUTs get depth range (default)
      //get min depths
      getline(infile, buf);
      buf_vec = HcalQIEManager::splitString(buf);
      for (std::vector<std::string>::const_iterator iter = buf_vec.begin(); iter != buf_vec.end(); iter++) {
        _lutset.depth_min.push_back(HcalLutManager::getInt(*iter));
      }

      //get max depths
      getline(infile, buf);
      buf_vec = HcalQIEManager::splitString(buf);
      for (std::vector<std::string>::const_iterator iter = buf_vec.begin(); iter != buf_vec.end(); iter++) {
        _lutset.depth_max.push_back(HcalLutManager::getInt(*iter));
      }
    }

    bool first_lut_entry = true;
    while (getline(infile, buf)) {
      buf_vec = HcalQIEManager::splitString(buf);
      for (unsigned int i = 0; i < buf_vec.size(); i++) {
        if (first_lut_entry) {
          std::vector<unsigned int> _l;
          _lutset.lut.push_back(_l);
        }
        _lutset.lut[i].push_back(HcalLutManager::getInt(buf_vec[i]));
      }
      first_lut_entry = false;
    }
  }

  edm::LogInfo("HcalLutManager") << "done.";

  return _lutset;
}

//
//_____ get HO from ASCII master here ___________________________________
//
std::map<int, std::shared_ptr<LutXml>> HcalLutManager::getLinearizationLutXmlFromAsciiMasterEmap(std::string _filename,
                                                                                                 std::string _tag,
                                                                                                 int _crate,
                                                                                                 bool split_by_crate) {
  edm::LogInfo("HcalLutManager") << "Generating linearization (input) LUTs from ascii master file...";
  std::map<int, std::shared_ptr<LutXml>> _xml;  // index - crate number

  EMap _emap(emap);
  std::vector<EMap::EMapRow>& _map = _emap.get_map();
  edm::LogInfo("HcalLutManager") << "EMap contains " << _map.size() << " entries";

  // read LUTs and their eta/phi/depth/subdet ranges
  HcalLutSet _set = getLutSetFromFile(_filename);
  int lut_set_size = _set.lut.size();  // number of different luts
  edm::LogInfo("HcalLutManager") << "  ==> " << lut_set_size << " sets of different LUTs read from the master file";

  // setup "zero" LUT for channel masking
  std::vector<unsigned int> zeroLut(128, 0);

  RooGKCounter _counter;
  //loop over all EMap channels
  for (std::vector<EMap::EMapRow>::const_iterator row = _map.begin(); row != _map.end(); row++) {
    if ((row->subdet.find("HB") != string::npos || row->subdet.find("HE") != string::npos ||
         row->subdet.find("HO") != string::npos || row->subdet.find("HF") != string::npos) &&
        row->subdet.size() == 2) {
      LutXml::Config _cfg;

      // search for the correct LUT for a given channel,
      // higher LUT numbers have priority in case of overlapping
      int lut_index = -1;
      for (int i = 0; i < lut_set_size; i++) {
        if ((row->crate == _crate || _crate == -1) &&  // -1 stands for all crates
            _set.eta_min[i] <= row->ieta && _set.eta_max[i] >= row->ieta && _set.phi_min[i] <= row->iphi &&
            _set.phi_max[i] >= row->iphi && _set.depth_min[i] <= row->idepth && _set.depth_max[i] >= row->idepth &&
            _set.subdet[i].find(row->subdet) != string::npos) {
          lut_index = i;
        }
      }
      if (lut_index >= 0) {
        if (_xml.count(row->crate) == 0 && split_by_crate) {
          _xml.insert(std::pair<int, std::shared_ptr<LutXml>>(row->crate, std::make_shared<LutXml>()));
        } else if (_xml.count(0) == 0 && !split_by_crate) {
          _xml.insert(std::pair<int, std::shared_ptr<LutXml>>(0, std::make_shared<LutXml>()));
        }
        _cfg.ieta = row->ieta;
        _cfg.iphi = row->iphi;
        _cfg.depth = row->idepth;
        _cfg.crate = row->crate;
        _cfg.slot = row->slot;
        if (row->topbottom.find('t') != std::string::npos)
          _cfg.topbottom = 1;
        else if (row->topbottom.find('b') != std::string::npos)
          _cfg.topbottom = 0;
        else if (row->topbottom.find('u') != std::string::npos)
          _cfg.topbottom = 2;
        else
          edm::LogWarning("HcalLutManager") << "fpga out of range...";
        _cfg.fiber = row->fiber;
        _cfg.fiberchan = row->fiberchan;
        _cfg.lut_type = 1;
        _cfg.creationtag = _tag;
        _cfg.creationstamp = get_time_stamp(time(nullptr));
        _cfg.targetfirmware = "1.0.0";
        _cfg.formatrevision = "1";  //???
        // "original" definition of GENERALIZEDINDEX from Mike Weinberger
        //    int generalizedIndex=id.ietaAbs()+1000*id.depth()+10000*id.iphi()+
        //        ((id.ieta()<0)?(0):(100))+((id.subdet()==HcalForward && id.ietaAbs()==29)?(4*10000):(0));
        _cfg.generalizedindex =
            _cfg.iphi * 10000 + _cfg.depth * 1000 + (row->ieta > 0) * 100 + abs(row->ieta) +
            (((row->subdet.find("HF") != string::npos) && abs(row->ieta) == 29) ? (4 * 10000) : (0));
        //
        // consider channel status here
        DetId _detId(row->rawId);
        uint32_t status_word = cq->getValues(_detId)->getValue();
        if ((status_word & status_word_to_mask) > 0) {
          _cfg.lut = zeroLut;
        } else {
          _cfg.lut = _set.lut[lut_index];
        }
        if (split_by_crate) {
          _xml[row->crate]->addLut(_cfg, lut_checksums_xml);
          _counter.count();
        } else {
          _xml[0]->addLut(_cfg, lut_checksums_xml);
          _counter.count();
        }
      }
    }
  }
  edm::LogInfo("HcalLutManager") << "LUTs generated: " << _counter.getCount() << std::endl
                                 << "Generating linearization (input) LUTs from ascii master file...DONE" << std::endl;
  return _xml;
}

std::map<int, std::shared_ptr<LutXml>> HcalLutManager::getLinearizationLutXmlFromCoder(const HcalTPGCoder& _coder,
                                                                                       std::string _tag,
                                                                                       bool split_by_crate) {
  edm::LogInfo("HcalLutManager") << "Generating linearization (input) LUTs from HcaluLUTTPGCoder...";
  std::map<int, std::shared_ptr<LutXml>> _xml;  // index - crate number

  //EMap _emap("../../../CondFormats/HcalObjects/data/official_emap_v6.03_080817.txt");
  //std::vector<EMap::EMapRow> & _map = _emap.get_map();
  //std::cout << "EMap contains " << _map . size() << " entries" << std::endl;

  LMap _lmap;
  _lmap.read("backup/HCALmapHBEF.txt", "HBEF");
  // HO is not part of trigger, so TPGCoder cannot generate LUTs for it
  //_lmap . read( "backup/HCALmapHO.txt", "HO" );
  std::map<int, LMapRow>& _map = _lmap.get_map();
  edm::LogInfo("HcalLutManager") << "LMap contains " << _map.size() << " channels";

  // read LUTs and their eta/phi/depth/subdet ranges
  //HcalLutSet _set = getLinearizationLutSetFromCoder();
  //int lut_set_size = _set.lut.size(); // number of different luts

  //loop over all HCAL channels
  RooGKCounter _counter;
  for (std::map<int, LMapRow>::const_iterator row = _map.begin(); row != _map.end(); row++) {
    LutXml::Config _cfg;

    if (_xml.count(row->second.crate) == 0 && split_by_crate) {
      _xml.insert(std::pair<int, std::shared_ptr<LutXml>>(row->second.crate, std::make_shared<LutXml>()));
    } else if (_xml.count(0) == 0 && !split_by_crate) {
      _xml.insert(std::pair<int, std::shared_ptr<LutXml>>(0, std::make_shared<LutXml>()));
    }
    _cfg.ieta = row->second.side * row->second.eta;
    _cfg.iphi = row->second.phi;
    _cfg.depth = row->second.depth;
    _cfg.crate = row->second.crate;
    _cfg.slot = row->second.htr;
    if (row->second.fpga.find("top") != std::string::npos)
      _cfg.topbottom = 1;
    else if (row->second.fpga.find("bot") != std::string::npos)
      _cfg.topbottom = 0;
    else
      edm::LogWarning("HcalLutManager") << "fpga out of range...";
    // FIXME: probably fixed. fiber==htr_fi, not rm_fi in LMAP notation.
    //_cfg.fiber = row->second.rm_fi;
    _cfg.fiber = row->second.htr_fi;
    _cfg.fiberchan = row->second.fi_ch;
    _cfg.lut_type = 1;
    _cfg.creationtag = _tag;
    _cfg.creationstamp = get_time_stamp(time(nullptr));
    _cfg.targetfirmware = "1.0.0";
    _cfg.formatrevision = "1";  //???
    // "original" definition of GENERALIZEDINDEX from Mike Weinberger
    //    int generalizedIndex=id.ietaAbs()+1000*id.depth()+10000*id.iphi()+
    //        ((id.ieta()<0)?(0):(100))+((id.subdet()==HcalForward && id.ietaAbs()==29)?(4*10000):(0));
    _cfg.generalizedindex = _cfg.iphi * 10000 + _cfg.depth * 1000 + (row->second.side > 0) * 100 + row->second.eta +
                            ((row->second.det == HcalForward && row->second.eta == 29) ? (4 * 10000) : (0));

    //HcalDetId _detid(row->first);
    HcalDetId _detid(row->second.det, row->second.side * row->second.eta, row->second.phi, row->second.depth);
    //std::cout << "### DEBUG: rawid = " << _detid.rawId() << std::endl;

    //std::cout << "### DEBUG: subdetector = " << row->second.det << std::endl;
    std::vector<unsigned short> coder_lut = _coder.getLinearizationLUT(_detid);
    for (std::vector<unsigned short>::const_iterator _i = coder_lut.begin(); _i != coder_lut.end(); _i++) {
      unsigned int _temp = (unsigned int)(*_i);
      //if (_temp!=0) std::cout << "DEBUG non-zero LUT!!!!!!!!!!!!!!!" << (*_i) << "     " << _temp << std::endl;
      //unsigned int _temp = 0;
      _cfg.lut.push_back(_temp);
    }
    if (split_by_crate) {
      _xml[row->second.crate]->addLut(_cfg, lut_checksums_xml);
      _counter.count();
    } else {
      _xml[0]->addLut(_cfg, lut_checksums_xml);
      _counter.count();
    }
  }
  edm::LogInfo("HcalLutManager") << "Generated LUTs: " << _counter.getCount() << std::endl
                                 << "Generating linearization (input) LUTs from HcaluLUTTPGCoder...DONE" << std::endl;
  return _xml;
}

std::map<int, std::shared_ptr<LutXml>> HcalLutManager::getMasks(int masktype, std::string _tag, bool split_by_crate) {
  edm::LogInfo("HcalLutManager") << "Generating TDC masks...";

  EMap _emap(emap);
  std::vector<EMap::EMapRow>& _map = _emap.get_map();
  edm::LogInfo("HcalLutManager") << "EMap contains new" << _map.size() << " entries";

  std::map<int, std::vector<uint64_t>> masks;

  for (const auto& row : _map) {
    std::string subdet = row.subdet;
    if (subdet != "HF")
      continue;
    int crate = row.crate;
    int slot = row.slot;
    int crot = 100 * crate + slot;
    int fiber = row.fiber;
    int channel = row.fiberchan;
    unsigned int finel = 4 * fiber + channel;
    if (masks.count(crot) == 0)
      masks[crot] = {};
    if (finel >= masks[crot].size())
      masks[crot].resize(finel + 1);

    if (masktype == 0) {
      HcalSubdetector _subdet;
      if (row.subdet.find("HB") != string::npos)
        _subdet = HcalBarrel;
      else if (row.subdet.find("HE") != string::npos)
        _subdet = HcalEndcap;
      else if (row.subdet.find("HO") != string::npos)
        _subdet = HcalOuter;
      else if (row.subdet.find("HF") != string::npos)
        _subdet = HcalForward;
      else
        _subdet = HcalOther;
      HcalDetId _detid(_subdet, row.ieta, row.iphi, row.idepth);
      masks[crot][finel] = conditions->getHcalTPChannelParameter(_detid)->getMask();
    } else {
      auto parameters = conditions->getHcalTPParameters();
      masks[crot][finel] = masktype == 1 ? parameters->getADCThresholdHF() : parameters->getTDCMaskHF();
    }
  }

  std::map<int, std::shared_ptr<LutXml>> _xml;  // index - crate number
  RooGKCounter _counter;

  for (const auto& i : masks) {
    int crot = i.first;
    int crate = crot / 100;

    LutXml::Config _cfg;
    _cfg.lut_type = 5 + masktype;
    _cfg.crate = crate;
    _cfg.slot = crot % 100;
    _cfg.generalizedindex = crot;
    _cfg.mask = i.second;
    _cfg.creationtag = _tag;
    _cfg.targetfirmware = "1.0.0";
    _cfg.formatrevision = "1";

    int c = split_by_crate ? crate : 0;
    if (_xml.count(c) == 0)
      _xml[c] = std::make_shared<LutXml>();

    _xml[c]->addLut(_cfg);
    _counter.count();
  }

  edm::LogInfo("HcalLutManager") << "Generated LUTs: " << _counter.getCount() << std::endl
                                 << "Generating Masks...DONE" << std::endl;
  return _xml;
}

std::map<int, std::shared_ptr<LutXml>> HcalLutManager::getLinearizationLutXmlFromCoderEmap(const HcalTPGCoder& _coder,
                                                                                           std::string _tag,
                                                                                           bool split_by_crate) {
  edm::LogInfo("HcalLutManager") << "Generating linearization (input) LUTs from HcaluLUTTPGCoder...";
  std::map<int, std::shared_ptr<LutXml>> _xml;  // index - crate number

  EMap _emap(emap);
  std::vector<EMap::EMapRow>& _map = _emap.get_map();
  edm::LogInfo("HcalLutManager") << "EMap contains " << _map.size() << " entries";

  RooGKCounter _counter;
  //loop over all EMap channels
  for (std::vector<EMap::EMapRow>::const_iterator row = _map.begin(); row != _map.end(); row++) {
    if ((row->subdet.find("HB") != string::npos || row->subdet.find("HE") != string::npos ||
         row->subdet.find("HF") != string::npos) &&
        row->subdet.size() == 2) {
      LutXml::Config _cfg;

      if (_xml.count(row->crate) == 0 && split_by_crate) {
        _xml.insert(std::pair<int, std::shared_ptr<LutXml>>(row->crate, std::make_shared<LutXml>()));
      } else if (_xml.count(0) == 0 && !split_by_crate) {
        _xml.insert(std::pair<int, std::shared_ptr<LutXml>>(0, std::make_shared<LutXml>()));
      }
      _cfg.ieta = row->ieta;
      _cfg.iphi = row->iphi;
      _cfg.depth = row->idepth;
      _cfg.crate = row->crate;
      _cfg.slot = row->slot;
      if (row->topbottom.find('t') != std::string::npos)
        _cfg.topbottom = 1;
      else if (row->topbottom.find('b') != std::string::npos)
        _cfg.topbottom = 0;
      else if (row->topbottom.find('u') != std::string::npos)
        _cfg.topbottom = 2;
      else
        edm::LogWarning("HcalLutManager") << "fpga out of range...";
      _cfg.fiber = row->fiber;
      _cfg.fiberchan = row->fiberchan;
      _cfg.lut_type = 1;
      _cfg.creationtag = _tag;
      _cfg.creationstamp = get_time_stamp(time(nullptr));
      _cfg.targetfirmware = "1.0.0";
      _cfg.formatrevision = "1";  //???
      // "original" definition of GENERALIZEDINDEX from Mike Weinberger
      //    int generalizedIndex=id.ietaAbs()+1000*id.depth()+10000*id.iphi()+
      //        ((id.ieta()<0)?(0):(100))+((id.subdet()==HcalForward && id.ietaAbs()==29)?(4*10000):(0));
      _cfg.generalizedindex = _cfg.iphi * 10000 + _cfg.depth * 1000 + (row->ieta > 0) * 100 + abs(row->ieta) +
                              (((row->subdet.find("HF") != string::npos) && abs(row->ieta) == 29) ? (4 * 10000) : (0));
      HcalSubdetector _subdet;
      if (row->subdet.find("HB") != string::npos)
        _subdet = HcalBarrel;
      else if (row->subdet.find("HE") != string::npos)
        _subdet = HcalEndcap;
      else if (row->subdet.find("HO") != string::npos)
        _subdet = HcalOuter;
      else if (row->subdet.find("HF") != string::npos)
        _subdet = HcalForward;
      else
        _subdet = HcalOther;
      HcalDetId _detid(_subdet, row->ieta, row->iphi, row->idepth);

      for (const auto i : _coder.getLinearizationLUT(_detid))
        _cfg.lut.push_back(i);

      if (split_by_crate) {
        _xml[row->crate]->addLut(_cfg, lut_checksums_xml);
        _counter.count();
      } else {
        _xml[0]->addLut(_cfg, lut_checksums_xml);
        _counter.count();
      }
    }
  }
  edm::LogInfo("HcalLutManager") << "Generated LUTs: " << _counter.getCount() << std::endl
                                 << "Generating linearization (input) LUTs from HcaluLUTTPGCoder...DONE" << std::endl;
  return _xml;
}

std::map<int, std::shared_ptr<LutXml>> HcalLutManager::getCompressionLutXmlFromCoder(
    const CaloTPGTranscoderULUT& _coder, std::string _tag, bool split_by_crate) {
  edm::LogInfo("HcalLutManager") << "Generating compression (output) LUTs from CaloTPGTranscoderULUT," << std::endl
                                 << "initialized from Event Setup" << std::endl;
  std::map<int, std::shared_ptr<LutXml>> _xml;  // index - crate number

  EMap _emap(emap);

  std::map<int, unsigned int> maxsize;

  std::vector<EMap::EMapRow>& _map = _emap.get_map();
  edm::LogInfo("HcalLutManager") << "EMap contains " << _map.size() << " channels";

  //need to equalize compression LUT size in each crate-slot, needed for mixed uHTR
  for (const auto& row : _map) {
    if (row.subdet.find("HT") == std::string::npos)
      continue;
    HcalTrigTowerDetId _detid(row.rawId);
    if (!cq->topo()->validHT(_detid))
      continue;
    int crot = 100 * row.crate + row.slot;
    unsigned int size = _coder.getCompressionLUT(_detid).size();
    if (maxsize.count(crot) == 0 || size > maxsize[crot])
      maxsize[crot] = size;
  }

  RooGKCounter _counter;
  for (std::vector<EMap::EMapRow>::const_iterator row = _map.begin(); row != _map.end(); row++) {
    LutXml::Config _cfg;

    if (row->subdet.find("HT") == std::string::npos)
      continue;

    HcalTrigTowerDetId _detid(row->rawId);

    if (!cq->topo()->validHT(_detid))
      continue;

    if (_xml.count(row->crate) == 0 && split_by_crate) {
      _xml.insert(std::pair<int, std::shared_ptr<LutXml>>(row->crate, std::make_shared<LutXml>()));
    } else if (_xml.count(0) == 0 && !split_by_crate) {
      _xml.insert(std::pair<int, std::shared_ptr<LutXml>>(0, std::make_shared<LutXml>()));
    }

    _cfg.ieta = row->ieta;
    _cfg.iphi = row->iphi;
    _cfg.depth = row->idepth;
    _cfg.crate = row->crate;
    _cfg.slot = row->slot;
    if (row->topbottom.find('t') != std::string::npos)
      _cfg.topbottom = 1;
    else if (row->topbottom.find('b') != std::string::npos)
      _cfg.topbottom = 0;
    else if (row->topbottom.find('u') != std::string::npos)
      _cfg.topbottom = 2;
    else
      edm::LogWarning("HcalLutManager") << "fpga out of range...";
    _cfg.fiber = row->fiber;
    _cfg.fiberchan = row->fiberchan;
    _cfg.lut_type = 2;
    _cfg.creationtag = _tag;
    _cfg.creationstamp = get_time_stamp(time(nullptr));
    _cfg.targetfirmware = "1.0.0";
    _cfg.formatrevision = "1";                                                           //???
    _cfg.generalizedindex = _cfg.iphi * 10000 + (row->ieta > 0) * 100 + abs(row->ieta);  //is this used for anything?

    _cfg.lut = _coder.getCompressionLUT(_detid);
    auto pWeight = conditions->getHcalTPChannelParameter(_detid, false);
    if (pWeight) {
      _cfg.weight = pWeight->getauxi1();
      _cfg.codedvetothreshold = pWeight->getauxi2();
    }

    int crot = 100 * row->crate + row->slot;
    unsigned int size = _cfg.lut.size();
    if (size < maxsize[crot]) {
      edm::LogWarning("HcalLutManager") << " resizing LUT for " << _detid << ", channel=[" << _cfg.crate << ":"
                                        << _cfg.slot << ":" << _cfg.fiber << ":" << _cfg.fiberchan
                                        << "], using value=" << _cfg.lut[size - 1] << std::endl;
      for (unsigned int i = size; i < maxsize[crot]; ++i)
        _cfg.lut.push_back(_cfg.lut[size - 1]);
    }

    if (split_by_crate) {
      _xml[row->crate]->addLut(_cfg, lut_checksums_xml);
      _counter.count();
    } else {
      _xml[0]->addLut(_cfg, lut_checksums_xml);
      _counter.count();
    }
  }
  edm::LogInfo("HcalLutManager") << "LUTs generated: " << _counter.getCount() << std::endl
                                 << "Generating compression (output) LUTs from CaloTPGTranscoderULUT...DONE"
                                 << std::endl;

  return _xml;
}

int HcalLutManager::writeLutXmlFiles(std::map<int, std::shared_ptr<LutXml>>& _xml,
                                     std::string _tag,
                                     bool split_by_crate) {
  for (std::map<int, std::shared_ptr<LutXml>>::const_iterator cr = _xml.begin(); cr != _xml.end(); cr++) {
    std::stringstream output_file_name;
    if (split_by_crate) {
      output_file_name << _tag << "_" << cr->first << ".xml";
    } else {
      output_file_name << _tag << ".xml";
    }
    cr->second->write(output_file_name.str());
  }
  return 0;
}

void HcalLutManager::addLutMap(std::map<int, std::shared_ptr<LutXml>>& result,
                               const std::map<int, std::shared_ptr<LutXml>>& other) {
  for (std::map<int, std::shared_ptr<LutXml>>::const_iterator lut = other.begin(); lut != other.end(); lut++) {
    edm::LogInfo("HcalLutManager") << "Added LUTs for crate " << lut->first;
    if (result.count(lut->first) == 0) {
      result.insert(*lut);
    } else {
      *(result[lut->first]) += *(lut->second);
    }
  }
}

string HcalLutManager::get_time_stamp(time_t _time) {
  char timebuf[50];
  //strftime( timebuf, 50, "%c", gmtime( &_time ) );
  strftime(timebuf, 50, "%Y-%m-%d %H:%M:%S", gmtime(&_time));
  std::string creationstamp = timebuf;

  return creationstamp;
}

int HcalLutManager::read_lmap(std::string lmap_hbef_file, std::string lmap_ho_file) {
  delete lmap;
  lmap = new LMap();
  lmap->read(lmap_hbef_file, "HBEF");
  lmap->read(lmap_ho_file, "HO");
  edm::LogInfo("HcalLutManager") << "LMap contains " << lmap->get_map().size()
                                 << " channels (compare to 9072 of all HCAL channels)";
  return 0;
}

int HcalLutManager::read_luts(std::string lut_xml_file) {
  delete db;
  db = new HCALConfigDB();
  db->connect(lut_xml_file);
  return 0;
}

int HcalLutManager::local_connect(std::string lut_xml_file, std::string lmap_hbef_file, std::string lmap_ho_file) {
  read_lmap(lmap_hbef_file, lmap_ho_file);
  read_luts(lut_xml_file);
  return 0;
}

int HcalLutManager::create_lut_loader(std::string file_list,
                                      std::string _prefix,
                                      std::string tag_name,
                                      std::string comment,
                                      std::string version,
                                      int subversion) {
  edm::LogInfo("HcalLutManager") << "Generating XML loader for LUTs...";
  //std::cout << _prefix << "..." << tag_name << std::endl;

  XMLLUTLoader::loaderBaseConfig baseConf;
  XMLLUTLoader::lutDBConfig conf;
  XMLLUTLoader::checksumsDBConfig CSconf;

  baseConf.tag_name = tag_name;
  //baseConf . comment_description = tag_name;
  baseConf.comment_description = comment;
  baseConf.iov_begin = "1";
  baseConf.iov_end = "-1";

  conf.version = version;

  std::stringstream _subversion;
  _subversion << subversion;
  conf.subversion = _subversion.str();

  CSconf.version = conf.version;
  CSconf.subversion = conf.subversion;
  CSconf.trig_prim_lookuptbl_data_file = _prefix + "_checksums.xml.dat";
  CSconf.comment_description = tag_name;

  XMLLUTLoader doc(&baseConf);

  std::vector<int> crate_number;
  std::vector<std::string> file_name = HcalQIEManager::splitString(file_list);
  for (std::vector<std::string>::const_iterator _f = file_name.begin(); _f != file_name.end(); _f++) {
    int crate_begin = _f->rfind("_");
    int crate_end = _f->rfind(".xml.dat");
    crate_number.push_back(getInt(_f->substr(crate_begin + 1, crate_end - crate_begin - 1)));
  }
  //
  //_____ fix due to the new convention: version/subversion combo must be unique for every payload
  //
  char _buf[128];
  time_t _offset = time(nullptr);
  sprintf(_buf, "%d", (uint32_t)_offset);
  conf.version.append(".");
  conf.version.append(_buf);
  CSconf.version = conf.version;
  //
  for (std::vector<std::string>::const_iterator _file = file_name.begin(); _file != file_name.end(); _file++) {
    conf.trig_prim_lookuptbl_data_file = *_file;
    //conf . trig_prim_lookuptbl_data_file += ".dat";
    conf.crate = crate_number[_file - file_name.begin()];
    //
    //_____ fix due to the new convention: version/subversion combo must be unique for every payload
    //
    sprintf(_buf, "%.2d", conf.crate);
    conf.subversion.clear();
    conf.subversion.append(_buf);
    sprintf(_buf, "CRATE%.2d", conf.crate);
    std::string _namelabel;
    _namelabel.append(_buf);
    conf.name_label = _namelabel;
    doc.addLUT(&conf);
  }

  doc.addChecksums(&CSconf);
  //doc . write( _prefix + "_Loader.xml" );
  doc.write(tag_name + "_Loader.xml");

  edm::LogInfo("HcalLutManager") << "Generating XML loader for LUTs... done.";

  return 0;
}

//
//_____ attempt to include ZDC LUTs _____________________________________
//
int HcalLutManager::createLutXmlFiles_HBEFFromCoder_HOFromAscii_ZDC(std::string _tag,
                                                                    const HcalTPGCoder& _coder,
                                                                    const CaloTPGTranscoderULUT& _transcoder,
                                                                    std::string _lin_file,
                                                                    bool split_by_crate) {
  std::map<int, std::shared_ptr<LutXml>> xml;
  if (!lut_checksums_xml) {
    lut_checksums_xml = new XMLDOMBlock("CFGBrick", 1);
  }

  if (!_lin_file.empty()) {
    const std::map<int, std::shared_ptr<LutXml>> _lin_lut_ascii_xml =
        getLinearizationLutXmlFromAsciiMasterEmap(_lin_file, _tag, -1, split_by_crate);
    addLutMap(xml, _lin_lut_ascii_xml);
  }
  const std::map<int, std::shared_ptr<LutXml>> _lin_lut_xml =
      getLinearizationLutXmlFromCoderEmap(_coder, _tag, split_by_crate);
  addLutMap(xml, _lin_lut_xml);
  //
  const std::map<int, std::shared_ptr<LutXml>> _comp_lut_xml =
      getCompressionLutXmlFromCoder(_transcoder, _tag, split_by_crate);
  addLutMap(xml, _comp_lut_xml);

  for (auto masktype : {0, 1, 2}) {
    const auto masks = getMasks(masktype, _tag, split_by_crate);
    addLutMap(xml, masks);
  }
  //
  const auto _zdc_lut_xml = getZdcLutXml(_coder, _tag, split_by_crate, false);
  addLutMap(xml, _zdc_lut_xml);

  const auto _zdc_ootpu_lut_xml = getZdcLutXml(_coder, _tag, split_by_crate, true);
  addLutMap(xml, _zdc_ootpu_lut_xml);

  writeLutXmlFiles(xml, _tag, split_by_crate);

  std::string checksums_file = _tag + "_checksums.xml";
  lut_checksums_xml->write(checksums_file);

  return 0;
}

std::map<int, std::shared_ptr<LutXml>> HcalLutManager::getZdcLutXml(const HcalTPGCoder& _coder,
                                                                    std::string _tag,
                                                                    bool split_by_crate,
                                                                    bool ootpu_lut) {
  edm::LogInfo("HcalLutManager") << "Generating ZDC LUTs ...may the Force be with us...";
  std::map<int, std::shared_ptr<LutXml>> _xml;  // index - crate number

  EMap _emap(emap);

  std::vector<EMap::EMapRow>& _map = _emap.get_map();
  edm::LogInfo("HcalLutManager") << "EMap contains " << _map.size() << " channels";

  const auto lutMetaDataChannels = conditions->getHcalLutMetadata()->getAllChannels();

  //loop over all EMap channels
  RooGKCounter _counter;
  for (std::vector<EMap::EMapRow>::const_iterator row = _map.begin(); row != _map.end(); row++) {
    LutXml::Config _cfg;

    // only ZDC channels
    if (row->zdc_section.find("ZDC") != std::string::npos) {
      if (_xml.count(row->crate) == 0 && split_by_crate) {
        _xml.insert(std::pair<int, std::shared_ptr<LutXml>>(row->crate, std::make_shared<LutXml>()));
      } else if (_xml.count(0) == 0 && !split_by_crate) {
        _xml.insert(std::pair<int, std::shared_ptr<LutXml>>(0, std::make_shared<LutXml>()));
      }
      //  FIXME: introduce proper tag names in ZDC bricks for logical channel info
      _cfg.ieta = row->zdc_channel;  // int
      //_cfg.ieta = row->zdc_zside; // int
      //_cfg.iphi = row->zdc_section; // string
      _cfg.depth = row->zdc_zside;  // int
      _cfg.crate = row->crate;
      _cfg.slot = row->slot;
      if (row->topbottom.find('t') != std::string::npos)
        _cfg.topbottom = 1;
      else if (row->topbottom.find('b') != std::string::npos)
        _cfg.topbottom = 0;
      else if (row->topbottom.find('u') != std::string::npos)
        _cfg.topbottom = 2;
      else
        edm::LogWarning("HcalLutManager") << "fpga out of range...";

      if (ootpu_lut)
        _cfg.fiber = row->fiber + 6;
      else
        _cfg.fiber = row->fiber;

      _cfg.fiberchan = row->fiberchan;
      _cfg.lut_type = 1;
      _cfg.creationtag = _tag;
      _cfg.creationstamp = get_time_stamp(time(nullptr));
      _cfg.targetfirmware = "1.0.0";
      _cfg.formatrevision = "1";  //???
      _cfg.generalizedindex = 0;

      HcalZDCDetId::Section section = HcalZDCDetId::Unknown;
      if (row->zdc_section == "ZDC EM") {
        section = HcalZDCDetId::EM;
        _cfg.iphi = 1;
      } else if (row->zdc_section == "ZDC HAD") {
        section = HcalZDCDetId::HAD;
        _cfg.iphi = 2;
      } else {
        continue;
      }
      HcalZDCDetId _zdcdetid(section, (row->zdc_zside > 0), row->zdc_channel);

      bool isInLutMetadata = false;
      for (const auto& detid : lutMetaDataChannels) {
        if (detid.det() != DetId::Calo or detid.subdetId() != HcalZDCDetId::SubdetectorId)
          continue;

        HcalZDCDetId zdcdetid(detid.rawId());
        if (_zdcdetid == zdcdetid) {
          isInLutMetadata = true;
          break;
        }
      }

      if (!isInLutMetadata)
        continue;

      for (const auto i : _coder.getLinearizationLUT(_zdcdetid, ootpu_lut)) {
        _cfg.lut.push_back(i);
      }

      if (split_by_crate) {
        _xml[row->crate]->addLut(_cfg, lut_checksums_xml);
        _counter.count();
      } else {
        _xml[0]->addLut(_cfg, lut_checksums_xml);
        _counter.count();
      }
      //size of lut
    }
  }

  edm::LogInfo("HcalLutManager") << "LUTs generated: " << _counter.getCount() << std::endl
                                 << "Generating ZDC LUTs...DONE" << std::endl;

  return _xml;
}
