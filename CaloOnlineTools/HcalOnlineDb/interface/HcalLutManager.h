#ifndef HcalLutManager_h
#define HcalLutManager_h

/**

   \class HcalLutManager
   \brief Various manipulations with trigger Lookup Tables
   \author Gena Kukartsev, Brown University, March 14, 2008

*/

#include "CalibCalorimetry/HcalTPGAlgos/interface/LutXml.h"
#include "CalibCalorimetry/CaloTPG/interface/CaloTPGTranscoderULUT.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalTPGCoder.h"
#include "CaloOnlineTools/HcalOnlineDb/interface/HCALConfigDB.h"
#include "CaloOnlineTools/HcalOnlineDb/interface/HcalAssistant.h"
#include "CaloOnlineTools/HcalOnlineDb/interface/HcalChannelIterator.h"
#include "CaloOnlineTools/HcalOnlineDb/interface/LMap.h"
#include "CondFormats/HcalObjects/interface/AllObjects.h"
#include "DataFormats/HcalDetId/interface/HcalGenericDetId.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "SimCalorimetry/HcalTrigPrimAlgos/interface/HcalFinegrainBit.h"

#include <iostream>
#include <map>
#include <string>
#include <vector>

class XMLDOMBlock;

class HcalLutSet {
public:
  std::string label;
  std::vector<std::string> subdet;
  std::vector<int> eta_min, eta_max, phi_min, phi_max, depth_min, depth_max;
  std::vector<std::vector<unsigned int> > lut;
};

class HcalLutManager {
public:
  HcalLutManager();
  HcalLutManager(std::vector<HcalGenericDetId>& map);
  HcalLutManager(const HcalElectronicsMap* _emap,
                 const HcalChannelQuality* _cq = nullptr,
                 uint32_t _status_word_to_mask = 0x0000);

  HcalLutManager(const HcalDbService* conditions,
                 const HcalChannelQuality* _cq = nullptr,
                 uint32_t _status_word_to_mask = 0x0000);

  ~HcalLutManager();

  void init(void);
  std::string& getLutXml(std::vector<unsigned int>& _lut);

  std::map<int, std::shared_ptr<LutXml> > getLinearizationLutXmlFromAsciiMasterEmap(std::string _filename,
                                                                                    std::string _tag,
                                                                                    int _crate,
                                                                                    bool split_by_crate = true);

  std::map<int, std::shared_ptr<LutXml> > getLinearizationLutXmlFromCoder(const HcalTPGCoder& _coder,
                                                                          std::string _tag,
                                                                          bool split_by_crate = true);

  std::map<int, std::shared_ptr<LutXml> > getMasks(int var, std::string _tag, bool split_by_crate = true);

  std::map<int, std::shared_ptr<LutXml> > getLinearizationLutXmlFromCoderEmap(const HcalTPGCoder& _coder,
                                                                              std::string _tag,
                                                                              bool split_by_crate = true);

  std::map<int, std::shared_ptr<LutXml> > getCompressionLutXmlFromCoder(const CaloTPGTranscoderULUT& _coder,
                                                                        std::string _tag,
                                                                        bool split_by_crate = true);

  std::map<int, std::shared_ptr<LutXml> > getZdcLutXml(const HcalTPGCoder& _coder,
                                                       std::string _tag,
                                                       bool split_by_crate = true,
                                                       bool ootpu_lut = false);

  // add two std::map<s with LUTs. Designed mainly for joining compression LUTs to linearization ones.
  void addLutMap(std::map<int, std::shared_ptr<LutXml> >& result, const std::map<int, std::shared_ptr<LutXml> >& other);

  // read LUTs from ASCII master file.
  HcalLutSet getLutSetFromFile(std::string _filename, int _type = 1);  // _type = 1 - linearization, 2 - compression

  int writeLutXmlFiles(std::map<int, std::shared_ptr<LutXml> >& _xml,
                       std::string _tag = "default_tag",
                       bool split_by_crate = true);

  int createLutXmlFiles_HBEFFromCoder_HOFromAscii_ZDC(std::string _tag,
                                                      const HcalTPGCoder& _coder,
                                                      const CaloTPGTranscoderULUT& _transcoder,
                                                      std::string _lin_file,
                                                      bool split_by_crate = true);

  // connect to local XML file with LUTs and local ASCII file with LMAP
  // connection interface through protected members db and lmap
  int read_lmap(std::string lmap_hbef_file, std::string lmap_ho_file);
  int read_luts(std::string lut_xml_file);
  int local_connect(std::string lut_xml_file, std::string lmap_hbef_file, std::string lmap_ho_file);

  int create_lut_loader(std::string file_list,
                        std::string _prefix,
                        std::string tag_name,
                        std::string comment = "default comment",
                        std::string version = "V00-01-01",
                        int subversion = 1);

  // get md5 checksums for LUTs
  std::string get_checksum(std::vector<unsigned int>& lut);

  static int getInt(std::string number);
  static HcalSubdetector get_subdetector(std::string _subdet);
  static std::string get_time_stamp(time_t _time);

  // gives the iterator a list of channels
  int initChannelIterator(std::vector<HcalGenericDetId>& map);

protected:
  LutXml* lut_xml;
  XMLDOMBlock* lut_checksums_xml;
  HCALConfigDB* db;
  LMap* lmap;
  HcalChannelIterator _iter;
  HcalAssistant _ass;
  const HcalElectronicsMap* emap;
  const HcalChannelQuality* cq;
  const HcalDbService* conditions;
  uint32_t status_word_to_mask;
};

#endif
