#ifndef CALIBFORMATS_HCALOBJECTS_HCALTPGCODER_H
#define CALIBFORMATS_HCALOBJECTS_HCALTPGCODER_H 1

#include "CalibFormats/CaloObjects/interface/IntegerCaloSamples.h"
#include "DataFormats/HcalDetId/interface/HcalZDCDetId.h"
#include "DataFormats/HcalDigi/interface/HcalTriggerPrimitiveDigi.h"
#include "DataFormats/HcalDigi/interface/QIE10DataFrame.h"
#include "DataFormats/HcalDigi/interface/QIE11DataFrame.h"

// forward declaration of EventSetup is all that is needed here
namespace edm {
  class EventSetup;
}

/** \class HcalTPGCoder
  *  
  * Converts ADC to linear E or ET for use in the TPG path
  * Also compresses linear scale for transmission to RCT
  * 
  * Note : whether the coder produces E or ET is determined by the specific
  * implementation of the coder.
  *
  * \author J. Mans - Minnesota
  */
class HcalTPGCoder {
public:
  virtual ~HcalTPGCoder() = default;
  virtual void adc2Linear(const QIE10DataFrame& df, IntegerCaloSamples& ics, bool ootpu_lut) const = 0;
  virtual void adc2Linear(const QIE11DataFrame& df, IntegerCaloSamples& ics) const = 0;
  virtual void compress(const IntegerCaloSamples& ics,
                        const std::vector<bool>& featureBits,
                        HcalTriggerPrimitiveDigi& tp) const = 0;

  virtual std::vector<unsigned short> getLinearizationLUT(HcalDetId id) const = 0;
  virtual std::vector<unsigned short> getLinearizationLUT(HcalZDCDetId id, bool ootput_lut) const = 0;
};

#endif
