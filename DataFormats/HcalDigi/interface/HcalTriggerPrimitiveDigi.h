#ifndef HCALTRIGGERPRIMITIVEDIGI_H
#define HCALTRIGGERPRIMITIVEDIGI_H 1

#include <ostream>
#include <vector>
#include <array>
#include "CalibFormats/CaloObjects/interface/IntegerCaloSamples.h"
#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"
#include "DataFormats/HcalDigi/interface/HcalTriggerPrimitiveSample.h"

/** \class HcalTriggerPrimitiveDigi
    
\author J. Mans - Minnesota
*/
class HcalTriggerPrimitiveDigi {
public:
  typedef HcalTrigTowerDetId key_type;  ///< For the sorted collection

  HcalTriggerPrimitiveDigi();  // for persistence
  explicit HcalTriggerPrimitiveDigi(const HcalTrigTowerDetId& id);

  static const int MAXSAMPLES = 10;
  static const int TPSAMPLES = 4;

  const HcalTrigTowerDetId& id() const { return id_; }
  int size() const { return (size_ & 0xF); }
  int presamples() const { return hcalPresamples_ & 0xF; }

  /// was ZS MarkAndPass?
  bool zsMarkAndPass() const { return (hcalPresamples_ & 0x10); }
  /// was ZS unsuppressed?
  bool zsUnsuppressed() const { return (hcalPresamples_ & 0x20); }

  void setZSInfo(bool unsuppressed, bool markAndPass);

  const HcalTriggerPrimitiveSample& operator[](int i) const { return data_[i]; }
  const HcalTriggerPrimitiveSample& sample(int i) const { return data_[i]; }

  /// Full "Sample of Interest"
  const HcalTriggerPrimitiveSample& t0() const { return data_[presamples()]; }
  /// Fine-grain bit for the "Sample of Interest"
  bool SOI_fineGrain(int i = 0) const { return t0().fineGrain(i); }
  /// Compressed ET for the "Sample of Interest"
  int SOI_compressedEt() const { return t0().compressedEt(); }

  void setSize(int size);
  void setPresamples(int ps);
  void setSample(int i, const HcalTriggerPrimitiveSample& sam) { data_[i] = sam; }
  void setInputLinearFrame(const IntegerCaloSamples& inputLinearFrame);
  void setOutputLinearFrame(const IntegerCaloSamples& outputLinearFrame);
  void setVetoedTPs(const std::vector<uint8_t>& vetoedTPs);

  const std::array<uint16_t, MAXSAMPLES>& getInputLinearFrame() const { return inputLinearFrame_; }
  const std::array<uint16_t, TPSAMPLES>& getOutputLinearFrame() const { return outputLinearFrame_; }
  const std::array<uint8_t, TPSAMPLES>& getVetoedTPs() const { return vetoedTPs_; }

private:
  HcalTrigTowerDetId id_;
  int size_;
  int hcalPresamples_;
  HcalTriggerPrimitiveSample data_[MAXSAMPLES];

  std::array<uint16_t, MAXSAMPLES> inputLinearFrame_{};
  std::array<uint16_t, TPSAMPLES> outputLinearFrame_{};
  std::array<uint8_t, TPSAMPLES> vetoedTPs_{};
};

std::ostream& operator<<(std::ostream& s, const HcalTriggerPrimitiveDigi& digi);

#endif
