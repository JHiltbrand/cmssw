#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "FWCore/Utilities/interface/HostDeviceConstant.h"

#include "KernelHelpers.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::ecal::reconstruction {

  namespace internal::barrel {

    ALPAKA_FN_ACC ALPAKA_FN_INLINE bool positiveZ(uint32_t id) { return id & 0x10000; }

    ALPAKA_FN_ACC ALPAKA_FN_INLINE uint32_t ietaAbs(uint32_t id) { return (id >> 9) & 0x7F; }

    ALPAKA_FN_ACC ALPAKA_FN_INLINE uint32_t iphi(uint32_t id) { return id & 0x1FF; }

    ALPAKA_FN_ACC int dccFromSm(int ism) {
      int idcc = 9 + ism;
      if (ism > 18)
        idcc -= 18;
      else
        idcc += 18;
      return idcc;
    }

    ALPAKA_FN_ACC int sm(int ieta, int iphi) {
      if (iphi > 360)
        iphi -= 360;
      int ism = (iphi - 1) / 20 + 1;
      if (ieta < 0)
        ism += 18;
      return ism;
    }

    ALPAKA_FN_ACC int dcc(int ieta, int iphi) {
      int const ism = sm(ieta, iphi);
      return dccFromSm(ism);
    }

    ALPAKA_FN_ACC int lm_channel(int iX, int iY) {
      static const int idx_[] = {
          // clang-format off
     // 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
        1, 2, 2, 2, 2, 4, 4, 4, 4, 6, 6, 6, 6, 8, 8, 8, 8,  // 3
        1, 2, 2, 2, 2, 4, 4, 4, 4, 6, 6, 6, 6, 8, 8, 8, 8,  // 2
        1, 3, 3, 3, 3, 5, 5, 5, 5, 7, 7, 7, 7, 9, 9, 9, 9,  // 1
        1, 3, 3, 3, 3, 5, 5, 5, 5, 7, 7, 7, 7, 9, 9, 9, 9  // 0
     // 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
          // clang-format on
      };

      constexpr int iym = 4;
      constexpr int ixm = 17;
      int const il = iym - iY + 1;
      int const ic = iX;
      int const ii = il * ixm + ic;
      if (ii < 0 || ii > (int)(sizeof(idx_) / sizeof(int))) {
        return -1;
      };
      return idx_[ii];
    }

    ALPAKA_FN_ACC int localCoord_x(int ieta) {
      int iz = 1;
      if (ieta < 0) {
        iz = -1;
      }
      ieta *= iz;
      int ix = ieta - 1;

      return ix;
    }

    ALPAKA_FN_ACC int localCoord_y(int ieta, int iphi) {
      if (iphi > 360) {
        iphi -= 360;
      }
      int iy = (iphi - 1) % 20;
      if (ieta < 0) {
        iy = 19 - iy;
      }

      return iy;
    }

    ALPAKA_FN_ACC int lmmod(int ieta, int iphi) {
      int const ix = localCoord_x(ieta);
      int const iy = localCoord_y(ieta, iphi);

      return lm_channel(ix / 5, iy / 5);
    }

    ALPAKA_FN_ACC int side(int ieta, int iphi) {
      int const ilmmod = lmmod(ieta, iphi);
      return (ilmmod % 2 == 0) ? 1 : 0;
    }

  }  // namespace internal::barrel

  ALPAKA_FN_ACC uint32_t hashedIndexEB(uint32_t id) {
    using namespace internal::barrel;
    return (EBDetId::MAX_IETA + (positiveZ(id) ? ietaAbs(id) - 1 : -ietaAbs(id))) * EBDetId::MAX_IPHI + iphi(id) - 1;
  }

  //
  // https://cmssdt.cern.ch/lxr/source/CalibCalorimetry/EcalLaserAnalyzer/src/MEEBGeom.cc
  //  function: "lmr"

  ALPAKA_FN_ACC int32_t laserMonitoringRegionEB(uint32_t id) {
    using namespace internal::barrel;

    int ieta;
    if (positiveZ(id)) {
      ieta = ietaAbs(id);
    } else {
      ieta = -ietaAbs(id);
    }

    int const idcc = dcc(ieta, (int)(iphi(id)));
    int const ism = idcc - 9;

    int const iside = side(ieta, (int)(iphi(id)));

    return (1 + 2 * (ism - 1) + iside);
  }

  namespace internal::endcap {

    ALPAKA_FN_ACC ALPAKA_FN_INLINE uint32_t ix(uint32_t id) { return (id >> 7) & 0x7F; }

    ALPAKA_FN_ACC ALPAKA_FN_INLINE uint32_t iy(uint32_t id) { return id & 0x7F; }

    ALPAKA_FN_ACC ALPAKA_FN_INLINE bool positiveZ(uint32_t id) { return id & 0x4000; }

    // these constants come from EE Det Id
    HOST_DEVICE_CONSTANT unsigned short kxf[] = {
        41, 51, 41, 51, 41, 51, 36, 51, 36, 51, 26, 51, 26, 51, 26, 51, 21, 51, 21, 51, 21, 51, 21, 51, 21,
        51, 16, 51, 16, 51, 14, 51, 14, 51, 14, 51, 14, 51, 14, 51, 9,  51, 9,  51, 9,  51, 9,  51, 9,  51,
        6,  51, 6,  51, 6,  51, 6,  51, 6,  51, 6,  51, 6,  51, 6,  51, 6,  51, 6,  51, 4,  51, 4,  51, 4,
        51, 4,  51, 4,  56, 1,  58, 1,  59, 1,  60, 1,  61, 1,  61, 1,  62, 1,  62, 1,  62, 1,  62, 1,  62,
        1,  62, 1,  62, 1,  62, 1,  62, 1,  62, 1,  61, 1,  61, 1,  60, 1,  59, 1,  58, 4,  56, 4,  51, 4,
        51, 4,  51, 4,  51, 6,  51, 6,  51, 6,  51, 6,  51, 6,  51, 6,  51, 6,  51, 6,  51, 6,  51, 6,  51,
        9,  51, 9,  51, 9,  51, 9,  51, 9,  51, 14, 51, 14, 51, 14, 51, 14, 51, 14, 51, 16, 51, 16, 51, 21,
        51, 21, 51, 21, 51, 21, 51, 21, 51, 26, 51, 26, 51, 26, 51, 36, 51, 36, 51, 41, 51, 41, 51, 41, 51};

    HOST_DEVICE_CONSTANT unsigned short kdi[] = {
        0,    10,   20,   30,   40,   50,   60,   75,   90,   105,  120,  145,  170,  195,  220,  245,  270,
        300,  330,  360,  390,  420,  450,  480,  510,  540,  570,  605,  640,  675,  710,  747,  784,  821,
        858,  895,  932,  969,  1006, 1043, 1080, 1122, 1164, 1206, 1248, 1290, 1332, 1374, 1416, 1458, 1500,
        1545, 1590, 1635, 1680, 1725, 1770, 1815, 1860, 1905, 1950, 1995, 2040, 2085, 2130, 2175, 2220, 2265,
        2310, 2355, 2400, 2447, 2494, 2541, 2588, 2635, 2682, 2729, 2776, 2818, 2860, 2903, 2946, 2988, 3030,
        3071, 3112, 3152, 3192, 3232, 3272, 3311, 3350, 3389, 3428, 3467, 3506, 3545, 3584, 3623, 3662, 3701,
        3740, 3779, 3818, 3857, 3896, 3935, 3974, 4013, 4052, 4092, 4132, 4172, 4212, 4253, 4294, 4336, 4378,
        4421, 4464, 4506, 4548, 4595, 4642, 4689, 4736, 4783, 4830, 4877, 4924, 4969, 5014, 5059, 5104, 5149,
        5194, 5239, 5284, 5329, 5374, 5419, 5464, 5509, 5554, 5599, 5644, 5689, 5734, 5779, 5824, 5866, 5908,
        5950, 5992, 6034, 6076, 6118, 6160, 6202, 6244, 6281, 6318, 6355, 6392, 6429, 6466, 6503, 6540, 6577,
        6614, 6649, 6684, 6719, 6754, 6784, 6814, 6844, 6874, 6904, 6934, 6964, 6994, 7024, 7054, 7079, 7104,
        7129, 7154, 7179, 7204, 7219, 7234, 7249, 7264, 7274, 7284, 7294, 7304, 7314};

    ALPAKA_FN_ACC int quadrant(int iX, int iY) {
      bool const near = iX >= 11;
      bool const far = !near;
      bool const top = iY >= 11;
      bool const bot = !top;

      int iquad = 0;
      if (near && top)
        iquad = 1;
      else if (far && top)
        iquad = 2;
      else if (far && bot)
        iquad = 3;
      else
        iquad = 4;

      return iquad;
    }

    ALPAKA_FN_ACC int sector(int iX, int iY) {
      //  Y (towards the surface)
      //  T
      //  |
      //  |
      //  |
      //  o---------| X  (towards center of LHC)
      //
      static const int idx_[] = {
          // clang-format off
         // 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
            0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 9, 9, 9, 0, 0, 0, 0, 0, 0, 0,  // 20
            0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 9, 9, 9, 9, 9, 9, 0, 0, 0, 0,  // 19
            0, 0, 0, 2, 1, 1, 1, 1, 1, 1, 9, 9, 9, 9, 9, 9, 8, 0, 0, 0,  // 18
            0, 0, 2, 2, 2, 1, 1, 1, 1, 1, 9, 9, 9, 9, 9, 8, 8, 8, 0, 0,  // 17
            0, 2, 2, 2, 2, 1, 1, 1, 1, 1, 9, 9, 9, 9, 9, 8, 8, 8, 8, 0,  // 16
            0, 2, 2, 2, 2, 2, 1, 1, 1, 1, 9, 9, 9, 9, 8, 8, 8, 8, 8, 0,  // 15
            0, 2, 2, 2, 2, 2, 2, 1, 1, 1, 9, 9, 9, 8, 8, 8, 8, 8, 8, 0,  // 14
            2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 9, 9, 8, 8, 8, 8, 8, 8, 8, 8,  // 13
            3, 3, 2, 2, 2, 2, 2, 2, 2, 0, 0, 8, 8, 8, 8, 8, 8, 8, 7, 7,  // 12
            3, 3, 3, 3, 3, 3, 3, 2, 0, 0, 0, 0, 8, 7, 7, 7, 7, 7, 7, 7,  // 11
            3, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, 7, 7, 7, 7, 7, 7, 7, 7,  // 10
            3, 3, 3, 3, 3, 3, 3, 4, 4, 0, 0, 6, 6, 7, 7, 7, 7, 7, 7, 7,  // 9
            3, 3, 3, 3, 3, 3, 4, 4, 4, 5, 5, 6, 6, 6, 7, 7, 7, 7, 7, 7,  // 8
            0, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 6, 6, 6, 6, 6, 7, 7, 7, 0,  // 7
            0, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 0,  // 6
            0, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 0,  // 5
            0, 0, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 0, 0,  // 4
            0, 0, 0, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 0, 0, 0,  // 3
            0, 0, 0, 0, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 0, 0, 0, 0,  // 2
            0, 0, 0, 0, 0, 0, 0, 5, 5, 5, 5, 5, 5, 0, 0, 0, 0, 0, 0, 0   // 1
         // 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
          // clang-format on
      };

      constexpr int iym = 20;
      constexpr int ixm = 20;
      int const il = iym - iY;
      int const ic = iX - 1;
      int const ii = il * ixm + ic;

      if (ii < 0 || ii > (int)(sizeof(idx_) / sizeof(int)) || idx_[ii] == 0) {
        return -1;
      };
      return idx_[ii];
    }

  }  // namespace internal::endcap

  ALPAKA_FN_ACC uint32_t hashedIndexEE(uint32_t id) {
    using namespace internal::endcap;

    const uint32_t jx(ix(id));
    const uint32_t jd(2 * (iy(id) - 1) + (jx - 1) / 50);
    return ((positiveZ(id) ? EEDetId::kEEhalf : 0) + kdi[jd] + jx - kxf[jd]);
  }

  //
  // https://cmssdt.cern.ch/lxr/source/CalibCalorimetry/EcalLaserAnalyzer/src/MEEEGeom.cc
  // https://github.com/cms-sw/cmssw/blob/master/CalibCalorimetry/EcalLaserCorrection/src/EcalLaserDbService.cc
  //

  ALPAKA_FN_ACC int32_t laserMonitoringRegionEE(uint32_t id) {
    using namespace internal::endcap;

    // SuperCrysCoord
    uint32_t const iX = (ix(id) - 1) / 5 + 1;
    uint32_t const iY = (iy(id) - 1) / 5 + 1;

    // Correct convention
    //   * @param iz iz/zside index: -1 for EE-, +1 for EE+
    //   https://github.com/cms-sw/cmssw/blob/master/DataFormats/EcalDetId/interface/EEDetId.h#L68-L71
    //   zside in https://github.com/cms-sw/cmssw/blob/master/CalibCalorimetry/EcalLaserCorrection/src/EcalLaserDbService.cc#L63
    //
    int const iz = positiveZ(id) ? 1 : -1;

    int const iquad = quadrant(iX, iY);
    int const isect = sector(iX, iY);
    if (isect < 0)
      return -1;

    int ilmr = 0;
    ilmr = isect - 6;
    if (ilmr <= 0)
      ilmr += 9;
    if (ilmr == 9)
      ilmr++;
    else if (ilmr == 8 && iquad == 4)
      ilmr++;
    if (iz == +1)
      ilmr += 72;
    else
      ilmr += 82;

    return ilmr;
  }

  ALPAKA_FN_ACC int32_t rechitSetMasked(uint32_t value, uint32_t x, uint32_t offset, uint32_t width) {
    const uint32_t mask = ((1 << width) - 1) << offset;
    value &= ~mask;
    value |= (x & ((1U << width) - 1)) << offset;
    return value;
  }

  HOST_DEVICE_CONSTANT float p10[] = {1.e-2f, 1.e-1f, 1.f, 1.e1f, 1.e2f, 1.e3f, 1.e4f, 1.e5f, 1.e6f};

  ALPAKA_FN_ACC int32_t rechitGetPower10(float e) {
    int b = e < p10[4] ? 0 : 5;
    for (; b < 9; ++b)
      if (e < p10[b])
        break;
    return b;
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::ecal::reconstruction
