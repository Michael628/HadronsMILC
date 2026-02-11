#pragma once

#include <array>
#include <iostream>

NAMESPACE_BEGIN(Grid)

// Spin taste parameters for modules that require StagGamma objects.
struct SpinTasteParams : Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(SpinTasteParams, std::string, gammas,
                                  std::string, gauge, bool, applyG5);
  SpinTasteParams(void) : gammas(""), gauge(""), applyG5(false) {}
};

class StagGamma {
public:
  // XYZT convention
  /*    GRID_SERIALIZABLE_ENUM(StagAlgebra, undef,
                             G1     , 0,
                             GT     , 1,
                             GZ     , 2,
                             GZT    , 3,
                             GY     , 4,
                             GYT    , 5,
                             GYZ    , 6,
                             G5X    , 7,
                             GX     , 8,
                             GXT    , 9,
                             GZX    , 10,
                             G5Y    , 11,
                             GXY    , 12,
                             G5Z    , 13,
                             G5T    , 14,
                             G5     , 15);*/

  // TXYZ convention
  GRID_SERIALIZABLE_ENUM(StagAlgebra, undef, G1, 0, GZ, 1, GY, 2, GYZ, 3, GX, 4,
                         GZX, 5, GXY, 6, G5T, 7, GT, 8, GZT, 9, GYT, 10, G5X,
                         11, GXT, 12, G5Y, 13, G5Z, 14, G5, 15);

  typedef std::pair<StagAlgebra, StagAlgebra> SpinTastePair;

public:
  StagGamma() : _spin(0), _taste(0) {}

  StagGamma(StagAlgebra spin, StagAlgebra taste) {
    _spin = spin;
    _taste = taste;
    calculatePhase();
  }

  StagGamma(SpinTastePair initg) { StagGamma(initg.first, initg.second); }

  void setGaugeField(LatticeGaugeField &U_) { U = &U_; }

  inline void setSpin(StagAlgebra g) {
    _spin = g;
    calculatePhase();
  }

  inline void setTaste(StagAlgebra g) {
    _taste = g;
    calculatePhase();
  }

  inline void setSpinTaste(StagAlgebra spin, StagAlgebra taste) {
    _spin = spin;
    _taste = taste;
    calculatePhase();
  }

  inline void setSpinTaste(SpinTastePair g) { setSpinTaste(g.first, g.second); }

  static std::vector<StagGamma::SpinTastePair>
  ParseSpinTasteString(std::string str, bool applyG5 = false) {
    auto gammas = strToVec<StagGamma::SpinTastePair>(str);

    if (applyG5) {
      StagGamma st;
      StagGamma g5(StagGamma::StagAlgebra::G5, StagGamma::StagAlgebra::G5);

      for (auto &g : gammas) {
        st.setSpinTaste(g);
        st = st * g5;
        g.first = st._spin;
        g.second = st._taste;
      }
    }

    return gammas;
  }
  static std::string GetName(StagAlgebra spin, StagAlgebra taste) {

    std::string name = StagGamma::name[spin];
    name = (name + "_") + StagGamma::name[taste];

    return name;
  }

  static std::string GetName(SpinTastePair g) {
    return StagGamma::GetName(g.first, g.second);
  }

  std::string getName() const { return StagGamma::GetName(_spin, _taste); }

  template <typename obj>
  void applyGamma(Lattice<obj> &lhs, const Lattice<obj> &rhs) const;

  template <typename obj>
  void applyPhase(Lattice<obj> &lhs, const Lattice<obj> &rhs) const;

  template <typename obj>
  void oneLink(Lattice<obj> &lhs, const Lattice<obj> &rhs, int shift_dir) const;

  template <typename obj>
  inline void operator()(Lattice<obj> &lhs, const Lattice<obj> &rhs) const {
    applyGamma(lhs, rhs);
  }

private:
  // Calculate the < and > operations as defined in Follana (2007) eqns A5 and
  // A7.
  static inline StagAlgebra LessThan(StagAlgebra g);
  static inline StagAlgebra GreaterThan(StagAlgebra g);

  // Assign negative orientations to StagAlgebra gammas according to txyz (or
  // xyzt?) oriented euclidean space.
  inline int getOrientation(StagAlgebra g);

  // Implements eqn. E3 of Follana (2007)
  inline void calculatePhase();

  // Implements (-1)^(x[mu] * ( _taste^< + _spin^> ) ( see eqn. E3 of Follana
  // (2007) )
  inline void calculateOscillation();

  // Implements (-1)^(_spin * (_spin + _taste)^<) ( see eqn. E3 of Follana
  // (2007) )
  inline void calculateNegation();

  inline void toggleNegation() { _negated = !_negated; }

public:
  static constexpr unsigned int nGamma = 16;
  // TXYZ convention
  static inline const std::array<const char *, nGamma> name = {
      {"G1", "GZ", "GY", "GYZ", "GX", "GZX", "GXY", "G5T", "GT", "GZT", "GYT",
       "G5X", "GXT", "G5Y", "G5Z", "G5"}};
  static inline const std::array<const StagAlgebra, 4> gmu = {{
      StagGamma::StagAlgebra::GX,
      StagGamma::StagAlgebra::GY,
      StagGamma::StagAlgebra::GZ,
      StagGamma::StagAlgebra::GT,
  }};
  friend inline StagGamma operator*(const StagGamma &g1, const StagGamma &g2);

public:
  StagAlgebra _spin, _taste;
  LatticeGaugeField *U = nullptr;

private:
  StagAlgebra _oscillateDirs = 0;
  bool _negated = false;
  RealD _scaling;
};

inline StagGamma::StagAlgebra StagGamma::LessThan(StagAlgebra g) {
  uint8_t ret = 0;
  uint8_t mask = g;

  for (int i = 0; i < Nd - 1; i++) {
    mask = mask >> 1;
    ret = ret ^ mask; // each bit will toggle for each 1 that passes over it
  }
  return ret;
}

inline StagGamma::StagAlgebra StagGamma::GreaterThan(StagAlgebra g) {
  uint8_t ret = 0;
  uint8_t mask = g;

  for (int i = 0; i < Nd - 1; i++) {
    mask = mask << 1;
    ret = ret ^ mask; // each bit will toggle for each 1 that passes over it
  }
  return (0x0F & ret);
}

template <class obj>
void StagGamma::applyGamma(Lattice<obj> &lhs, const Lattice<obj> &rhs) const {
  uint16_t shift = _spin ^ _taste;

  // Dir index according to Grid convention, XYZT
  int dir = 0;

  if (shift != 0) {
    assert(U != nullptr);
  }

  switch (shift) {
  case StagAlgebra::G1:
    applyPhase(lhs, rhs);
    break;
  case StagAlgebra::GT:
    dir++;
  case StagAlgebra::GZ:
    dir++;
  case StagAlgebra::GY:
    dir++;
  case StagAlgebra::GX:
    oneLink(lhs, rhs, dir);
    applyPhase(lhs, lhs);
    break;
  default:
    assert(0);
  }
}

template <class obj>
void StagGamma::oneLink(Lattice<obj> &lhs, const Lattice<obj> &rhs,
                        int shift_dir) const {

  Lattice<obj> temp(rhs.Grid());
  LatticeColourMatrix Umu(rhs.Grid());

  if (rhs.Grid()->_isCheckerBoarded) {
    LatticeColourMatrix Umu_full(U->Grid());
    Umu_full = PeekIndex<LorentzIndex>(*U, shift_dir);
    pickCheckerboard(rhs.Checkerboard(), Umu, Umu_full);
    temp = adj(Umu) * rhs;
    pickCheckerboard(lhs.Checkerboard(), Umu, Umu_full);
  } else {
    Umu = PeekIndex<LorentzIndex>(*U, shift_dir);
    temp = adj(Umu) * rhs;
  }
  lhs = Cshift(temp, shift_dir, -1);
  temp = Cshift(rhs, shift_dir, 1);
  lhs += Umu * temp;
}

inline int StagGamma::getOrientation(StagAlgebra g) {
  switch (g) {
    // XYZT convention
  case StagAlgebra::GZX:
  // case StagAlgebra::G5X:
  // case StagAlgebra::G5Z:
  case StagAlgebra::G5Y:
  case StagAlgebra::G5T:
  case StagAlgebra::G5:

    // TXYZ convention
    // case StagAlgebra::GZX:
    // case StagAlgebra::GXT:
    // case StagAlgebra::GYT:
    // case StagAlgebra::GZT:
    // case StagAlgebra::G5Y:
    // case StagAlgebra::G5T:
    return -1;
    break;
  }
  return 1;
}

inline void StagGamma::calculateNegation() {
  StagAlgebra result = _spin & LessThan(_spin ^ _taste);

  for (auto &dir : StagGamma::gmu) {
    if (dir & result) {
      toggleNegation();
    }
  }
}

inline void StagGamma::calculateOscillation() {

  _oscillateDirs = LessThan(_taste) ^ GreaterThan(_spin);
}

inline void StagGamma::calculatePhase() {

  // scale down according to number of terms in symmetric shift
  switch (_spin ^ _taste) {
  case StagAlgebra::G1:
    _scaling = 1.0;
    break;
  case StagAlgebra::GT:
  case StagAlgebra::GZ:
  case StagAlgebra::GY:
  case StagAlgebra::GX:
    _scaling = 0.5;
    break;
  default:
    assert(0);
  }

  _negated = false;
  calculateOscillation();
  calculateNegation();

  // Include sign flip for consistent orientation of gammas ( see eqn. A4 of
  // Follana (2007) ) if (getOrientation(_spin) != getOrientation(_taste)) {
  // toggleNegation();
  // }
}

template <class obj>
void StagGamma::applyPhase(Lattice<obj> &lhs, const Lattice<obj> &rhs) const {

  GridBase *grid = lhs.Grid();

  Lattice<obj> temp(grid);
  Lattice<iScalar<vInteger>> coor(grid), stag_dirs(grid);
  iScalar<vInteger> one = 1;

  if (_negated) {
    stag_dirs = one;
  } else {
    stag_dirs = Zero();
  }

  for (int dir = 0; dir < gmu.size(); dir++) {
    if (gmu[dir] & _oscillateDirs) { // gmu[dir] maps Grid XYZT convention to
                                     // our current binary convention
      LatticeCoordinate(coor, dir);
      stag_dirs += coor;
    }
  }

  temp = where(mod(stag_dirs, 2) == 0, _scaling * rhs, -_scaling * rhs);

  lhs = std::move(temp);
}

inline StagGamma operator*(const StagGamma &g1, const StagGamma &g2) {

  StagGamma ret(g1._spin ^ g2._spin, g1._taste ^ g2._taste);

  if (g1.U != nullptr) {
    assert(g2.U == g1.U);
    ret.setGaugeField(*(g1.U));
  }

  if (g1._negated != g2._negated) {
    ret.toggleNegation();
  }

  // Following eqn. A4 of Follana (2007)
  uint8_t negate = ((g1._spin & StagGamma::LessThan(g2._spin)) ^
                    (g1._taste & StagGamma::LessThan(g2._taste)));

  for (auto &dir : StagGamma::gmu) {
    if (dir & negate) {
      ret.toggleNegation();
    }
  }

  return ret;
}

template <class obj>
inline Lattice<obj> operator*(const StagGamma &g1, const Lattice<obj> &lat) {

  Lattice<obj> temp(lat.Grid());
  g1.applyGamma(temp, lat);
  return temp;
}

template <class obj>
inline Lattice<obj> operator*(const Lattice<obj> &lat, const StagGamma &g1) {
  return g1 * lat;
}

// Array used to keep Grid ordering, XYZT

// XYZT convention
/*const std::array<const char *, StagGamma::nGamma> StagGamma::name = {{
     "G1" ,
     "GT" ,
     "GZ" ,
     "GZT",
     "GY" ,
     "GYT",
     "GYZ",
     "G5X",
     "GX" ,
     "GXT",
     "GZX",
     "G5Y",
     "GXY",
     "G5Z",
     "G5T",
     "G5" }};
     */

NAMESPACE_END(Grid)
