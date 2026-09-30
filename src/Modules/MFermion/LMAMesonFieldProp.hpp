/*
 * LMAMesonFieldProp.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2026
 *
 * Author: Michael Lynch <michaellynch628@gmail.com>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 */

/*  END LEGAL */
#ifndef HadronsMILC_MFermion_LMAMesonFieldProp_hpp_
#define HadronsMILC_MFermion_LMAMesonFieldProp_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <A2AMatrix.hpp>
#include <EigenPack.hpp>
#include <Modules/MContraction/MesonField.hpp>
#include <Modules/MFermion/SpinTaste.hpp> // TGammaMap + the gammas module
#include <GridMilc/GridMilc.h>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *        Low mode propagators from precomputed meson fields                *
 ******************************************************************************/
/*  GaugeProp-style replacement of the former file-driven LMA solver
    family: instead of registering lazy Solver closures consumed
    one at a time through GaugeProp middlemen, this module EAGERLY produces
    one scalar PropagatorField per (timeslice, gamma) output name during
    its own execute(); downstream modules fetch the field objects directly
    by name (no middleman, no solver indirection).

    The file's eigenvector rows come in |e+o>/|e-o> pairs (rows 2k/2k+1 of
    eigenpair k) -- produced identically by EigenPackFullPairs (full-volume)
    or by MesonField's checkerboarded cbPairs path. Their sum and
    difference reconstruct the parity-split inner products LowModeProj
    computes live:

      SUM_k(t,j)  ~ <e_cb^k|eta_j,cb>(t)
      DIFF_k(t,j) ~ (i/lam_k) <Meooe(e_cb^k)|eta_j,!cb>(t)

    where cb is the eigenvectors' checkerboard. The rows are |E+O>/|E-O>
    in physical parity, so M[t][2k][j] + M[t][2k+1][j] is the even-site
    partial and M[t][2k][j] - M[t][2k+1][j] the odd-site one: SUM is the
    former for even packs and the latter for odd packs (evenEigen=false),
    and DIFF is the other one.

    giving (see the 2026-09-17 design artifact, decision D8, for the full
    derivation):

      ferm_c = (norm_cb / pairScale) *
               [ SUM_k e_k  -  i * Meooe( SUM_k (DIFF_k/lam_k) e_k ) ]

    with ferm_c placed in COLOR SLOT c of noise n's output propagator:
    each color-diluted noise source occupies 3 adjacent table columns
    (one per color index, TimeDilutedSpinColorDiagonal color-fast
    layout), global noise g = noiseIndex + n in window [3g, 3g + 3);
    the module reconstructs one FermionField per column j = 3g + c and
    assembles them FermToProp-style -- the GaugeProp solveField
    color-loop pattern applied to meson-field columns, run once per
    noise window (per-noise eigenpass loops: each noise's propagator is
    isolated; averaging over noises is deferred to the
    MContraction::Meson vector contract, whose srcCorrs hold the
    per-source correlators and corr their average).

    action      Staggered action module (Meooe parity move)
    lowModes    MassShiftEigenPack module with the CHECKERBOARDED
                eigenvectors/eigenvalues (row pair k <-> evec[k], eval[k])
    gammas      name of an MFermion::SpinTaste module; this module consumes
                the SpinTaste module's `_map` companion output
                (par().gammas + "_map", a TGammaMap<StagGamma>) for
                per-label StagGamma lookup -- NOT the bare vector. The
                module applies no spin-taste operator: each requested
                label's physics pair (_spin/_taste = P, fold-invariant
                under the gammas module's applyG5) is the FILE KEY that
                must match each loaded file's metadata (the producing
                StagA2AMesonField writes gamma_spin/gamma_taste as its own
                P pair). The module's gauge binding is irrelevant here (U
                is only read by appliers)
    labels      REQUIRED whitespace-separated list, positionally parallel
                ONLY to mesonField (labels[i] <-> mesonField[i]); each
                label is looked up by key in the gammas module's `_map` at
                setup (fatal if missing, with a dynamic available-keys
                diagnostic) -- fully replaces the former
                gammas[i]<->mesonField[i] positional contract. Naming
                follows the label (keys the per-timeslice TGammaMap
                outputs), physics follows the resolved StagGamma object. A
                subset of the gammas module's full published list is now
                legal (labels need not cover every entry the gammas
                module publishes)
    mesonField  whitespace-separated LoadMesonField module names, one per
                label (positional parallel list, labels[i] <-> entry i);
                each loader's published MesonFieldMILCMetadata side object
                ("<loader>_metadata") is cross-checked against the
                physics pair (_spin/_taste) of labels[i]'s resolved
                StagGamma object at execute time (mesonField[i] must load
                the file produced under that gamma's physics pair), so a
                miswired label/file pairing fails loudly instead of
                silently mislabeling output names
    noiseIndex  index of the FIRST noise window (unit = NOISES, not
                columns): the module reconstructs noises
                noiseIndex..noiseIndex+nNoise-1, window n reading
                columns (noiseIndex+n)*3 .. (noiseIndex+n)*3+2 (one
                column per color). The table must COVER the windows:
                cols >= (noiseIndex+nNoise)*3 (execute-time check)
    nNoise      number of noise sources reconstructed per output name
                (>= 1). nNoise == 1 publishes a scalar PropagatorField
                per output name (legacy grammar); nNoise > 1 publishes
                a std::vector<PropagatorField> of length nNoise,
                noise-major (entry n = noise n's propagator) -- the
                RandomWall scalar-collapse precedent (ef6d6d0)
    tA          first timeslice to produce (inclusive)
    tB          last timeslice to produce (inclusive, must be < nt)
    tStep       timeslice stride (>= 1); ONE TGammaMap output per
                t in [tA, tB] with stride tStep, named "<name>_t<t>"
                (gamma-free: the label is the map key, resolved
                by consumers such as a StagGaugeProp's guess/source
                lookup at execute time)
    eigStart    first eigenpair to include (pair space)
    nEigs       number of eigenpairs (< 1: all)
    projector   ""/"false" (default): reconstruct LowModeProj's
                project=false (invmag-weighted, deflated-solution)
                formula -- the ONLY form safe to use as a solver guess.
                "true": reconstruct the bare projection (LowModeProj's
                project=true formula, LowModeProj.hpp:206-210) -- same
                units as the source, NOT a solution; do not feed to a
                solver as guess=. std::string (not bool), mirroring
                negFirst: Grid's XmlReader::readDefault aborts the
                process if a bool-typed serializable member's XML node
                is entirely absent (BaseIO.h:518-535's fromString hits
                failbit on an empty stream); a std::string member
                degrades a missing node to "" with only a warning, so
                XMLs that never set this param stay safe
    negFirst    ""/"false" (default): row 2k is |e+o>; "true": |e-o> comes
                first (flips the DIFF sign)
    pairScale   production normalization constant P (default sqrt(2); use
                1.0 for unit-norm full-volume Lanczos files)
    noise       optional name of the noise-vector environment object
                (std::vector<FermionField>) used for a one-shot pairing/
                normalization self-check at execute time, run once per
                noise window n against noise[(noiseIndex + n)*3]
                (full-time-extent sum; only the (G1,G1) table is
                checkable against the live, gamma-independent reference;
                other gammas log a skip)
    a2a_batch   ""/"false" (default): per-timeslice outputs, unchanged
                behavior. "true": publish ONE A2A batch instead -- a
                single-key TGammaMap<std::vector<FermionField>> named
                "<name>" (no per-timeslice "_t<t>" outputs), keyed by the
                single required label, holding nNoise*nSlices*3 columns
                in the shared column order 3*(n*nSlices + j) + c
                (RandomWall order, docs/CONTEXT.md): column (n, j, c) is
                the reconstruction for timeslice tA + j*tStep, noise
                noiseIndex + n, color c -- the same values the
                per-timeslice mode writes into that propagator's color
                slot (same reconstruction path, only the destination
                differs; bit-identical). tA plays RandomWall's t0: the
                USER is responsible for matching tA/tStep/nNoise to the
                RandomWall's t0/tStep/nSrc; when the batch feeds a
                StagGaugeProp guess= the key must equal the gammas
                module's effective label and the guess length is never
                checked (cross-module contracts stay user
                responsibility, like tA/tStep/nNoise matching). Requires
                exactly one label and projector off (setup errors
                otherwise)

    setup() only parses parameters, checks container sizes and allocates
    and zeroes the outputs (the scheduler dry-runs it for memory profiling
    while the loader tables are still empty 0x0 matrices); ALL table-content
    checks (row layout, color-window bounds, metadata cross-check,
    self-check) live at the top of execute(). The _subtract variants of
    the former solver family are dropped (they were never end-to-end
    validated and the producer contract has no caller-supplied source).
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MFermion)

class LMAMesonFieldPropMILCPar : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(LMAMesonFieldPropMILCPar,
                                  std::string,   action,
                                  std::string,   lowModes,
                                  std::string,   mesonField,
                                  std::string,   gammas,
                                  std::string,   labels,
                                  unsigned int,  noiseIndex,
                                  unsigned int,  nNoise,
                                  unsigned int,  tA,
                                  unsigned int,  tB,
                                  unsigned int,  tStep,
                                  unsigned int,  eigStart,
                                  int,           nEigs,
                                  std::string,   projector,
                                  std::string,   negFirst,
                                  std::string,   pairScale,
                                  std::string,   noise,
                                  std::string,   a2a_batch);
  LMAMesonFieldPropMILCPar(void)
      : tStep(1), nNoise(1), projector(""), negFirst(""), pairScale(""),
        a2a_batch("") {}
};
// gammas: name of an MFermion::SpinTaste module; this module consumes its
//         `_map` companion output (par().gammas + "_map") for per-label
//         StagGamma lookup -- NOT the bare vector.
// labels: REQUIRED whitespace-separated list, positionally parallel ONLY
//         to mesonField (labels[i] <-> mesonField[i]); each label is
//         looked up by key in the gammas module's `_map` at setup (fatal
//         if missing) -- fully replaces the former
//         gammas[i]<->mesonField[i] positional contract.

template <typename FImpl, typename Pack>
class TLMAMesonFieldPropMILC : public Module<LMAMesonFieldPropMILCPar> {
public:
  FERM_TYPE_ALIASES(FImpl, );

private:
  // one LoadMesonField module name per label (positional parallel list,
  // strToVec<std::string>); count mismatch vs the labels list is fatal
  // before any positional access (checked at setup, where the objects
  // exist; getInput only wires edges)
  std::vector<std::string> mesonFieldList(
      const std::vector<std::string> &labels) const;
  // REQUIRED labels parameter: one entry per mesonField, positionally
  // parallel (labels[i] <-> mesonField[i]); each looked up by key in the
  // gammas module's `_map` at setup (fatal if missing).
  std::vector<std::string> parseLabels(void) const;
  // the timeslices this instance materializes: t in [tA, tB] stride
  // tStep, clamped to the lattice time extent (single source shared by
  // getOutput()/setup()/execute() -- the name family cannot diverge
  // between the three sites)
  std::vector<unsigned int> sliceTimes(void) const;
  // the per-timeslice TGammaMap output name: "<name>_t<t>" (gamma-free;
  // the label is the map key)
  std::string mapName(const unsigned int t) const;

public:
  // constructor
  TLMAMesonFieldPropMILC(const std::string name);
  // destructor
  virtual ~TLMAMesonFieldPropMILC(void){};
  // dependency relation
  virtual std::vector<std::string> getInput(void);
  virtual std::vector<std::string> getOutput(void);
  // setup
  virtual void setup(void);
  // execution
  virtual void execute(void);
};

MODULE_REGISTER_TMP(StagLMAMesonFieldProp,
                    ARG(TLMAMesonFieldPropMILC<STAGIMPL,
                                               MassShiftEigenPack<STAGIMPL>>),
                    MFermion);

/******************************************************************************
 *                TLMAMesonFieldPropMILC implementation                      *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
TLMAMesonFieldPropMILC<FImpl, Pack>::TLMAMesonFieldPropMILC(
    const std::string name)
    : Module<LMAMesonFieldPropMILCPar>(name) {}

// labels / meson-field parallel lists //////////////////////////////////////////
template <typename FImpl, typename Pack>
std::vector<std::string> TLMAMesonFieldPropMILC<FImpl, Pack>::parseLabels(
    void) const {
  if (par().labels.empty()) {
    HADRONS_ERROR(Argument,
                  "LMAMesonFieldProp requires the 'labels' parameter "
                  "(one label per mesonField entry)");
  }
  return strToVec<std::string>(par().labels);
}

template <typename FImpl, typename Pack>
std::vector<std::string> TLMAMesonFieldPropMILC<FImpl, Pack>::mesonFieldList(
    const std::vector<std::string> &labels) const {
  auto mfs = strToVec<std::string>(par().mesonField);
  if (mfs.size() != labels.size()) {
    HADRONS_ERROR(Argument,
                  "'labels' has " + std::to_string(labels.size()) +
                      " entries but 'mesonField' names " +
                      std::to_string(mfs.size()) +
                      " module(s): labels and mesonField must be "
                      "positionally parallel lists of equal length");
  }

  return mfs;
}

// output name family /////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
std::vector<unsigned int> TLMAMesonFieldPropMILC<FImpl, Pack>::sliceTimes(
    void) const {
  std::vector<unsigned int> ts;
  // defensive: a zero tStep would loop forever below; setup() rejects
  // it loudly, but getOutput() may run before setup(), so guard here
  // too and simply enumerate nothing
  if (par().tStep < 1) {
    return ts;
  }
  const unsigned int nt = static_cast<unsigned int>(env().getDim().back());
  // clamp to the time extent: getOutput() runs before setup() validates
  // the range
  const unsigned int tLast = std::min(par().tB, nt - 1);
  for (unsigned int i = 0;; ++i) {
    const unsigned int t = par().tA + i * par().tStep;
    if ((t < par().tA) || (t > tLast)) {
      break;
    }
    ts.push_back(t);
  }

  return ts;
}

template <typename FImpl, typename Pack>
std::string TLMAMesonFieldPropMILC<FImpl, Pack>::mapName(
    const unsigned int t) const {
  return getName() + "_t" + std::to_string(t);
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
std::vector<std::string> TLMAMesonFieldPropMILC<FImpl, Pack>::getInput(void) {
  if (par().gammas.empty()) {
    HADRONS_ERROR(Argument,
                  "LMAMesonFieldProp requires the 'gammas' SpinTaste "
                  "module name");
  }
  auto labels = parseLabels();
  std::vector<std::string> in{par().action, par().lowModes,
                              par().gammas + "_map"};
  for (auto &mf : mesonFieldList(labels)) {
    in.push_back(mf);
    // explicit edge on the loader's metadata side object (GaugeProp
    // guess-object pattern): naming an env object in getInput() keeps
    // it alive until this module executes, so the execute-time
    // envGet(<loader>_metadata) cross-check cannot hit a GC-freed object
    in.push_back(mf + "_metadata");
  }
  if (!par().noise.empty()) {
    in.push_back(par().noise);
  }

  return in;
}

template <typename FImpl, typename Pack>
std::vector<std::string> TLMAMesonFieldPropMILC<FImpl, Pack>::getOutput(void) {
  if (par().a2a_batch == "true") {
    // batch mode: ONE single-key map output (header); no per-timeslice
    // outputs. An invalid a2a_batch value falls through to the
    // per-timeslice list here and is rejected by setup's three-value
    // check immediately after
    return {getName()};
  }
  std::vector<std::string> out;
  for (auto &t : sliceTimes()) {
    out.push_back(mapName(t));
  }

  return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TLMAMesonFieldPropMILC<FImpl, Pack>::setup(void) {
  int Ls = env().getObjectLs(par().action);
  if (Ls > 1) {
    HADRONS_ERROR(Argument, "Ls > 1 not implemented");
  }

  // optional string parameters (missing XML nodes are tolerated for
  // strings only -- useStencil precedent, commit 297b724)
  if ((par().negFirst != "") && (par().negFirst != "false") &&
      (par().negFirst != "true")) {
    HADRONS_ERROR(Argument, "negFirst must be '', 'false' or 'true' (got '" +
                                par().negFirst + "')");
  }
  if ((par().a2a_batch != "") && (par().a2a_batch != "false") &&
      (par().a2a_batch != "true")) {
    HADRONS_ERROR(Argument, "a2a_batch must be '', 'false' or 'true' (got '" +
                                par().a2a_batch + "')");
  }
  if (!par().pairScale.empty()) {
    auto scale = strToVec<RealD>(par().pairScale);
    if (scale.size() != 1) {
      HADRONS_ERROR(Argument, "pairScale must be a single number (got '" +
                                  par().pairScale + "')");
    }
  }

  // timeslice range: inclusive [tA, tB] with stride tStep >= 1 (the
  // SeqGamma tA/tB idiom, plus the first step parameter in this
  // codebase)
  int nt = env().getDim().back();
  if (par().tStep < 1) {
    HADRONS_ERROR(Argument, "tStep must be >= 1 (got " +
                                std::to_string(par().tStep) + ")");
  }
  // noise count: one propagator per noise window; nNoise == 1 keeps the
  // legacy scalar outputs (RandomWall scalar-collapse precedent)
  if (par().nNoise < 1) {
    HADRONS_ERROR(Argument, "nNoise must be >= 1 (got " +
                                std::to_string(par().nNoise) + ")");
  }
  if (par().tA > par().tB) {
    HADRONS_ERROR(Argument, "tA (" + std::to_string(par().tA) +
                                ") must not exceed tB (" +
                                std::to_string(par().tB) + ")");
  }
  if (par().tB >= static_cast<unsigned int>(nt)) {
    HADRONS_ERROR(Argument, "tB (" + std::to_string(par().tB) +
                                ") out of range: the lattice has " +
                                std::to_string(nt) + " timeslices");
  }

  if (par().gammas.empty()) {
    HADRONS_ERROR(Argument,
                  "LMAMesonFieldProp requires the 'gammas' SpinTaste "
                  "module name");
  }
  auto labels = parseLabels();
  auto mfs = mesonFieldList(labels);

  const auto &gammaMap = envGet(TGammaMap<StagGamma>, par().gammas + "_map");
  // Setup-time key-set validation (GammaMapElement/Meson checkKeys
  // precedent): every requested label must exist in the gammas module's
  // _map, with a dynamic available-keys diagnostic on miss. `selected`
  // restricts the per-timeslice TGammaMap output to exactly the
  // requested labels (a SUBSET of the gammas module's full list is now
  // legal), built via the new TGammaMap(const TGammaMap<StagGamma>&,
  // Args...) ctor overload at each envCreate call below.
  TGammaMap<StagGamma> selected;
  for (auto &label : labels) {
    auto it = gammaMap.find(label);
    if (it == gammaMap.end()) {
      std::string available;
      for (auto &q : gammaMap) {
        available += (available.empty() ? "" : ", ") + q.first;
      }
      HADRONS_ERROR(Argument, "gammas module '" + par().gammas +
                                  "' has no entry for label '" + label +
                                  "' (available: " + available + ")");
    }
    selected.emplace(label, it->second);
  }

  const bool a2aBatch = (par().a2a_batch == "true");
  if (a2aBatch) {
    if (labels.size() != 1) {
      HADRONS_ERROR(Argument,
                    "a2a_batch=\"true\" requires exactly one label (got " +
                        std::to_string(labels.size()) + ": '" + par().labels +
                        "') -- the batch output is a single-key TGammaMap");
    }
    if (par().projector == "true") {
      HADRONS_ERROR(Argument,
                    "a2a_batch=\"true\" is incompatible with "
                    "projector=\"true\": the bare projection is not a "
                    "solution and must not be used as a solver guess "
                    "(leave projector empty for the invmag-weighted "
                    "guess form)");
    }
  }

  LOG(Message) << "Setting up meson-field driven low mode propagator '"
               << getName() << "' for action '" << par().action
               << "' using eigenvectors from '" << par().lowModes
               << "' (noise index " << par().noiseIndex << ", nNoise "
               << par().nNoise << ", gammas module '" << par().gammas
               << "', one map per timeslice in [tA=" << par().tA
               << ", tB=" << par().tB << "] with stride " << par().tStep
               << "):" << std::endl;
  for (unsigned int i = 0; i < labels.size(); ++i) {
    const auto &g = selected.at(labels[i]);
    LOG(Message) << "  label '" << labels[i] << "' (file key '"
                 << StagGamma::GetName(g._spin, g._taste)
                 << "') from '" << mfs[i] << "'" << std::endl;
  }

  auto &epack = envGet(Pack, par().lowModes);

  // eigenpair range: the fixed three-clause bounds check of LowModeProj
  // (commit a64d61b)
  unsigned int eigStart = par().eigStart;
  int nEigs = par().nEigs;
  if (nEigs < 1) {
    nEigs = epack.evec.size();
  }
  if (eigStart > nEigs || eigStart > epack.evec.size() ||
      nEigs - eigStart > epack.evec.size() - eigStart) {
    HADRONS_ERROR(Argument,
                  "Requested eigs (parameters eigStart and nEigs) out of "
                  "bounds.");
  }

  // meson-field tables: one per label; container-size check only (the
  // loader fills the nt-length vector in its own setup, so this is
  // dry-run safe). Row/column/metadata CONTENT checks live in execute()
  for (auto &mfName : mfs) {
    auto &mf = envGet(std::vector<A2AMatrix<HADRONS_A2AM_IO_TYPE>>, mfName);
    if (static_cast<int>(mf.size()) != nt) {
      HADRONS_ERROR(Size, "meson field '" + mfName + "' has " +
                              std::to_string(mf.size()) +
                              " timeslices, expected " +
                              std::to_string(nt));
    }
  }

  // temps: the per-color-column checkerboard accumulators of the
  // reconstruction (SUM channels rbTemp0..2, DIFF channels
  // rbTempNeg0..2) and three Meooe targets rbFermNeg0..2 (one per
  // color column: the three applications run before the single
  // assembly pass, which writes the output propagator directly).
  // envTmp, not the former envCache: the eager module consumes them only
  // inside its own execute(). The six accumulators exist so that ONE
  // element-wise pass per eigenvector can update every color column at
  // once (the columns differ only in their table coefficients)
  static_assert(FImpl::Dimension == 3,
                "color-fused reconstruction assumes three color columns");
  envTmp(FermionField, "rbTemp0", 1, envGetRbGrid(FermionField));
  envTmp(FermionField, "rbTemp1", 1, envGetRbGrid(FermionField));
  envTmp(FermionField, "rbTemp2", 1, envGetRbGrid(FermionField));
  envTmp(FermionField, "rbTempNeg0", 1, envGetRbGrid(FermionField));
  envTmp(FermionField, "rbTempNeg1", 1, envGetRbGrid(FermionField));
  envTmp(FermionField, "rbTempNeg2", 1, envGetRbGrid(FermionField));
  envTmp(FermionField, "rbFermNeg0", 1, envGetRbGrid(FermionField));
  envTmp(FermionField, "rbFermNeg1", 1, envGetRbGrid(FermionField));
  envTmp(FermionField, "rbFermNeg2", 1, envGetRbGrid(FermionField));

  // allocation-only output creation: ONE TGammaMap per timeslice --
  // TGammaMap<PropagatorField> when nNoise == 1 (a scalar field per
  // label), TGammaMap<std::vector<PropagatorField>> of length
  // nNoise otherwise (noise-major). Entries are constructed from
  // `selected` INSIDE the envCreate window (profiler-visible) and
  // zeroed; keys are the SELECTED labels. Sizing from par().nNoise --
  // never from table contents (dry-run safety)
  if (a2aBatch) {
    // ONE batch output: single-key TGammaMap<std::vector<FermionField>>
    // of length nNoise*nSlices*3, full-grid FermionField columns in the
    // shared column order 3*(n*nSlices + j) + c (nSlices = sliceTimes()
    // count). Constructed from `selected` (the one validated label)
    // with a pre-built prototype vector as ONE ctor argument -- the
    // nNoise>1 envCreate precedent -- so construction stays inside the
    // envCreate window; zeroed like every other output (sizing from
    // par(), never from table contents: dry-run safety)
    const unsigned int nSlices = sliceTimes().size();
    const unsigned int nCol =
        par().nNoise * nSlices * FImpl::Dimension;
    envCreate(TGammaMap<std::vector<FermionField>>, getName(), 1, selected,
              std::vector<FermionField>(nCol, envGetGrid(FermionField)));
    for (auto &p :
         envGet(TGammaMap<std::vector<FermionField>>, getName())) {
      for (auto &f : p.second) {
        f = Zero();
      }
    }
  } else {
  for (auto &t : sliceTimes()) {
    if (par().nNoise == 1) {
      envCreate(TGammaMap<PropagatorField>, mapName(t), 1, selected,
                envGetGrid(PropagatorField));
      for (auto &p : envGet(TGammaMap<PropagatorField>, mapName(t))) {
        p.second = Zero();
      }
    } else {
      // A pre-built vector is passed as ONE constructor argument
      // (GaugeProp setupHelper precedent): TGammaMap's ctor loops
      // `this->emplace(label, args...)`, and std::map::emplace forwards
      // args FLAT to pair's constructor (no automatic key/value split
      // without std::piecewise_construct) -- pair has no 3-argument
      // overload, so passing size_t and GridBase* separately does not
      // compile. Bundling them here keeps construction inside the
      // envCreate window (each entry still deep-copies from this
      // prototype during TGammaMap's own emplace loop)
      envCreate(TGammaMap<std::vector<PropagatorField>>, mapName(t), 1,
                selected,
                std::vector<PropagatorField>(par().nNoise,
                                             envGetGrid(PropagatorField)));
      for (auto &p :
           envGet(TGammaMap<std::vector<PropagatorField>>, mapName(t))) {
        for (auto &f : p.second) {
          f = Zero();
        }
      }
    }
  }
  }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TLMAMesonFieldPropMILC<FImpl, Pack>::execute(void) {
  // content checks run here, not in setup(): the scheduler's memory-
  // profiling pass dry-runs every module's setup() BEFORE anything
  // executes, when the loader tables are still nt empty 0x0 matrices
  auto labels = parseLabels();
  auto mfs = mesonFieldList(labels);
  const auto &gammaMap = envGet(TGammaMap<StagGamma>, par().gammas + "_map");
  auto &mat = envGet(FMat, par().action);
  auto &epack = envGet(Pack, par().lowModes);
  int nt = env().getDim().back();

  for (unsigned int i = 0; i < labels.size(); ++i) {
    const StagGamma &gamma = gammaMap.at(labels[i]);
    auto &mf =
        envGet(std::vector<A2AMatrix<HADRONS_A2AM_IO_TYPE>>, mfs[i]);
    if (static_cast<unsigned int>(mf[0].rows()) !=
        2 * static_cast<unsigned int>(epack.evec.size())) {
      HADRONS_ERROR(Size, "meson field '" + mfs[i] + "' has " +
                              std::to_string(mf[0].rows()) +
                              " rows, expected 2*" +
                              std::to_string(epack.evec.size()) +
                              " (|e+o>/|e-o> pair layout)");
    }
    // the table must COVER the requested noise windows: window n reads
    // columns (noiseIndex + n)*Nc .. (noiseIndex + n)*Nc + Nc-1, so the
    // noise index needs (noiseIndex + nNoise)*Nc columns. Superset
    // tables are legal (e.g. an all-timeslices time-diluted table:
    // window 0 is its t=0 noise)
    if (static_cast<unsigned int>(mf[0].cols()) <
        (par().noiseIndex + par().nNoise) * FImpl::Dimension) {
      HADRONS_ERROR(Size,
                    "meson field '" + mfs[i] + "' has " +
                        std::to_string(mf[0].cols()) +
                        " columns, too few for noise windows [noiseIndex=" +
                        std::to_string(par().noiseIndex) +
                        ", noiseIndex+nNoise=" +
                        std::to_string(par().noiseIndex + par().nNoise) +
                        ") (needs " +
                        std::to_string((par().noiseIndex + par().nNoise) *
                                       FImpl::Dimension) +
                        " columns: (noiseIndex+nNoise)*3)");
    }
    // label/file pairing cross-check against the loader's published
    // metadata, on the resolved object's physics pair (_spin/_taste = P;
    // the producing StagA2AMesonField writes gamma_spin/gamma_taste as
    // its own P pair)
    auto &md = envGet(MContraction::MesonFieldMILCMetadata,
                      mfs[i] + "_metadata");
    if ((md.gamma_spin != gamma._spin) || (md.gamma_taste != gamma._taste)) {
      // distinguish a never-filled side object (non-HDF5 build or
      // legacy file) from a genuine label/file miswire
      const std::string fileSt =
          StagGamma::GetName(md.gamma_spin, md.gamma_taste);
      if (fileSt.find("undef") != std::string::npos) {
        HADRONS_ERROR(Argument,
                      "meson field '" + mfs[i] + "' carries undefined "
                      "spin-taste metadata (non-HDF5 build or legacy "
                      "file): cannot cross-check label '" +
                          labels[i] + "'");
      }
      HADRONS_ERROR(Argument,
                    "label '" + labels[i] + "' (file key '" +
                        StagGamma::GetName(gamma._spin, gamma._taste) +
                        "' from gammas module '" + par().gammas +
                        "') is configured for meson field '" + mfs[i] +
                        "', but the file holds spin-taste '" + fileSt +
                        "'");
    }
  }

  // optional pairing/normalization self-check (unchanged from the former
  // solver family): the file pair sum of the first eigenpair, summed
  // over ALL timeslices, equals P * <e|eta_j,E> with the live
  // checkerboarded eigenvector; their ratio exposes the production
  // constants and catches pairing/order/normalization mistakes. The
  // live reference is gamma-independent, so the checkable table is the
  // one whose file CONTENT is the identity pairing: scan the physics
  // pairs (_spin/_taste; with applyG5=true the identity table is
  // reached through its conjugated file key); other labels get a skip
  // notice
  if (!par().noise.empty()) {
    auto &noise = envGet(std::vector<FermionField>, par().noise);
    if ((par().noiseIndex + par().nNoise) * FImpl::Dimension >
        noise.size()) {
      HADRONS_ERROR(Size, "noise windows [noiseIndex, noiseIndex+nNoise) "
                          "out of range for noise object '" +
                              par().noise + "'");
    }
    int gIdentity = -1;
    for (unsigned int i = 0; i < labels.size(); ++i) {
      const StagGamma &gamma = gammaMap.at(labels[i]);
      if ((gamma._spin == StagGamma::StagAlgebra::G1) &&
          (gamma._taste == StagGamma::StagAlgebra::G1)) {
        gIdentity = i;
        break;
      }
    }
    if (gIdentity < 0) {
      LOG(Message) << "Self-check skipped: no (G1,G1) gamma VALUE in "
                      "the list (the live reference <e|eta> is "
                      "gamma-independent)"
                   << std::endl;
    } else {
      auto &mf = envGet(std::vector<A2AMatrix<HADRONS_A2AM_IO_TYPE>>,
                        mfs[gIdentity]);
      unsigned int eigStart = par().eigStart;
      bool negFirst = (par().negFirst == "true");
      RealD pairScale = std::sqrt(2.0);
      if (!par().pairScale.empty()) {
        pairScale = strToVec<RealD>(par().pairScale)[0];
      }

      FermionField rbNoise(envGetRbGrid(FermionField));
      FermionField rbNoiseOdd(envGetRbGrid(FermionField));
      FermionField MeooeE(envGetRbGrid(FermionField));
      int cb = epack.evec[0].Checkerboard();
      int cbNeg = (cb == Even) ? Odd : Even;
      const RealD lam_D = epack.eval[eigStart].imag();

      // Meooe(e_eigStart) doesn't depend on the noise window -- hoisted
      // out of the loop below. Checkerboard labeled BEFORE the call,
      // matching EigenPackCBPairs.hpp:110-117's documented convention
      // and this file's own :1200-1202 precedent
      MeooeE.Checkerboard() = cbNeg;
      mat.Meooe(epack.evec[eigStart], MeooeE);

      // one self-check per noise window (D5): window n checks
      // noise[(noiseIndex + n)*3] against the same table column pair,
      // keeping the object-index == column convention (window-unit
      // noiseIndex: global noise g = noiseIndex + n)
      for (unsigned int n = 0; n < par().nNoise; ++n) {
        const unsigned int j =
            (par().noiseIndex + n) * FImpl::Dimension;

        rbNoise = Zero();
        rbNoise.Checkerboard() = cb;
        pickCheckerboard(cb, rbNoise, noise[j]);

        rbNoiseOdd = Zero();
        rbNoiseOdd.Checkerboard() = cbNeg;
        pickCheckerboard(cbNeg, rbNoiseOdd, noise[j]);

        ComplexD ipFull =
            TensorRemove(innerProduct(epack.evec[eigStart], rbNoise));
        ComplexD ipNegFull = TensorRemove(innerProduct(MeooeE, rbNoiseOdd));
        ComplexD sumFile = 0.;
        ComplexD diffFile = 0.;
        for (int t = 0; t < nt; ++t) {
          // physical-parity partials; swapped onto the eigenvector's
          // checkerboard exactly as in the reconstruction's coeffs
          ComplexD evenPart = ComplexD(mf[t](2 * eigStart, j)) +
                              ComplexD(mf[t](2 * eigStart + 1, j));
          ComplexD oddPart = ComplexD(mf[t](2 * eigStart, j)) -
                             ComplexD(mf[t](2 * eigStart + 1, j));
          if (negFirst) {
            oddPart = -oddPart;
          }
          sumFile += (cb == Even) ? evenPart : oddPart;
          diffFile += (cb == Even) ? oddPart : evenPart;
        }
        // |ip| through .real()/.imag(): std::abs has no overload for the
        // ComplexD (thrust::complex) of GPU builds
        if (std::hypot(ipFull.real(), ipFull.imag()) > 1.e-12) {
          ComplexD pLive = sumFile / ipFull;
          LOG(Message) << "Self-check (label '"
                       << labels[gIdentity]
                       << "', file key '"
                       << StagGamma::GetName(gammaMap.at(labels[gIdentity])._spin,
                                             gammaMap.at(labels[gIdentity])._taste)
                       << "', noise window " << n << "): file-derived "
                       << "production constant P = " << pLive
                       << " (configured pairScale = " << pairScale
                       << ")" << std::endl;
          // modulus via .real()/.imag(): std::abs has no ComplexD
          // overload on GPU builds (see EigenPackCheck)
          const ComplexD dP = pLive - static_cast<RealD>(pairScale);
          if (std::hypot(dP.real(), dP.imag()) >
              0.05 * std::abs(pairScale)) {
            LOG(Warning) << "Meson-field pair normalization mismatch "
                            "(noise window " << n << "): derived P = "
                         << pLive << " but pairScale = " << pairScale
                         << " -- check the |e+o>/|e-o> pair ordering and "
                            "the production normalization"
                         << std::endl;
          }
        } else {
          LOG(Warning) << "Self-check skipped (noise window " << n
                       << "): live inner product <e_eigStart|eta_" << j
                       << "> vanishes" << std::endl;
        }

        // DIFF-channel self-check: validates the file's DIFF row
        // against the module's own documented DIFF_k ~ (i/lam_k)
        // <Meooe(e_E^k)|eta_j,O> contract (module header, :54) --
        // closes the blind spot where only the SUM channel above was
        // ever exercised, even though this is exactly the row the
        // weighted-branch phase-drop defect lived in
        if ((std::abs(lam_D) > 1.e-12) &&
            (std::hypot(ipNegFull.real(), ipNegFull.imag()) > 1.e-12)) {
          const ComplexD iOverLam(0., 1. / lam_D);
          ComplexD pLiveDiff = diffFile / (iOverLam * ipNegFull);
          LOG(Message) << "Self-check DIFF channel (label '"
                       << labels[gIdentity]
                       << "', noise window " << n << "): file-derived "
                       << "production constant P = " << pLiveDiff
                       << " (configured pairScale = " << pairScale
                       << ")" << std::endl;
          const ComplexD dPDiff = pLiveDiff - static_cast<RealD>(pairScale);
          if (std::hypot(dPDiff.real(), dPDiff.imag()) >
              0.05 * std::abs(pairScale)) {
            LOG(Warning) << "Meson-field DIFF-channel normalization/phase "
                            "mismatch (noise window " << n << "): derived "
                            "P = " << pLiveDiff << " but pairScale = "
                         << pairScale << " -- check the Meooe(e_E^k)/"
                            "eta_O phase relation" << std::endl;
          }
        } else {
          LOG(Warning) << "DIFF-channel self-check skipped (noise window "
                       << n << "): lam_D vanishes or live "
                       << "<Meooe(e_eigStart)|eta_" << j << "_odd> vanishes"
                       << std::endl;
        }
      }
    }
  }

  // reconstruction: one output per (t, label) name -- a scalar
  // PropagatorField when nNoise == 1, nNoise propagators otherwise.
  // Color slot c of noise n is reconstructed from the color-diluted
  // source's table column j = (noiseIndex + n)*3 + c (adjacent columns
  // = adjacent colors, TimeDilutedSpinColorDiagonal color-fast layout)
  // and assembled FermToProp-style -- the GaugeProp solveField
  // color-loop pattern applied to meson-field columns, once per noise
  // window. This entire kernel (accelerator_for batching, Meooe
  // application, site-matrix assembly) is UNCHANGED from before this
  // slice -- only the outer per-(label,file) driver loop below changed
  // from vector-index `g` to `labels[i]`/`gammaMap.at(labels[i])`.
  unsigned int eigStart = par().eigStart;
  int nEigs = par().nEigs;
  if (nEigs < 1) {
    nEigs = epack.evec.size();
  }
  bool negFirst = (par().negFirst == "true");
  bool project = (par().projector == "true");
  bool a2aBatch = (par().a2a_batch == "true");
  RealD pairScale = std::sqrt(2.0);
  if (!par().pairScale.empty()) {
    pairScale = strToVec<RealD>(par().pairScale)[0];
  }

  envGetTmp(FermionField, rbTemp0);
  envGetTmp(FermionField, rbTemp1);
  envGetTmp(FermionField, rbTemp2);
  envGetTmp(FermionField, rbTempNeg0);
  envGetTmp(FermionField, rbTempNeg1);
  envGetTmp(FermionField, rbTempNeg2);
  envGetTmp(FermionField, rbFermNeg0);
  envGetTmp(FermionField, rbFermNeg1);
  envGetTmp(FermionField, rbFermNeg2);

  // SUM/DIFF accumulator columns (entry c reconstructs table column
  // nBase + c of the current noise window); array sugar over the named
  // env temps for the per-column post-processing loop below
  FermionField *rbTempC[FImpl::Dimension] = {&rbTemp0, &rbTemp1, &rbTemp2};
  FermionField *rbTempNegC[FImpl::Dimension] = {&rbTempNeg0, &rbTempNeg1,
                                                &rbTempNeg2};
  // one Meooe target per color column (the three applications must
  // coexist: the single assembly kernel below reads all three)
  FermionField *rbFermNegC[FImpl::Dimension] = {&rbFermNeg0, &rbFermNeg1,
                                                &rbFermNeg2};

  int cb = epack.evec[0].Checkerboard();
  RealD norm = 1. / ::sqrt(norm2(epack.evec[0]));

  auto ts = sliceTimes();
  for (unsigned int i = 0; i < labels.size(); ++i) {
    auto &mf =
        envGet(std::vector<A2AMatrix<HADRONS_A2AM_IO_TYPE>>, mfs[i]);
    const std::string &label = labels[i];
    // batch destination: fetched ONCE per label (single-key map, one
    // entry); the per-slice destinations below do not exist in batch
    // mode and must not be envGot
    std::vector<FermionField> *batchVec = nullptr;
    if (a2aBatch) {
      batchVec =
          &(envGet(TGammaMap<std::vector<FermionField>>, getName()))
               .at(label);
    }
    for (unsigned int j = 0; j < ts.size(); ++j) {
      const unsigned int t = ts[j];
      const A2AMatrix<HADRONS_A2AM_IO_TYPE> &mft = mf[t];
      std::vector<PropagatorField> *propVec = nullptr;
      PropagatorField *propScalar = nullptr;
      if (!a2aBatch) {
        if (par().nNoise == 1) {
          // parenthesized: envGet(...).at(label) unparenthesized binds
          // .at() to the raw getObject<T>() pointer BEFORE envGet's
          // leading dereference applies (macro-expansion precedence)
          propScalar =
              &(envGet(TGammaMap<PropagatorField>, mapName(t))).at(label);
        } else {
          propVec = &(envGet(TGammaMap<std::vector<PropagatorField>>,
                             mapName(t)))
                         .at(label);
        }
      }
      // per-noise reconstruction (per-noise eigenpass loops): noise
      // window n is global noise noiseIndex + n, occupying columns
      // (noiseIndex+n)*3 .. (noiseIndex+n)*3 + 2; its three color
      // columns reconstruct that noise's propagator in isolation -- the
      // accumulators, kernel passes and Meooe post-processing below
      // run once per noise
      for (unsigned int n = 0; n < par().nNoise; ++n) {
        // destination handle: the per-slice propagator (default modes)
        // or none (batch mode -- the assembly seam below addresses the
        // batch columns directly). Everything from the eigenpass to the
        // Meooe applications is destination-agnostic
        PropagatorField *propP = nullptr;
        if (!a2aBatch) {
          propP = (par().nNoise == 1) ? propScalar : &(*propVec)[n];
        }
        const unsigned int nBase =
            (par().noiseIndex + n) * FImpl::Dimension;
        // no zero-init of prop: the single assembly pass below writes
        // EVERY element of every site matrix (all rows, all three
        // color columns) from freshly computed values

        // set the six per-column accumulator checkerboard flags (SUM
        // channels rbTempC, DIFF channels rbTempNegC). The accumulators
        // are NOT zero-initialized when at least one eigenpair follows:
        // the eigenpass opens with a WRITE of the first eigenpair's
        // contribution (firstUnit below), deleting six zero kernels per
        // output. Arithmetic matches the former add-to-zero sequence
        // except that x + 0.0 becomes x, which differs only in the SIGN
        // of exact-zero components (IEEE: -0.0 + 0.0 = +0.0); every
        // nonzero value is bit-equal. Degenerate empty eigenpair range
        // (eigStart == nEigs): no kernel runs, so the accumulators keep
        // their explicit zero-init -- the tail would otherwise read
        // garbage
        const bool emptyEigenpass = (static_cast<int>(eigStart) >= nEigs);
        for (unsigned int c = 0; c < FImpl::Dimension; ++c) {
          if (emptyEigenpass) {
            *rbTempC[c] = Zero();
            *rbTempNegC[c] = Zero();
          }
          rbTempC[c]->Checkerboard() = cb;
          rbTempNegC[c]->Checkerboard() = cb;
        }

        // accumulate the two parity channels from the file row pairs; the
        // subtraction channel is accumulated as DIFF/lam so that the live
        // Meooe below completes the 1/lam_D^2 weighting of LowModeProj.
        // ALL color columns are built from the SAME eigenvectors, so one
        // fused element-wise pass per eigenpair updates all six
        // accumulators: the per-element expressions replicate the former
        // per-column axpy calls verbatim (bit-identical arithmetic) while
        // streaming e once for all three columns -- a third of the kernel
        // invocations and eigenvector traffic of the per-column version.
        // The hoisted aliased Read+Write view pairs mirror Grid's own
        // in-place axpy idiom; the extra scope closes every view before
        // the Meooe/setCheckerboard sequence below: those open CPU-mode
        // views, and the Grid memory manager asserts against mixing
        // concurrent view families on one buffer
        {
          autoView(t0W_v, rbTemp0, AcceleratorWrite);
          autoView(t0R_v, rbTemp0, AcceleratorRead);
          autoView(t1W_v, rbTemp1, AcceleratorWrite);
          autoView(t1R_v, rbTemp1, AcceleratorRead);
          autoView(t2W_v, rbTemp2, AcceleratorWrite);
          autoView(t2R_v, rbTemp2, AcceleratorRead);
          autoView(n0W_v, rbTempNeg0, AcceleratorWrite);
          autoView(n0R_v, rbTempNeg0, AcceleratorRead);
          autoView(n1W_v, rbTempNeg1, AcceleratorWrite);
          autoView(n1R_v, rbTempNeg1, AcceleratorRead);
          autoView(n2W_v, rbTempNeg2, AcceleratorWrite);
          autoView(n2R_v, rbTempNeg2, AcceleratorRead);
          // per-eigenpair table coefficients (expressions identical to the
          // former inline computation -- value- and order-identical)
          auto coeffs = [&](const int k, ComplexD sumC[], ComplexD negC[]) {
            const RealD mass = epack.eval[k].real();
            const RealD lam_D = epack.eval[k].imag();
            const RealD invmag = 1. / (mass * mass + lam_D * lam_D);
            // the file DIFF row carries an inherent i/lam_k phase
            // (module header, DIFF_k contract); the bare branch's
            // negC=diff/lam_D relies on the assembly's single downstream
            // negI (:1230) to supply it, but the weighted branch mixes
            // SUM and DIFF content together, so the cross-channel term
            // in each accumulator needs the phase applied explicitly --
            // negI reused here matches the assembly's constant exactly
            const ComplexD negI(0., -1.);
            for (unsigned int c = 0; c < FImpl::Dimension; ++c) {
              const unsigned int j = nBase + c;
              // rows 2k/2k+1 are |E+O>/|E-O> in PHYSICAL parity
              // (MesonField.hpp reconstructLegacy), so their sum is the
              // even-site partial and their difference the odd-site one.
              // SUM must be the eigenvector's own checkerboard: swap the
              // channels for odd-checkerboard packs (evenEigen=false)
              ComplexD evenPart =
                  ComplexD(mft(2 * k, j)) + ComplexD(mft(2 * k + 1, j));
              ComplexD oddPart =
                  ComplexD(mft(2 * k, j)) - ComplexD(mft(2 * k + 1, j));
              if (negFirst) {
                oddPart = -oddPart;
              }
              const ComplexD sum = (cb == Even) ? evenPart : oddPart;
              const ComplexD diff = (cb == Even) ? oddPart : evenPart;
              if (project) {
                sumC[c] = sum;
                negC[c] = diff / lam_D;
              } else {
                sumC[c] = invmag * (mass * sum + negI * lam_D * diff);
                negC[c] = invmag * (mass * diff / lam_D + negI * sum);
              }
            }
          };

          // batches of four eigenpairs per kernel launch: thread dispatch
          // dominates these small kernels, and batching leaves the
          // per-element accumulation order exactly as before (k strictly
          // descending; batch statements in k order), so results stay
          // bit-identical. The FIRST kernel (batch or single) WRITES the
          // highest-k contribution instead of accumulating onto zero
          // (the accumulators are deliberately not zero-initialized):
          // value-identical except for the sign of exact-zero
          // components. The trailing nEigs % 4 eigenpairs fall through
          // to the single-eigenpair kernel below
          bool firstUnit = true;
          int kHi = static_cast<int>(eigStart) + nEigs - 1;
          // batches of EIGHT eigenpairs per kernel launch: doubles the
          // batch-of-four scheme below. Every kernel re-reads and
          // re-writes all six accumulators, so accumulator traffic per
          // output scales inversely with batch size -- the covered
          // eigenpairs pay half the accumulator traffic of pairs of
          // four-batches, and thread dispatch amortizes twice as far.
          // Statement order per accumulator stays strictly
          // k-descending, so arithmetic is bit-identical to running
          // the same eigenpairs through batch-of-four kernels. A
          // remaining nEigs % 8 in [4,8) still batches four in the
          // loop below, then the singles finish
          for (; kHi - 7 >= int(eigStart); kHi -= 8) {
            const int kA = kHi,     kB = kHi - 1, kC = kHi - 2, kD = kHi - 3,
                      kE = kHi - 4, kF = kHi - 5, kG = kHi - 6, kH = kHi - 7;
            const FermionField &eA = epack.evec[kA];
            const FermionField &eB = epack.evec[kB];
            const FermionField &eC = epack.evec[kC];
            const FermionField &eD = epack.evec[kD];
            const FermionField &eE = epack.evec[kE];
            const FermionField &eF = epack.evec[kF];
            const FermionField &eG = epack.evec[kG];
            const FermionField &eH = epack.evec[kH];

            ComplexD sumC[FImpl::Dimension], negC[FImpl::Dimension];
            coeffs(kA, sumC, negC);
            const ComplexD sA0 = sumC[0], sA1 = sumC[1], sA2 = sumC[2];
            const ComplexD dA0 = negC[0], dA1 = negC[1], dA2 = negC[2];
            coeffs(kB, sumC, negC);
            const ComplexD sB0 = sumC[0], sB1 = sumC[1], sB2 = sumC[2];
            const ComplexD dB0 = negC[0], dB1 = negC[1], dB2 = negC[2];
            coeffs(kC, sumC, negC);
            const ComplexD sC0 = sumC[0], sC1 = sumC[1], sC2 = sumC[2];
            const ComplexD dC0 = negC[0], dC1 = negC[1], dC2 = negC[2];
            coeffs(kD, sumC, negC);
            const ComplexD sD0 = sumC[0], sD1 = sumC[1], sD2 = sumC[2];
            const ComplexD dD0 = negC[0], dD1 = negC[1], dD2 = negC[2];
            coeffs(kE, sumC, negC);
            const ComplexD sE0 = sumC[0], sE1 = sumC[1], sE2 = sumC[2];
            const ComplexD dE0 = negC[0], dE1 = negC[1], dE2 = negC[2];
            coeffs(kF, sumC, negC);
            const ComplexD sF0 = sumC[0], sF1 = sumC[1], sF2 = sumC[2];
            const ComplexD dF0 = negC[0], dF1 = negC[1], dF2 = negC[2];
            coeffs(kG, sumC, negC);
            const ComplexD sG0 = sumC[0], sG1 = sumC[1], sG2 = sumC[2];
            const ComplexD dG0 = negC[0], dG1 = negC[1], dG2 = negC[2];
            coeffs(kH, sumC, negC);
            const ComplexD sH0 = sumC[0], sH1 = sumC[1], sH2 = sumC[2];
            const ComplexD dH0 = negC[0], dH1 = negC[1], dH2 = negC[2];

            const bool first = firstUnit;
            firstUnit = false;
            autoView(eAR_v, eA, AcceleratorRead);
            autoView(eBR_v, eB, AcceleratorRead);
            autoView(eCR_v, eC, AcceleratorRead);
            autoView(eDR_v, eD, AcceleratorRead);
            autoView(eER_v, eE, AcceleratorRead);
            autoView(eFR_v, eF, AcceleratorRead);
            autoView(eGR_v, eG, AcceleratorRead);
            autoView(eHR_v, eH, AcceleratorRead);
            accelerator_for(ss, eAR_v.size(),
                            FermionField::vector_type::Nsimd(), {
              auto evA = coalescedRead(eAR_v[ss]);
              auto evB = coalescedRead(eBR_v[ss]);
              auto evC = coalescedRead(eCR_v[ss]);
              auto evD = coalescedRead(eDR_v[ss]);
              auto evE = coalescedRead(eER_v[ss]);
              auto evF = coalescedRead(eFR_v[ss]);
              auto evG = coalescedRead(eGR_v[ss]);
              auto evH = coalescedRead(eHR_v[ss]);
              // per accumulator the statements run k-descending (the A
              // group is hoisted into the write-or-accumulate branch);
              // accumulators are independent buffers, so the grouping
              // is sequence-preserving. first: the A statements OPEN
              // the accumulation (WRITE, not add-to-zero)
              if (first) {
                coalescedWrite(t0W_v[ss], sA0 * evA);
                coalescedWrite(t1W_v[ss], sA1 * evA);
                coalescedWrite(t2W_v[ss], sA2 * evA);
                coalescedWrite(n0W_v[ss], dA0 * evA);
                coalescedWrite(n1W_v[ss], dA1 * evA);
                coalescedWrite(n2W_v[ss], dA2 * evA);
              } else {
                coalescedWrite(t0W_v[ss],
                               sA0 * evA + coalescedRead(t0R_v[ss]));
                coalescedWrite(t1W_v[ss],
                               sA1 * evA + coalescedRead(t1R_v[ss]));
                coalescedWrite(t2W_v[ss],
                               sA2 * evA + coalescedRead(t2R_v[ss]));
                coalescedWrite(n0W_v[ss],
                               dA0 * evA + coalescedRead(n0R_v[ss]));
                coalescedWrite(n1W_v[ss],
                               dA1 * evA + coalescedRead(n1R_v[ss]));
                coalescedWrite(n2W_v[ss],
                               dA2 * evA + coalescedRead(n2R_v[ss]));
              }
              coalescedWrite(t0W_v[ss],
                             sB0 * evB + coalescedRead(t0R_v[ss]));
              coalescedWrite(t0W_v[ss],
                             sC0 * evC + coalescedRead(t0R_v[ss]));
              coalescedWrite(t0W_v[ss],
                             sD0 * evD + coalescedRead(t0R_v[ss]));
              coalescedWrite(t0W_v[ss],
                             sE0 * evE + coalescedRead(t0R_v[ss]));
              coalescedWrite(t0W_v[ss],
                             sF0 * evF + coalescedRead(t0R_v[ss]));
              coalescedWrite(t0W_v[ss],
                             sG0 * evG + coalescedRead(t0R_v[ss]));
              coalescedWrite(t0W_v[ss],
                             sH0 * evH + coalescedRead(t0R_v[ss]));
              coalescedWrite(t1W_v[ss],
                             sB1 * evB + coalescedRead(t1R_v[ss]));
              coalescedWrite(t1W_v[ss],
                             sC1 * evC + coalescedRead(t1R_v[ss]));
              coalescedWrite(t1W_v[ss],
                             sD1 * evD + coalescedRead(t1R_v[ss]));
              coalescedWrite(t1W_v[ss],
                             sE1 * evE + coalescedRead(t1R_v[ss]));
              coalescedWrite(t1W_v[ss],
                             sF1 * evF + coalescedRead(t1R_v[ss]));
              coalescedWrite(t1W_v[ss],
                             sG1 * evG + coalescedRead(t1R_v[ss]));
              coalescedWrite(t1W_v[ss],
                             sH1 * evH + coalescedRead(t1R_v[ss]));
              coalescedWrite(t2W_v[ss],
                             sB2 * evB + coalescedRead(t2R_v[ss]));
              coalescedWrite(t2W_v[ss],
                             sC2 * evC + coalescedRead(t2R_v[ss]));
              coalescedWrite(t2W_v[ss],
                             sD2 * evD + coalescedRead(t2R_v[ss]));
              coalescedWrite(t2W_v[ss],
                             sE2 * evE + coalescedRead(t2R_v[ss]));
              coalescedWrite(t2W_v[ss],
                             sF2 * evF + coalescedRead(t2R_v[ss]));
              coalescedWrite(t2W_v[ss],
                             sG2 * evG + coalescedRead(t2R_v[ss]));
              coalescedWrite(t2W_v[ss],
                             sH2 * evH + coalescedRead(t2R_v[ss]));
              coalescedWrite(n0W_v[ss],
                             dB0 * evB + coalescedRead(n0R_v[ss]));
              coalescedWrite(n0W_v[ss],
                             dC0 * evC + coalescedRead(n0R_v[ss]));
              coalescedWrite(n0W_v[ss],
                             dD0 * evD + coalescedRead(n0R_v[ss]));
              coalescedWrite(n0W_v[ss],
                             dE0 * evE + coalescedRead(n0R_v[ss]));
              coalescedWrite(n0W_v[ss],
                             dF0 * evF + coalescedRead(n0R_v[ss]));
              coalescedWrite(n0W_v[ss],
                             dG0 * evG + coalescedRead(n0R_v[ss]));
              coalescedWrite(n0W_v[ss],
                             dH0 * evH + coalescedRead(n0R_v[ss]));
              coalescedWrite(n1W_v[ss],
                             dB1 * evB + coalescedRead(n1R_v[ss]));
              coalescedWrite(n1W_v[ss],
                             dC1 * evC + coalescedRead(n1R_v[ss]));
              coalescedWrite(n1W_v[ss],
                             dD1 * evD + coalescedRead(n1R_v[ss]));
              coalescedWrite(n1W_v[ss],
                             dE1 * evE + coalescedRead(n1R_v[ss]));
              coalescedWrite(n1W_v[ss],
                             dF1 * evF + coalescedRead(n1R_v[ss]));
              coalescedWrite(n1W_v[ss],
                             dG1 * evG + coalescedRead(n1R_v[ss]));
              coalescedWrite(n1W_v[ss],
                             dH1 * evH + coalescedRead(n1R_v[ss]));
              coalescedWrite(n2W_v[ss],
                             dB2 * evB + coalescedRead(n2R_v[ss]));
              coalescedWrite(n2W_v[ss],
                             dC2 * evC + coalescedRead(n2R_v[ss]));
              coalescedWrite(n2W_v[ss],
                             dD2 * evD + coalescedRead(n2R_v[ss]));
              coalescedWrite(n2W_v[ss],
                             dE2 * evE + coalescedRead(n2R_v[ss]));
              coalescedWrite(n2W_v[ss],
                             dF2 * evF + coalescedRead(n2R_v[ss]));
              coalescedWrite(n2W_v[ss],
                             dG2 * evG + coalescedRead(n2R_v[ss]));
              coalescedWrite(n2W_v[ss],
                             dH2 * evH + coalescedRead(n2R_v[ss]));
            });
          }
          for (; kHi - 3 >= int(eigStart); kHi -= 4) {
            const int kA = kHi, kB = kHi - 1, kC = kHi - 2, kD = kHi - 3;
            const FermionField &eA = epack.evec[kA];
            const FermionField &eB = epack.evec[kB];
            const FermionField &eC = epack.evec[kC];
            const FermionField &eD = epack.evec[kD];

            ComplexD sumC[FImpl::Dimension], negC[FImpl::Dimension];
            coeffs(kA, sumC, negC);
            const ComplexD sA0 = sumC[0], sA1 = sumC[1], sA2 = sumC[2];
            const ComplexD dA0 = negC[0], dA1 = negC[1], dA2 = negC[2];
            coeffs(kB, sumC, negC);
            const ComplexD sB0 = sumC[0], sB1 = sumC[1], sB2 = sumC[2];
            const ComplexD dB0 = negC[0], dB1 = negC[1], dB2 = negC[2];
            coeffs(kC, sumC, negC);
            const ComplexD sC0 = sumC[0], sC1 = sumC[1], sC2 = sumC[2];
            const ComplexD dC0 = negC[0], dC1 = negC[1], dC2 = negC[2];
            coeffs(kD, sumC, negC);
            const ComplexD sD0 = sumC[0], sD1 = sumC[1], sD2 = sumC[2];
            const ComplexD dD0 = negC[0], dD1 = negC[1], dD2 = negC[2];

            const bool first = firstUnit;
            firstUnit = false;
            autoView(eAR_v, eA, AcceleratorRead);
            autoView(eBR_v, eB, AcceleratorRead);
            autoView(eCR_v, eC, AcceleratorRead);
            autoView(eDR_v, eD, AcceleratorRead);
            accelerator_for(ss, eAR_v.size(),
                            FermionField::vector_type::Nsimd(), {
              auto evA = coalescedRead(eAR_v[ss]);
              auto evB = coalescedRead(eBR_v[ss]);
              auto evC = coalescedRead(eCR_v[ss]);
              auto evD = coalescedRead(eDR_v[ss]);
              // statements in k order (A,B,C,D): each accumulator's
              // update sequence per element matches the unbatched kernel
              // exactly (the R and W views alias one buffer, so later
              // statements read the value written by earlier ones).
              // Reordering the six A statements into one block is
              // sequence-preserving: every accumulator is an independent
              // buffer, and each keeps its own A-then-B/C/D order.
              // first: the A statements OPEN the accumulation (WRITE,
              // not add-to-zero -- see firstUnit)
              if (first) {
                coalescedWrite(t0W_v[ss], sA0 * evA);
                coalescedWrite(t1W_v[ss], sA1 * evA);
                coalescedWrite(t2W_v[ss], sA2 * evA);
                coalescedWrite(n0W_v[ss], dA0 * evA);
                coalescedWrite(n1W_v[ss], dA1 * evA);
                coalescedWrite(n2W_v[ss], dA2 * evA);
              } else {
                coalescedWrite(t0W_v[ss],
                               sA0 * evA + coalescedRead(t0R_v[ss]));
                coalescedWrite(t1W_v[ss],
                               sA1 * evA + coalescedRead(t1R_v[ss]));
                coalescedWrite(t2W_v[ss],
                               sA2 * evA + coalescedRead(t2R_v[ss]));
                coalescedWrite(n0W_v[ss],
                               dA0 * evA + coalescedRead(n0R_v[ss]));
                coalescedWrite(n1W_v[ss],
                               dA1 * evA + coalescedRead(n1R_v[ss]));
                coalescedWrite(n2W_v[ss],
                               dA2 * evA + coalescedRead(n2R_v[ss]));
              }
              coalescedWrite(t0W_v[ss],
                             sB0 * evB + coalescedRead(t0R_v[ss]));
              coalescedWrite(t0W_v[ss],
                             sC0 * evC + coalescedRead(t0R_v[ss]));
              coalescedWrite(t0W_v[ss],
                             sD0 * evD + coalescedRead(t0R_v[ss]));
              coalescedWrite(t1W_v[ss],
                             sB1 * evB + coalescedRead(t1R_v[ss]));
              coalescedWrite(t1W_v[ss],
                             sC1 * evC + coalescedRead(t1R_v[ss]));
              coalescedWrite(t1W_v[ss],
                             sD1 * evD + coalescedRead(t1R_v[ss]));
              coalescedWrite(t2W_v[ss],
                             sB2 * evB + coalescedRead(t2R_v[ss]));
              coalescedWrite(t2W_v[ss],
                             sC2 * evC + coalescedRead(t2R_v[ss]));
              coalescedWrite(t2W_v[ss],
                             sD2 * evD + coalescedRead(t2R_v[ss]));
              coalescedWrite(n0W_v[ss],
                             dB0 * evB + coalescedRead(n0R_v[ss]));
              coalescedWrite(n0W_v[ss],
                             dC0 * evC + coalescedRead(n0R_v[ss]));
              coalescedWrite(n0W_v[ss],
                             dD0 * evD + coalescedRead(n0R_v[ss]));
              coalescedWrite(n1W_v[ss],
                             dB1 * evB + coalescedRead(n1R_v[ss]));
              coalescedWrite(n1W_v[ss],
                             dC1 * evC + coalescedRead(n1R_v[ss]));
              coalescedWrite(n1W_v[ss],
                             dD1 * evD + coalescedRead(n1R_v[ss]));
              coalescedWrite(n2W_v[ss],
                             dB2 * evB + coalescedRead(n2R_v[ss]));
              coalescedWrite(n2W_v[ss],
                             dC2 * evC + coalescedRead(n2R_v[ss]));
              coalescedWrite(n2W_v[ss],
                             dD2 * evD + coalescedRead(n2R_v[ss]));
            });
          }
          // trailing eigenpairs (nEigs % 4): the original single-eigenpair
          // kernel; its first launch WRITES when no batch preceded it
          // (nEigs < 4)
          for (; kHi >= int(eigStart); --kHi) {
            const FermionField &e = epack.evec[kHi];

            ComplexD sumC[FImpl::Dimension], negC[FImpl::Dimension];
            coeffs(kHi, sumC, negC);
            const ComplexD s0 = sumC[0], s1 = sumC[1], s2 = sumC[2];
            const ComplexD d0 = negC[0], d1 = negC[1], d2 = negC[2];

            const bool first = firstUnit;
            firstUnit = false;
            autoView(eR_v, e, AcceleratorRead);
            accelerator_for(ss, eR_v.size(),
                            FermionField::vector_type::Nsimd(), {
              auto ev = coalescedRead(eR_v[ss]);
              if (first) {
                coalescedWrite(t0W_v[ss], s0 * ev);
                coalescedWrite(t1W_v[ss], s1 * ev);
                coalescedWrite(t2W_v[ss], s2 * ev);
                coalescedWrite(n0W_v[ss], d0 * ev);
                coalescedWrite(n1W_v[ss], d1 * ev);
                coalescedWrite(n2W_v[ss], d2 * ev);
              } else {
                coalescedWrite(t0W_v[ss], s0 * ev + coalescedRead(t0R_v[ss]));
                coalescedWrite(t1W_v[ss], s1 * ev + coalescedRead(t1R_v[ss]));
                coalescedWrite(t2W_v[ss], s2 * ev + coalescedRead(t2R_v[ss]));
                coalescedWrite(n0W_v[ss], d0 * ev + coalescedRead(n0R_v[ss]));
                coalescedWrite(n1W_v[ss], d1 * ev + coalescedRead(n1R_v[ss]));
                coalescedWrite(n2W_v[ss], d2 * ev + coalescedRead(n2R_v[ss]));
              }
            });
          }
        }

        // ferm_c = (norm/pairScale) *
        //          [ SUM-channel - i * Meooe(DIFF/lam-channel) ]
        // per color column c, from the c-th accumulators above: the
        // three Meooe applications run FIRST (one per column, into
        // three dedicated buffers), then ONE full-grid assembly
        // kernel writes ALL THREE columns of every site matrix --
        // launches drop to four per output (3 Meooe + 1 assembly) and
        // the propagator is stored once instead of three times
        for (unsigned int c = 0; c < FImpl::Dimension; ++c) {
          // no zero-init: Meooe fully overwrites its output
          // (DhopImproved opens out AcceleratorWrite and
          // coalescedWrites every site, never reading it; proven
          // bit-identical against stale buffer content in a prior
          // experiment). Only the checkerboard flag is metadata-set
          rbFermNegC[c]->Checkerboard() = (cb == Even) ? Odd : Even;

          mat.Meooe(*rbTempNegC[c], *rbFermNegC[c]);
        }

        if (a2aBatch) {
          // batch assembly: three full-grid FermionField columns
          // batch[3*(n*nSlices + j) + c]. Per site, per column, the
          // value is EXACTLY what the site-matrix kernel stores into
          // pmat()()(r, c) -- the same ssh-indexed reads and the same
          // scale/negI scalar expressions; the FermionField site object
          // IS that column's 3-component color vector, stored whole.
          // One fused kernel writes all three columns (three
          // destination views), mirroring the site-matrix kernel's
          // fused structure. Every element is written (both parities):
          // no zero-init
          const unsigned int nSlices = ts.size();
          FermionField &outF0 = (*batchVec)[3 * (n * nSlices + j) + 0];
          FermionField &outF1 = (*batchVec)[3 * (n * nSlices + j) + 1];
          FermionField &outF2 = (*batchVec)[3 * (n * nSlices + j) + 2];
          const GridBase *halfGrid = rbFermNegC[0]->Grid();
          const Coordinate rdimFull = outF0.Grid()->_rdimensions;
          const Coordinate rdimHalf = halfGrid->_rdimensions;
          const Coordinate cbMask = halfGrid->_checker_dim_mask;
          const Coordinate ostride = halfGrid->_ostride;
          const int ndim = halfGrid->_ndimension;
          const RealD scale = norm / pairScale;
          const ComplexD negI(0., -1.);
          const int cbSum = cb;

          // named references first: autoView(n, *ptr[c], m) expands to
          // *ptr[c].View(m) -- '.' binds tighter than '*' (ledger)
          FermionField &sumF0 = *rbTempC[0];
          FermionField &sumF1 = *rbTempC[1];
          FermionField &sumF2 = *rbTempC[2];
          FermionField &negF0 = *rbFermNegC[0];
          FermionField &negF1 = *rbFermNegC[1];
          FermionField &negF2 = *rbFermNegC[2];
          autoView(out0W, outF0, AcceleratorWrite);
          autoView(out1W, outF1, AcceleratorWrite);
          autoView(out2W, outF2, AcceleratorWrite);
          autoView(sumR0, sumF0, AcceleratorRead);
          autoView(sumR1, sumF1, AcceleratorRead);
          autoView(sumR2, sumF2, AcceleratorRead);
          autoView(negR0, negF0, AcceleratorRead);
          autoView(negR1, negF1, AcceleratorRead);
          autoView(negR2, negF2, AcceleratorRead);
          accelerator_for(ss, outF0.Grid()->oSites(),
                          FermionField::vector_type::Nsimd(), {
            Coordinate coor;
            int linear = 0;
            Lexicographic::CoorFromIndex(coor, ss, rdimFull);
            for (int d = 0; d < ndim; ++d) {
              if (cbMask[d]) {
                linear += coor[d];
              }
            }
            int ssh = 0;
            for (int d = 0; d < ndim; ++d) {
              if (d == 0) {
                ssh += ostride[d] * ((coor[d] / 2) % rdimHalf[d]);
              } else {
                ssh += ostride[d] * (coor[d] % rdimHalf[d]);
              }
            }
            if ((linear & 0x1) == cbSum) {
              coalescedWrite(out0W[ss],
                             coalescedRead(sumR0[ssh]) * scale);
              coalescedWrite(out1W[ss],
                             coalescedRead(sumR1[ssh]) * scale);
              coalescedWrite(out2W[ss],
                             coalescedRead(sumR2[ssh]) * scale);
            } else {
              coalescedWrite(out0W[ss],
                             (negI * coalescedRead(negR0[ssh])) * scale);
              coalescedWrite(out1W[ss],
                             (negI * coalescedRead(negR1[ssh])) * scale);
              coalescedWrite(out2W[ss],
                             (negI * coalescedRead(negR2[ssh])) * scale);
            }
          });
        } else {
          PropagatorField &prop = *propP;
        // single assembly pass over the full-grid output propagator.
        // Per site, per column, the op sequence is IDENTICAL to the
        // former per-column kernel -- cb-parity sites: sum * scale;
        // other-parity sites: ((0,-1) * meooe) * scale -- with the
        // same scalar types and operand order the lattice-level ops
        // lowered to (Lattice *= lowers to (*this)*r, Lattice_base.h;
        // the complex scalar product lowers to the same tensor
        // operator* the eigenpass coefficients use). Columns are
        // distinct matrix elements, so interleaving the three
        // columns' statements is sequence-preserving. Every element
        // of every site matrix is written (all rows, all columns):
        // no zero-init and no read-modify-write -- a fresh site
        // object is filled and stored. The full-grid -> rb-grid site
        // mapping replicates Grid's acceleratorSetCheckerboard
        // (Lattice_transfer.h): coordinate from _rdimensions, parity
        // from _checker_dim_mask, rb index from _ostride with the
        // checker dim halved
        {
          const GridBase *halfGrid = rbFermNegC[0]->Grid();
          const Coordinate rdimFull = prop.Grid()->_rdimensions;
          const Coordinate rdimHalf = halfGrid->_rdimensions;
          const Coordinate cbMask = halfGrid->_checker_dim_mask;
          const Coordinate ostride = halfGrid->_ostride;
          const int ndim = halfGrid->_ndimension;
          const RealD scale = norm / pairScale;
          const ComplexD negI(0., -1.);
          const int cbSum = cb;

          autoView(propW, prop, AcceleratorWrite);
          autoView(propR, prop, AcceleratorRead);
          // named references first: autoView(n, *ptr[c], m) expands to
          // *ptr[c].View(m) -- '.' binds tighter than '*' (ledger)
          FermionField &sumF0 = *rbTempC[0];
          FermionField &sumF1 = *rbTempC[1];
          FermionField &sumF2 = *rbTempC[2];
          FermionField &negF0 = *rbFermNegC[0];
          FermionField &negF1 = *rbFermNegC[1];
          FermionField &negF2 = *rbFermNegC[2];
          autoView(sumR0, sumF0, AcceleratorRead);
          autoView(sumR1, sumF1, AcceleratorRead);
          autoView(sumR2, sumF2, AcceleratorRead);
          autoView(negR0, negF0, AcceleratorRead);
          autoView(negR1, negF1, AcceleratorRead);
          autoView(negR2, negF2, AcceleratorRead);
          accelerator_for(ss, prop.Grid()->oSites(),
                          FermionField::vector_type::Nsimd(), {
            Coordinate coor;
            int linear = 0;
            Lexicographic::CoorFromIndex(coor, ss, rdimFull);
            for (int d = 0; d < ndim; ++d) {
              if (cbMask[d]) {
                linear += coor[d];
              }
            }
            int ssh = 0;
            for (int d = 0; d < ndim; ++d) {
              if (d == 0) {
                ssh += ostride[d] * ((coor[d] / 2) % rdimHalf[d]);
              } else {
                ssh += ostride[d] * (coor[d] % rdimHalf[d]);
              }
            }
            // site-matrix carrier: read the object for its VALUE TYPE
            // only (Lattice::vector_type names the raw SIMD type, not
            // the site object, so coalescedRead is the portable
            // spelling); every element is overwritten below, so no
            // zero-init of prop is needed
            auto pmat = coalescedRead(propR[ss]);
            if ((linear & 0x1) == cbSum) {
              auto v0 = coalescedRead(sumR0[ssh]) * scale;
              auto v1 = coalescedRead(sumR1[ssh]) * scale;
              auto v2 = coalescedRead(sumR2[ssh]) * scale;
              pmat()()(0, 0) = v0()()(0);
              pmat()()(1, 0) = v0()()(1);
              pmat()()(2, 0) = v0()()(2);
              pmat()()(0, 1) = v1()()(0);
              pmat()()(1, 1) = v1()()(1);
              pmat()()(2, 1) = v1()()(2);
              pmat()()(0, 2) = v2()()(0);
              pmat()()(1, 2) = v2()()(1);
              pmat()()(2, 2) = v2()()(2);
            } else {
              auto w0 = (negI * coalescedRead(negR0[ssh])) * scale;
              auto w1 = (negI * coalescedRead(negR1[ssh])) * scale;
              auto w2 = (negI * coalescedRead(negR2[ssh])) * scale;
              pmat()()(0, 0) = w0()()(0);
              pmat()()(1, 0) = w0()()(1);
              pmat()()(2, 0) = w0()()(2);
              pmat()()(0, 1) = w1()()(0);
              pmat()()(1, 1) = w1()()(1);
              pmat()()(2, 1) = w1()()(2);
              pmat()()(0, 2) = w2()()(0);
              pmat()()(1, 2) = w2()()(1);
              pmat()()(2, 2) = w2()()(2);
            }
            coalescedWrite(propW[ss], pmat);
          });
        }
        }

      }
      if (a2aBatch) {
        LOG(Message) << "Reconstructed batch '" << getName() << "' label '"
                     << label << "' slice t=" << t << " (" << par().nNoise
                     << " noise window(s) x " << FImpl::Dimension
                     << " colors) from columns "
                     << par().noiseIndex * FImpl::Dimension << ".."
                     << ((par().noiseIndex + par().nNoise) *
                         FImpl::Dimension - 1)
                     << " of '" << mfs[i] << "'" << std::endl;
      } else {
        LOG(Message) << "Reconstructed '" << mapName(t) << "' label '"
                     << label << "' (" << par().nNoise
                     << " noise window(s)) from columns "
                     << par().noiseIndex * FImpl::Dimension << ".."
                     << ((par().noiseIndex + 1) *
                         FImpl::Dimension - 1)
                     << " of '" << mfs[i] << "'" << std::endl;
      }
    }
  }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // HadronsMILC_MFermion_LMAMesonFieldProp_hpp_
