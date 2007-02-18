#ifndef MULTIPLESCATTERINGUPDATOR_H
#define MULTIPLESCATTERINGUPDATOR_H

#include "FastSimulation/MaterialEffects/interface/MaterialEffectsUpdator.h"

/** 
 * This class computes the direction change by multiple scattering 
 * of a charged particle (under the form of a ParticlePropagator, 
 * i.e., a RawParticle) in the tracker layer, and returns the 
 * RawParticle with the modified momentum direction. The tracker 
 * material is assumed to be 100% Si and the Tracker layers are 
 * assumed infinitely thin. The fraction of radiation lengths 
 * traversed by the particle in this tracker layer is determined 
 * in MaterialEffectsUpdator.
 *
 * This version (a la PDG) of a multiple scattering simulator replaces 
 * the buggy GEANT3 Fortran -> C++ former version (up to FAMOS_0_8_0_pre7).
 *
 * \author Patrick Janot
 * $Date: 8-Jan-2004
 */ 

class ParticlePropagator;
class RandomEngine;

class MultipleScatteringUpdator : public MaterialEffectsUpdator
{
 public:

  /// Default Constructor
  MultipleScatteringUpdator(const RandomEngine* engine);

  /// Default Destructor
  ~MultipleScatteringUpdator() {} ;

 private:

  /// The real dE/dx generation and particle update
  void compute(ParticlePropagator &Particle);

};

#endif
