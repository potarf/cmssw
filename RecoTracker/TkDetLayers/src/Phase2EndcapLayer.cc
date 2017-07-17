#include "Phase2EndcapLayer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/GeometrySurface/interface/SimpleDiskBounds.h"

#include "TrackingTools/DetLayers/interface/DetLayerException.h"
#include "TrackingTools/GeomPropagators/interface/HelixForwardPlaneCrossing.h"

#include<array>
#include "DetGroupMerger.h"


//#include "CommonDet/DetLayout/src/DetLessR.h"


using namespace std;

typedef GeometricSearchDet::DetWithState DetWithState;

//hopefully is never called!
const std::vector<const GeometricSearchDet*>& Phase2EndcapLayer::components() const{
  if (not theComponents) {
    auto temp = std::make_unique<std::vector<const GeometricSearchDet*>>();
    temp->reserve(15);   // This number is just an upper bound
    for ( auto c: theComps) temp->push_back(c);
    std::vector<const GeometricSearchDet*>* expected = nullptr;
    if(theComponents.compare_exchange_strong(expected,temp.get())) {
      //this thread set the value
      temp.release();
    }
  }
  return *theComponents;
 }


void
Phase2EndcapLayer::fillRingPars(int i) {
  const BoundDisk& ringDisk = static_cast<const BoundDisk&>(theComps[i]->surface());
  float ringMinZ = std::abs( ringDisk.position().z()) - ringDisk.bounds().thickness()/2.;
  float ringMaxZ = std::abs( ringDisk.position().z()) + ringDisk.bounds().thickness()/2.; 
  RingPar tempPar; 
  tempPar.thetaRingMin =  ringDisk.innerRadius()/ ringMaxZ;
  tempPar.thetaRingMax =  ringDisk.outerRadius()/ ringMinZ;
  tempPar.theRingR=( ringDisk.innerRadius() +
			 ringDisk.outerRadius())/2.;
  ringPars.push_back(tempPar); 
}


Phase2EndcapLayer::Phase2EndcapLayer(vector<const Phase2EndcapRing*>& rings, const bool isOT):
  RingedForwardLayer(true),
  isOuterTracker(isOT),
  theComponents{nullptr}
{
  //They should be already R-ordered. TO BE CHECKED!!
  //sort( theRings.begin(), theRings.end(), DetLessR());

  theRingSize = rings.size();
  LogDebug("TkDetLayers") << "Number of rings in Phase2 OT EC layer is " << theRingSize << std::endl;
  setSurface( computeDisk( rings ) );

  for(unsigned int i=0; i!=rings.size(); ++i) {
    theComps.push_back(rings[i]);
    fillRingPars(i);
    theBasicComps.insert(theBasicComps.end(),	
			 (*rings[i]).basicComponents().begin(),
			 (*rings[i]).basicComponents().end());
  }

 
  LogDebug("TkDetLayers") << "==== DEBUG Phase2EndcapLayer =====" ; 
  LogDebug("TkDetLayers") << "r,zed pos  , thickness, innerR, outerR: " 
			  << this->position().perp() << " , "
			  << this->position().z() << " , "
			  << this->specificSurface().bounds().thickness() << " , "
			  << this->specificSurface().innerRadius() << " , "
			  << this->specificSurface().outerRadius() ;
}


BoundDisk* 
Phase2EndcapLayer::computeDisk( const vector<const Phase2EndcapRing*>& rings) const
{
  float theRmin = rings.front()->specificSurface().innerRadius();
  float theRmax = rings.front()->specificSurface().outerRadius();
  float theZmin = rings.front()->position().z() -
    rings.front()->surface().bounds().thickness()/2;
  float theZmax = rings.front()->position().z() +
    rings.front()->surface().bounds().thickness()/2;
  
  for (vector<const Phase2EndcapRing*>::const_iterator i = rings.begin(); i != rings.end(); i++) {
    float rmin = (**i).specificSurface().innerRadius();
    float rmax = (**i).specificSurface().outerRadius();
    float zmin = (**i).position().z() - (**i).surface().bounds().thickness()/2.;
    float zmax = (**i).position().z() + (**i).surface().bounds().thickness()/2.;
    theRmin = min( theRmin, rmin);
    theRmax = max( theRmax, rmax);
    theZmin = min( theZmin, zmin);
    theZmax = max( theZmax, zmax);
  }
  
  float zPos = (theZmax+theZmin)/2.;
  PositionType pos(0.,0.,zPos);
  RotationType rot;

  return new BoundDisk( pos, rot, new SimpleDiskBounds(theRmin, theRmax,    
				      theZmin-zPos, theZmax-zPos));

}


Phase2EndcapLayer::~Phase2EndcapLayer(){
  for (auto c : theComps) delete c;

  delete theComponents.load();
} 

  

void
Phase2EndcapLayer::groupedCompatibleDetsV( const TrajectoryStateOnSurface& startingState,
				 const Propagator& prop,
				 const MeasurementEstimator& est,
				 std::vector<DetGroup> & result) const
{
  std::array<int,3> const & ringIndices = ringIndicesByCrossingProximity(startingState,prop);
  if ( ringIndices[0]==-1 || ringIndices[1] ==-1 || ringIndices[2] == -1 ) {
    edm::LogError("TkDetLayers") << "TkRingedForwardLayer::groupedCompatibleDets : error in CrossingProximity";
    return;
  }

  //order is: rings in front = 0; rings in back = 1
  //rings should be already ordered in r
  //if the layer has 1 ring, it does not matter
  //FIXME: to be optimized once the geometry is stable
  std::vector<int> ringOrder(theRingSize);
  std::fill(ringOrder.begin(), ringOrder.end(), 1);
  if (theRingSize > 1){
    if(fabs(theComps[0]->position().z()) < fabs(theComps[1]->position().z())){
      for (int i=0; i<theRingSize; i++) {
        if(i % 2 == 0) ringOrder[i] = 0;
      }
    } else if(fabs(theComps[0]->position().z()) > fabs(theComps[1]->position().z())){
      std::fill(ringOrder.begin(), ringOrder.end(), 0);
      for (int i=0; i<theRingSize; i++) {
        if(i % 2 == 0) ringOrder[i] = 1;
      }
    } else {
      throw DetLayerException("Rings in Endcap Layer have same z position, no idea how to order them!");
    }
  }

  auto index = [&ringIndices,& ringOrder](int i) { return ringOrder[ringIndices[i]];};

  std::vector<DetGroup> closestResult;
  theComps[ringIndices[0]]->groupedCompatibleDetsV( startingState, prop, est, closestResult);		
  // if the closest is empty, use the next one and exit: inherited from TID !
  if ( closestResult.empty() ){
    theComps[ringIndices[1]]->groupedCompatibleDetsV( startingState, prop, est, result); 
    return;
  }

  DetGroupElement closestGel( closestResult.front().front());  
  float rWindow = computeWindowSize( closestGel.det(), closestGel.trajectoryState(), est); 

  // check if next ring and next next ring are found and if there is overlap

  bool ring1ok = ringIndices[1] != -1 && overlapInR(closestGel.trajectoryState(),ringIndices[1],rWindow);
  bool ring2ok = ringIndices[2] != -1 && overlapInR(closestGel.trajectoryState(),ringIndices[2],rWindow);

  // look for the two rings in the same plane (are they only two?)

  // determine if we are propagating from in to out (0) or from out to in (1)

  int direction = 0;
  if(startingState.globalPosition().z()*startingState.globalMomentum().z()>0) {
    if(prop.propagationDirection() == alongMomentum) direction=0;
    else direction=1;
  }
  else{
    if(prop.propagationDirection() == alongMomentum) direction=1;
    else direction=0;
  }

  if((index(0) == index(1)) && (index(0) == index(2))) {
    edm::LogInfo("AllRingsInOnePlane") << " All rings: " 
				       << ringIndices[0] << " " 
				       << ringIndices[1] << " " 
				       << ringIndices[2] << " in one plane. Only the first two will be considered";
    ring2ok=false;
  }

  if(index(0) == index(1)) {
    if(ring1ok) {
      std::vector<DetGroup> ring1res;
      theComps[ringIndices[1]]->groupedCompatibleDetsV( startingState, prop, est, ring1res);
      DetGroupMerger::addSameLevel(std::move(ring1res),closestResult);
    }
    if(ring2ok) {
      std::vector<DetGroup> ring2res;
      theComps[ringIndices[2]]->groupedCompatibleDetsV( startingState, prop, est, ring2res);
      DetGroupMerger::orderAndMergeTwoLevels(std::move(closestResult),std::move(ring2res),result,index(0),direction);
      return;
    }
    else {
      result.swap(closestResult);
      return;
    }
  }
  else if(index(0) == index(2)) {
    if(ring2ok) {
      std::vector<DetGroup> ring2res;
      theComps[ringIndices[2]]->groupedCompatibleDetsV( startingState, prop, est, ring2res);
      DetGroupMerger::addSameLevel(std::move(ring2res),closestResult);
    }
    if(ring1ok) {
      std::vector<DetGroup> ring1res;
      theComps[ringIndices[1]]->groupedCompatibleDetsV( startingState, prop, est, ring1res);
      DetGroupMerger::orderAndMergeTwoLevels(std::move(closestResult),std::move(ring1res),result,index(0),direction);
      return;
    }
    else {
      result.swap(closestResult);
      return;
    }
  }
  else {
    std::vector<DetGroup> ring12res;
    if(ring1ok) {
      std::vector<DetGroup> ring1res;
      theComps[ringIndices[1]]->groupedCompatibleDetsV( startingState, prop, est, ring1res);
      ring12res.swap(ring1res);
    } 
    if(ring2ok) {
      std::vector<DetGroup> ring2res;
      theComps[ringIndices[2]]->groupedCompatibleDetsV( startingState, prop, est, ring2res);
      DetGroupMerger::addSameLevel(std::move(ring2res),ring12res);
    }
    if(!ring12res.empty()) {
      DetGroupMerger::orderAndMergeTwoLevels(std::move(closestResult),std::move(ring12res),result,index(0),direction);
      return;
    }
    else {
      result.swap(closestResult);
      return;
    }
  }
}


std::array<int,3> 
Phase2EndcapLayer::ringIndicesByCrossingProximity(const TrajectoryStateOnSurface& startingState,
					 const Propagator& prop ) const
{
  typedef HelixForwardPlaneCrossing Crossing; 
  typedef MeasurementEstimator::Local2DVector Local2DVector;

  HelixPlaneCrossing::PositionType startPos( startingState.globalPosition());
  HelixPlaneCrossing::DirectionType startDir( startingState.globalMomentum());
  PropagationDirection propDir( prop.propagationDirection());
  float rho( startingState.transverseCurvature());

  // calculate the crossings with the ring surfaces
  // rings are assumed to be sorted in R !
  
  Crossing myXing(  startPos, startDir, rho, propDir );

  std::vector<GlobalPoint> ringCrossings;
  ringCrossings.reserve(theRingSize);
  // vector<GlobalVector>  ringXDirections;

  for (int i = 0; i < theRingSize ; i++ ) {
    const BoundDisk & theRing  = static_cast<const BoundDisk &>(theComps[i]->surface());
    pair<bool,double> pathlen = myXing.pathLength( theRing);
    if ( pathlen.first ) { 
      ringCrossings.push_back(GlobalPoint( myXing.position(pathlen.second )));
      // ringXDirections.push_back( GlobalVector( myXing.direction(pathlen.second )));
    } else {
      // TO FIX.... perhaps there is something smarter to do
      //throw DetLayerException("trajectory doesn't cross TID rings");
      ringCrossings.push_back(GlobalPoint( 0.,0.,0.));
      //  ringXDirections.push_back( GlobalVector( 0.,0.,0.));
    }
  }

  //find three closest rings to the crossing

  std::array<int,3> closests = findThreeClosest(ringCrossings);

  return closests;
}




float 
Phase2EndcapLayer::computeWindowSize( const GeomDet* det, 
			     const TrajectoryStateOnSurface& tsos, 
			     const MeasurementEstimator& est) const
{
  const Plane& startPlane = det->surface();  
  MeasurementEstimator::Local2DVector maxDistance = 
    est.maximalLocalDisplacement( tsos, startPlane);
  return maxDistance.y();
}

std::array<int,3>
Phase2EndcapLayer::findThreeClosest(std::vector<GlobalPoint> ringCrossing ) const
{
  std::array<int,3> theBins={{-1,-1,-1}};
  theBins[0] = 0;
  float initialR =  ringPars[0].theRingR;
  float rDiff0 = std::abs( ringCrossing[0].perp() - initialR);
  float rDiff1 = -1.;
  float rDiff2 = -1.;
  for (int i = 1; i < theRingSize ; i++){
    float ringR =  ringPars[i].theRingR;
    float testDiff = std::abs( ringCrossing[i].perp() - ringR);
    if ( testDiff<rDiff0 ) {
      rDiff2 = rDiff1;
      rDiff1 = rDiff0;
      rDiff0 = testDiff;
      theBins[2] = theBins[1];
      theBins[1] = theBins[0];
      theBins[0] = i;
    } 
    else 
    if ( rDiff1 < 0 || testDiff<rDiff1 ) {
      rDiff2 = rDiff1;
      rDiff1 = testDiff;
      theBins[2] = theBins[1];
      theBins[1] = i;
    } 
    else 
    if ( rDiff2 < 0 || testDiff<rDiff2 ) {
      rDiff2 = testDiff;
      theBins[2] = i;
    } 
  }

  return theBins;
}

bool
Phase2EndcapLayer::overlapInR( const TrajectoryStateOnSurface& tsos, int index, double ymax ) const 
{
  // assume "fixed theta window", i.e. margin in local y = r is changing linearly with z
  float tsRadius = tsos.globalPosition().perp();
  float thetamin = ( max(0.,tsRadius-ymax))/(std::abs(tsos.globalPosition().z())+10.f); // add 10 cm contingency 
  float thetamax = ( tsRadius + ymax)/(std::abs(tsos.globalPosition().z())-10.f);
  
  // do the theta regions overlap ?

  return !( thetamin > ringPars[index].thetaRingMax || ringPars[index].thetaRingMin > thetamax);
}


