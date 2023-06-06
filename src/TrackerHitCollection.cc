#include "TrackerHitCollection.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4THitsCollection.hh"


TrackerHitCollection::~TrackerHitCollection() {}

TrackerHitCollection::TrackerHitCollection(G4String detName, G4String colName):G4THitsCollection<TrackerHit>(detName, colName){
    //fCF_Value1= 0.;
    //fCF_Value2=0.;
    fMean1=0.;
    fMean2=0.;
    fnumberOfCounts1 =0;
    fnumberOfCounts2 =0;
   /* fmedian1=0.;
    fmedian2=0.;
    fmax1 = 0.;
    fmax2 = 0.;*/

}
