#pragma once
#ifndef TrackerSD_h
#define TrackerSD_h 1

#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "TrackerHit.hh"
#include "G4AnalysisManager.hh"
#include <vector>
#include <cmath>
#include <numeric>
#include "DetectorInfo.hh"
#include "EventInfo.hh"
#include "TrackerHitCollection.hh"
class G4Step;
class G4HCofThisEvent;



/// Tracker sensitive detector class
///
/// The hits are accounted in hits in ProcessHits() function which is called
/// by Geant4 kernel at each step. A hit is created with each step with non zero
/// energy deposit.

class TrackerSD : public G4VSensitiveDetector
{
  public:
    TrackerSD(G4String name);
    ~TrackerSD() override;

    // methods from base class
    
    void   Initialize(G4HCofThisEvent* hitCollection) override;
    G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;
    void   BeginOfEvent(G4HCofThisEvent* hitCollection);
    void   EndOfEvent(G4HCofThisEvent* hitCollection) override;
    
  
  const float bucket_size1 = 0.100;
  int number_of_buckets1 = (int)ceil(6000/ bucket_size1); // requires <cmath>

    //std::vector<double> Arrivaltime1;
    //std::vector<double> Arrivaltime2;
   //std::vector<double> X;
    std::vector<double> MeanArrivalTime1;
    std::vector<double> MeanArrivalTime2;

  double getAverage1(std::vector<double> const& MeanArrivalTime1){
  if(MeanArrivalTime1.empty()){
    return 0;
  }
  return std::accumulate(MeanArrivalTime1.begin(),MeanArrivalTime1.end(),0.0)/MeanArrivalTime1.size();
}

double getAverage2(std::vector<double> const& MeanArrivalTime2){
  if(MeanArrivalTime2.empty()){
    return 0;
  }
  return std::accumulate(MeanArrivalTime2.begin(),MeanArrivalTime2.end(),0.0)/MeanArrivalTime2.size();
}

  private:
    TrackerHitCollection* fhitCollection;
    
    
    DetectorInfo m_detector1;
    DetectorInfo m_detector2;
    EventInfo m_eventInfo;

  double X_CF_time1;
  double X_CF_time2;
  int x_value1;
  int x_value2;
  int x_value3;
  int x_value4;
  int smallest_element1;
  int largest_element1;
  int smallest_element2;
  int largest_element2;
  double avg1;
  double avg2;
  double Median1;
  double Median2;
  double Max1;
  double Max2;
  int numberOfCounts1;
  int numberOfCounts2;
};



#endif
