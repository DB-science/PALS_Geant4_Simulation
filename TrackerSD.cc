
#include "TrackerSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include <G4ITStepProcessor.hh>
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4SystemOfUnits.hh"
#include "G4VTouchable.hh"
#include "G4OpticalPhoton.hh"
#include "G4Gamma.hh"
#include "DetectorConstruction.hh"
#include "G4ParticleDefinition.hh"
#include <iterator>
#include "TrackerHit.hh"

using namespace std;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerSD::TrackerSD(G4String name)
 : G4VSensitiveDetector(name), 
  fhitCollection(0)
{

  G4String hitsCollectionName;
  collectionName.insert(hitsCollectionName="HitCollect");
 // Arrivaltime1.resize(number_of_buckets1);
  //Arrivaltime2.resize(number_of_buckets1);
//for(int i=0; i< number_of_buckets1; i++){
  //X.push_back(i*bucket_size1);
//}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerSD::~TrackerSD()
{}


void TrackerSD::BeginOfEvent(G4HCofThisEvent*hitCollection )
{
  
/*Arrivaltime1.clear();
Arrivaltime2.clear();
std::fill(Arrivaltime1.begin(),Arrivaltime1.end(),0.0);
  
std::fill(Arrivaltime2.begin(),Arrivaltime2.end(),0.0);

X_CF_time1 = 0;
X_CF_time2 = 0;
x_value1 = 0;
x_value2 = 0;
x_value3 = 0;
x_value4 = 0;
smallest_element1 = 1;
smallest_element2 = 1;
largest_element1 = 0;
largest_element2 = 0;*/
}


void TrackerSD::Initialize(G4HCofThisEvent* hitCollection)
{

  
  // Create hits collection

   fhitCollection 
    = new TrackerHitCollection(SensitiveDetectorName, collectionName[0]);

  // Add this collection in hce

  G4int hcID
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hitCollection->AddHitsCollection( hcID, fhitCollection );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool TrackerSD::ProcessHits(G4Step* aStep,
                                     G4TouchableHistory *ROhist)
{

	static G4ParticleDefinition* gamma =
       G4Gamma::GammaDefinition();
	static G4ParticleDefinition* opticalPhoton =
       G4OpticalPhoton::OpticalPhotonDefinition();
  G4AnalysisManager *AnalysisManager = G4AnalysisManager::Instance();
  
  G4Track *track = aStep->GetTrack();
  track->SetTrackStatus(fStopAndKill);
  G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
  G4StepPoint *postStepPoint = aStep->GetPostStepPoint();


 // TrackerHit* newHit = new TrackerHit();
  

  // Get Detector number 
  const G4VTouchable *touchable = aStep->GetPreStepPoint()->GetTouchable();
  G4int copyNo = touchable->GetVolume()->GetCopyNo();
  
 // fhitCollection->insert( newHit );

  //G4int nofHits = fhitCollection->entries();
  

  //____Detector 0
  if((opticalPhoton)&&(copyNo==0)){
  G4double ArrTime1 =
   aStep->GetPreStepPoint()-> GetGlobalTime()/ns;
    //AnalysisManager->FillNtupleDColumn(0,0,ArrTime1);
    //AnalysisManager->AddNtupleRow(0);
    MeanArrivalTime1.push_back(ArrTime1);
    //G4int nofHits1 = fhitCollection->entries();
   /* int bucket1 = (int)round(ArrTime1 / bucket_size1);
    if(bucket1>=number_of_buckets1){
      bucket1= number_of_buckets1-1;
    }
    Arrivaltime1[bucket1] += 1;
    if((Arrivaltime1[bucket1]>largest_element1)&&(Arrivaltime1[bucket1]>2))
      {
        largest_element1 = Arrivaltime1[bucket1];
        x_value2 = bucket1;
      }*/
      
    }
    
    
    
  
//_____Detector 1
   if((opticalPhoton)&&(copyNo==1)){
   G4double timeOfArrival = aStep->GetPreStepPoint()-> GetGlobalTime()/ns;
    //G4int nofHits2 = fhitCollection->entries();
    //AnalysisManager->FillNtupleDColumn(0,1,timeOfArrival);
    //AnalysisManager->AddNtupleRow(0);
    MeanArrivalTime2.push_back(timeOfArrival);
   /* int bucket2 = (int)round(timeOfArrival / bucket_size1);
    if(bucket2>=number_of_buckets1){
      bucket2= number_of_buckets1-1;
    }
    Arrivaltime2[bucket2] += 1;
    if((Arrivaltime2[bucket2]>largest_element2)&&(Arrivaltime2[bucket2]>2))
      {
        largest_element2 = Arrivaltime2[bucket2];
        x_value4 = bucket2;
      }*/
    
   
   
    }
    
    //AnalysisManager->FillNtupleIColumn(2,6,copyNo);
    //AnalysisManager->AddNtupleRow(2);

   

    //AnalysisManager->FillNtupleIColumn(2,5,m_id);
    //AnalysisManager->AddNtupleRow(2);
  
  
  
  
  
  
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerSD::EndOfEvent(G4HCofThisEvent*hitCollection)
{
//smallest_element1 = 1;
//smallest_element2 = 1;

//Mean Value of Arrivaltime at detector0
if(MeanArrivalTime1.size()>=2){
avg1 = getAverage1(MeanArrivalTime1);
numberOfCounts1 =MeanArrivalTime1.size();
}
else{
  avg1=0;
  numberOfCounts1 =0;
}
//Mean Value of Arrivaltime at detector1
if(MeanArrivalTime2.size()>=2){
avg2 = getAverage2(MeanArrivalTime2);
numberOfCounts2=MeanArrivalTime2.size();
}
else{
  avg2=0;
  numberOfCounts2=0;
}
//std::cout<<"Mean1 is :" <<avg1<<std::endl;
//std::cout<<"Mean2 is :" <<avg2<<std::endl;
  
/*if(MeanArrivalTime1.size()>=2){
    Max1=x_value2*0.01;
}
else{
  Max1=0;
}
if(MeanArrivalTime2.size()>=2){
    Max2=x_value4*0.01;
}
else{
  Max2=0;
}*/
  fhitCollection->SetMean1(avg1);
  fhitCollection->SetMean2(avg2);
  fhitCollection->SetnumberOfCounts1(numberOfCounts1);
  fhitCollection->SetnumberOfCounts2(numberOfCounts2);
 /* fhitCollection->SetMax1(Max1);
  fhitCollection->SetMax2(Max2);
  
  //median 
if(MeanArrivalTime1.size()>=2){
  std::sort(MeanArrivalTime1.begin(),MeanArrivalTime1.end());
  int index1 =MeanArrivalTime1.size()/2;
  Median1 = MeanArrivalTime1[index1];
  
  //if size is pair
  if(!MeanArrivalTime1.size()%2)
  {
    Median1=(Median1+MeanArrivalTime1[index1-1])/2;
  }
}
else{
  Median1 = 0.0;
}
  
  
if(MeanArrivalTime2.size()>=2){
  std::sort(MeanArrivalTime2.begin(),MeanArrivalTime2.end());
  int index2 =MeanArrivalTime2.size()/2;
  Median2 = MeanArrivalTime2[index2];
  
  //if size is pair
  if(!MeanArrivalTime2.size()%2)
  {
    Median2=(Median1+MeanArrivalTime2[index2-1])/2;
  }
}
else{
  Median2 = 0.0;
}

fhitCollection->SetMedian1(Median1);
fhitCollection->SetMedian2(Median2);
*/
    //niedrigesten und höchsten Wert finden

 /*   for (int i=x_value2;i>=0; i--)
    {
      
      if((Arrivaltime1[i]<smallest_element1)&&(Arrivaltime1[i-1]==0)) 
      {
        
        smallest_element1 = Arrivaltime1[i];
        x_value1 = i;
      }
    }
  

    //Mathematik für die Steigung
    double slope1;
    
    //wenn keine Puls dann muss wert 0 sein
    bool isPuls1 = ((largest_element1>0)&&(smallest_element1==0)&&(x_value2>0));
    
    if(isPuls1){
      slope1 = (largest_element1 - smallest_element1)/((x_value2*bucket_size1) - (x_value1*bucket_size1));

      X_CF_time1 = ((0.66*largest_element1)/slope1) + (x_value1*bucket_size1);
      }
    else{
      X_CF_time1 =0.;
    }

    //cout<<"X_CF_1: "<<X_CF_time1<<endl;
    //cout<<"kleinster Wert: "<<smallest_element1<<endl;
    //cout<<"größter Wert: "<< largest_element1<<endl;

 
    for (int i=x_value4; i>=0; i--)
    {
      if((Arrivaltime2[i]<smallest_element2)&&(Arrivaltime2[i-1]==0)) 
      {
        
        smallest_element2 = Arrivaltime2[i];
        x_value3 = i;
        
      }
    }
    //cout<<"größter Wert2: "<< largest_element2<<endl;
    //cout<<"kleinster Wert2: "<<smallest_element2<<endl;
    //cout<<"größter Wert2 X:"<<(x_value4*bucket_size1)<<endl;
    //cout<<"kleinster Wert2 X: "<<(x_value3*bucket_size1)<<endl;

    //Mathematik für die Steigung
    double slope2;
    
    //wenn keine Puls dann muss wert 0 sein
    bool isPuls2 = ((largest_element2>0)&&(smallest_element2==0)&&(x_value4>0));
    
    if (isPuls2){
      slope2 = (largest_element2 - smallest_element2)/((x_value4*bucket_size1) - (x_value3*bucket_size1));

      X_CF_time2 = ((0.66*largest_element2)/slope2) + (x_value3*bucket_size1);
    }
    else
    {
      X_CF_time2 = 0.;
    }

    //cout<<"X_CF_2: "<<X_CF_time2<<endl;

    
    
    
    fhitCollection->SetCF1(X_CF_time1);
    fhitCollection->SetCF2(X_CF_time2);*/


    fhitCollection->insert(new TrackerHit);
    
    

 /* G4int nofHits = fHitsCollection->entries();
  if( nofHits >= 3){
  G4AnalysisManager *AnalysisManager = G4AnalysisManager::Instance();
  AnalysisManager->FillNtupleIColumn(0,0,nofHits);
  AnalysisManager->AddNtupleRow(0);
  AnalysisManager->FillNtupleIColumn(2,0,nofHits);
  AnalysisManager->AddNtupleRow(2);}*/

/*std::fill(Arrivaltime1.begin(),Arrivaltime1.end(),0.0);
  
std::fill(Arrivaltime2.begin(),Arrivaltime2.end(),0.0);

X_CF_time1 = 0;
X_CF_time2 = 0;
x_value1 = 0;
x_value2 = 0;
x_value3 = 0;
x_value4 = 0;
smallest_element1 = 1;
smallest_element2 = 1;
largest_element1 = 0;
largest_element2 = 0;*/
MeanArrivalTime1.clear();
MeanArrivalTime2.clear();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



