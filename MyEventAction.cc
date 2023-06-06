

#include "MyEventAction.hh"
#include "EventInfo.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"
#include "G4RunManager.hh"
#include "Randomize.hh"
#include <iomanip>
#include "G4AnalysisManager.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "TrackerHit.hh"
#include "TrackerHitCollection.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyEventAction::MyEventAction(RunAction *runAction, PrimaryGeneratorAction *primaryGenAction): G4UserEventAction(), fHitID(-1)
{
	fEdep1 = 0. ;
	fEdep2 = 0. ;
	fTime0 = 0.*s;
	
	m_primaryGenAction = primaryGenAction;
	m_runAction = runAction;
	
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyEventAction::~MyEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MyEventAction::BeginOfEventAction(const G4Event* evt)
{
	G4int evtNb = evt->GetEventID();
	
	fEdep1 = fEdep2 =0.;
    fTime0 = 0.*s;
    //ArrCF1 = 0;
    //ArrCF2 = 0;
    m_eventInfo.clear();
    m_detector1.clear();
    m_detector2.clear();
    
    m_eventInfo.setInfo(m_primaryGenAction->getTheoLifetime(), m_primaryGenAction->getStartTime(), m_primaryGenAction->getStartAngleX(), m_primaryGenAction->getStartAngleY(),m_primaryGenAction->getStartAngleZ(),
    m_primaryGenAction->getStoppAngle1X(),m_primaryGenAction->getStoppAngle1Y(),m_primaryGenAction->getStoppAngle1Z(),
    m_primaryGenAction->getStoppAngle2X(),m_primaryGenAction->getStoppAngle2Y(), m_primaryGenAction->getStoppAngle2Z());

     G4SDManager * SDman = G4SDManager::GetSDMpointer();
    G4String colNam;
    fHitID = SDman->GetCollectionID(colNam = "HitCollect");
   
}

void MyEventAction::AddEdep(G4int iVol, G4double edep,
                                      G4double time, G4double interactionX, G4double interactionY, G4double interactionZ)           
{
	
    if (fTime0 < 0.) fTime0 = time;
  
  
    if (iVol == 1) { fEdep1 += edep;}
    if (iVol == 2) { fEdep2 += edep;}
    
    //interaction position
    if (iVol == 1) { finterAction1X = interactionX;}
    if (iVol == 2) { finterAction2X = interactionX;}
    if (iVol == 1) { finterAction1Y = interactionY;}
    if (iVol == 2) { finterAction2Y = interactionY;}
    if (iVol == 1) { finterAction1Z = interactionZ;}
    if (iVol == 2) { finterAction2Z = interactionZ;}
	
	
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void MyEventAction::EndOfEventAction(const G4Event* event)
{
    G4AnalysisManager *AnalysisManager = G4AnalysisManager::Instance();
 
   
    G4HCofThisEvent* hitCollection = event->GetHCofThisEvent(); 
    TrackerHitCollection* EHCO =0;
    if (hitCollection){EHCO = (TrackerHitCollection*)(hitCollection->GetHC(fHitID));}
    if(EHCO){int n_hits = EHCO->entries();}
    if(EHCO){
       
        /*ArrCF1 = EHCO->GetCF1();

        ArrCF2 = 0;
        ArrCF2 = EHCO->GetCF2();*/
        
        MeanArr1 = EHCO->GetMean1();
        MeanArr2 = EHCO->GetMean2();
        NumberOfCounts1 = EHCO->GetnumberOfCounts1();
        NumberOfCounts2 = EHCO->GetnumberOfCounts2();
/*
        MedianArr1 = EHCO->GetMedian1();
        MedianArr2 = EHCO->GetMedian2();

        MaxPulse1 = EHCO->GetMax1();
        MaxPulse2 = EHCO->GetMax2();

        */
        //cout<<"MYEVENTACTION Mean1 :"<<MeanArr1<<endl;
        //cout<<"MYEVENTACTION Mean2 :"<<MeanArr2<<endl;
    }
      
  

 
    // threshold in target and detector        
    const G4double Threshold1(40*keV), Threshold2(40*keV);
  
    //coincidence, anti-coincidences 
    //  
    G4bool coincidence       = ((fEdep1 >= Threshold1) && (fEdep2 >= Threshold2));
    
 		 
    G4int evtNb = event->GetEventID();
    

    m_detector1.AddEnergyDep(fEdep1, fTime0, fEdep1 >= Threshold1, evtNb, MeanArr1, finterAction1X, finterAction1Y, finterAction1Z, NumberOfCounts1/*, MedianArr1,MaxPulse1*/); // G4double fEdep, G4double globalTime, G4double ArrCF, bool bHit, G4int evtNb (EventNumber), MeanValue of pulse, Median value of the pulse, Peakvalue X of Pulse
    m_detector2.AddEnergyDep(fEdep2, fTime0, fEdep2 >= Threshold2, evtNb, MeanArr2, finterAction2X, finterAction2Y, finterAction2Z, NumberOfCounts2/*, MedianArr2,MaxPulse2*/);
    
    m_eventInfo.attach(m_detector1,m_detector2);
    
    m_runAction->addEventInfo(m_eventInfo);
    
   
}


G4THitsMap<G4double>*
MyEventAction::GetHitsCollection(G4int hcID,
		const G4Event* event) const
{
	G4THitsMap<G4double>* hitsCollection
		= static_cast<G4THitsMap<G4double>*>(
				event->GetHCofThisEvent()->GetHC(hcID));

	if ( ! hitsCollection ) {
							G4ExceptionDescription msg;
							msg << "Cannot access hitsCollection ID " << hcID;
							G4Exception("B4dEventAction::GetHitsCollection()",
							"MyCode0003", FatalException, msg);
							}

	return hitsCollection;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

