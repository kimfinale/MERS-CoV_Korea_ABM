import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;

/**
 * Parameters.java
 * defines input values
 */


public class Parameters {
	private int 				debug = 0;
	private int                 numSampleRun = 1;
	private int                 randomSeed = 3;
	// unit time = day
	private double              stopTime = 60;
	private double              stepSize = 	0.2;
	private double              reportFreq = 1;
	
	// parameter to be estimated
	private double              rateTransmit = 1.5; //  infectious person can transmit infection to a susceptible person with probability of 0.2 
	private double              dayDelayBeforeMovingToAnotherHospital = 1.8; //
		
	private String 				outbreakScenario = "2015";//2015, importation,
	// vaccination related
	private boolean 			underVaccinationScenario = false;
	private String 				vaccinationScenario = "Distance";//Distance, Region, or Hospital
	private double              vaccCoverage = 0.9;
	private double              vaccEfficacy = 0.7;
	//vaccProb* are updated according to vaccEfficacy, vaccCoverage, relativeVaccEfficacyPostExposure, and stepsize
	private double              vaccProbPerStepForSusc = 0.0;
	private double              vaccProbPerStepForExp = 0.0;
	private double              vaccProbPerStep = 0.0;
	
	private double              meanDelayVaccineInducedImmunity = 14.0;//14 
	private double              relativeVaccEfficacyPostExposure = 0.5;
	private int               	thresholdNumberCaseForVaccinationInitiation = 1;
	private int               	thresholdDayVaccinationInitiation = 14;
	private boolean             dayVaccinationStartAdjusted = false; //used as a switch to ensure that vaccination start date is adjusted only once (in response to the number of cases detected)
	private double              vaccinationTargetRadius = 30;//km
	private boolean 			preEmptiveVaccination = false; 
	private boolean             areaTargetedVaccination = false;
	private boolean 			hospitalTargetedVaccination = false;
	private int 				numHospitalsForTargetedVaccination = 5;
	private double              hospitalCoverage = 1.0;
	private double              timeIndexCaseConfirmation = 9.0; // reported time for confirmation of the index case(May 20 after May 11 (onset of symptoms)
	
	private double              shapeGammaOffspring = 0.2; // NegBin (Poisson-Gamma) distribution 
	private double              propSeekingCareFromOtherHospitals = 0.1189189; // 0.1189189; 
	private double              factorHighRiskTransmissibility = 79.0;
	
	// Parameters were estimated from Korean data using coarseDataTools package in R
	private double              shapeGammaDurationOfIncubation = 3.489;	// 
	private double              rateGammaDurationOfIncubation = 0.502;	//
	private double              shapeGammaDurationOfInfectiousness = 5.7;// 
	private double              rateGammaDurationOfInfectiousness = 0.5;//
	
	private double              meanDurationOfInfectiousness = 10.0;//
	private double              maxDurationOfInfectiousness = 13;//
	private double              minDurationOfIncubation = 0;//
	private double              maxDurationOfIncubation = 17;//
	private int              	numStates = 5; // s, e, i, j, r, 
	private int             	numInitPop = 10000;	// total population size
	private int                 numInitInfectious = 1; // number of people infected with MERS-CoV and infectious
	private int                 numInitExposed = 0; // 28 number of people infected with MERS-CoV but not yet infectious
	private int                 numInitIsolated = 0; // number of people isolated
	
	private int 				cumulInc = 0 ; //cumulative incidence is useful to computer the daily incidence 
	private int 				cumulVaccDose = 0 ; //cumulative incidence is useful to computer the daily incidence
	private int 				cumulVaccProtected = 0 ; //cumulative incidence is useful to computer the daily incidence
	
	private double              fracHighInfectivity = 0.11827957; // 22/186 cases, who visited multiple healthcare facilities
	private double              fracReductionLowInfectivity = 0.01265823; //0.1/7.9; // relative reproduction number who remain in one hospital	
	// 2015 outbreak
	private int                 firstGenerationOffspring2015 = 18;
	private int                 minDelaySymptomOnsetToTransmission2015 = 2; // May 11 (symptom onset of the index case), May 15 -  (Pyeongtaek St. Mary's Hospital)
	private int                 maxDelaySymptomOnsetToTransmission2015 = 3; //May 11 (symptom onset of the index case), - May 17 (Pyeongtaek St. Mary's Hospital)
	private int                 numberInfectedsImported = 1;
	// in the primary scenario, Levels 3 and 4 hospitals are included (i.e., cutoffHospitalLevel = 2)
	private int                 cutoffHospitalLevel = 2; // 1=primary, 2=secondary, 3=tertiary, 4=highest level
	private int               	dayIntervention = 29; // intense intervention (significant increase in isolation rate starts at day 29 after the first case developed symptoms
	
	private double              fracVaccTargetPopulation = 1.0;
	
	private double 				timeNeededForVaccination = 10.0;

	private double              dayVaccinationStart = 61.0;
	private int[]               vaccTargetRegionID = { 0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16 };
	// 0=Gangwon, 1=Gyeonggi, 2=Gyeongnam,
	// 3=Gyeongbuk, 4=Gwangju, 5=Daegu,
	// 6=Daejeon, 7=Busan, 8=Seoul,
	// 9=Sejong City, 10=Ulsan, 11=Incheon, 12=Jeonnam, 13=Jeonbuk,
	// 14=Jeju, 15=Chungnam, 16=Chungbuk

	private double              radiusHospitalSearch = 30; //km 
	private double 				timeImmunityDevelopment = 14.0; //time needed for the desired proportion of vaccine recipients, as dictated by vaccine efficacy, is protected 
	//time needed for the desired vaccine coverage is completed
	//patients(n=675), visitors (n=683), healthcare workers (n=218) Cho et al. (2017) Lancet
	private double              fracHCW = 0.14;
	private double              fracPatient = 0.43;
	private double              fracVisitor = 0.43;
	
	private double[] 			meanTimeToIsolation = { 4.259, 2.4, 0.5 }; 
	private double[] 			maxTimeToIsolation = { 13, 5, 0.5 };
	private int[] 				periodCutOff = { 18, 29 };// first 18 days (May 11-28) and the next 11 days (May 29-Jun 9)
	
	private String 				filePathHospSize = "C:/Users/jongHoon.kim/workspace/IVI_Projects/MERS/MERS-CoV_Korea_ABM_Analysis/java/MERS-CoV_Korea_ABM/MERS_Korea_IBM/pop_risk.csv";
	private String 				filePathHospLevel = "C:/Users/jongHoon.kim/workspace/IVI_Projects/MERS/MERS-CoV_Korea_ABM_Analysis/java/MERS-CoV_Korea_ABM/MERS_Korea_IBM/level.csv";
	private String 				filePathHospRegion = "C:/Users/jongHoon.kim/workspace/IVI_Projects/MERS/MERS-CoV_Korea_ABM_Analysis/java/MERS-CoV_Korea_ABM/MERS_Korea_IBM/region_id.csv";
	private String 				filePathHospLongitude = "C:/Users/jongHoon.kim/workspace/IVI_Projects/MERS/MERS-CoV_Korea_ABM_Analysis/java/MERS-CoV_Korea_ABM/MERS_Korea_IBM/longitude.csv";
	private String 				filePathHospLatitude = "C:/Users/jongHoon.kim/workspace/IVI_Projects/MERS/MERS-CoV_Korea_ABM_Analysis/java/MERS-CoV_Korea_ABM/MERS_Korea_IBM/latitude.csv";
	
	//	private double              rateMoveToOtherHospital = 0.05; // 22/186 cases, who visited multiple healthcare facilities, translated to rate of moving to other hospital
	private	int             	numHospital = 100;	// total number of hospitals
	private int             	numPopInHospitalFactor = 1;	// total number of hospitals
	private int                 numPopAtRiskFirstHospital = 716; // population at risk of infection for Pyeongtaik St. Mary's Hospital
	private double              maxNumTransmitPerDay = 200; //  infectious person can transmit infection to a susceptible person with probability of 0.2 
	//	private double              rateIsolate = 0.04; // 1 / 6.0 day since onset of symptoms to isolation
	private double              factorBetaReduceIsolated = 0.0;	// per 13.0 day, i.e,, remain 13.0 days in the hospital on average
//	private double              timelagSymptomOnsetToIsolationInitial = 5.659087;
//	private double              timelagSymptomOnsetToIsolationSlope = 0.133837; // timelag at day t (since the simulation starts) = 5.659087 - 0.133837 * t
//	private double[]            meanDelayFromSymptomOnsetToIsolation = { 3.10, 0.909 };
//	private double[]            maxDelayFromSymptomOnsetToIsolation = { 13, 5 };

	private double              fracRingIsolation = 0.0;
//	private double              timelagSymptomOnsetToTransmission = 1.0; // to reflect that symptom have to become severe before going to a hospital and transmit to others
	private double              timelagSymptomOnsetToTransmissionFirstCase = 4.0; // to reflect that symptom have to become severe before going to a hospital and transmit to others

	
	
	
	
	
	public Parameters(){}
	
	///////////////////////////////////////////////
	// readFile
	// read a file in the designated path
	public ArrayList<Integer> readFile( String location ) {
		ArrayList<Integer> list = new ArrayList<Integer>();
		Scanner s;
		try {
			s = new Scanner( new File( location ));
		    while( s.hasNext() ){
		        if( s.hasNextInt() ){
		            list.add( s.nextInt() );
		        } else {
		            s.next();
		        }
		    }
		    s.close();
//		    System.out.println(list);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	    return ( list );
	}
	
	
	///////////////////////////////////////////////
	// readFile
	// read a file in the designated path
	public ArrayList<Double> readFileDouble( String location ) {
		ArrayList<Double> list = new ArrayList<Double>();
		Scanner s;
		try {
			s = new Scanner( new File( location ));
		    while( s.hasNext() ){
		        if( s.hasNextDouble() ){
		            list.add( s.nextDouble() );
		        } else {
		            s.next();
		        }
		    }
		    s.close();
//		    System.out.println(list);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	    return ( list );
	}
	
	

	

	
	
	// setters and getters
//	public double getTimelagSymptomOnsetToIsolationInitial() {
//		return timelagSymptomOnsetToIsolationInitial;
//	}
//	public void setTimelagSymptomOnsetToIsolationInitial(double timelagSymptomOnsetToIsolationInitial) {
//		this.timelagSymptomOnsetToIsolationInitial = timelagSymptomOnsetToIsolationInitial;
//	}
//	public double getTimelagSymptomOnsetToIsolationSlope() {
//		return timelagSymptomOnsetToIsolationSlope;
//	}
//	public void setTimelagSymptomOnsetToIsolationSlope(double timelagSymptomOnsetToIsolationSlope) {
//		this.timelagSymptomOnsetToIsolationSlope = timelagSymptomOnsetToIsolationSlope;
//	}
	public double[] getMeanTimeToIsolation() {
		return meanTimeToIsolation;
	}
	public void setMeanTimeToIsolation( double[] meanTimeToIsolation ) {
		this.meanTimeToIsolation = meanTimeToIsolation;
	}
	public double[] getMaxTimeToIsolation() {
		return maxTimeToIsolation;
	}
	public void setMaxTimeToIsolation( double[] maxTimeToIsolation ) {
		this.maxTimeToIsolation = maxTimeToIsolation;
	}
	public int[] getPeriodCutOff() {
		return periodCutOff;
	}
	public void setPeriodCutOff(int[] timePeriod) {
		this.periodCutOff = timePeriod;
	}
	public void setRateTransmit( double d ){
		rateTransmit = d;
	}
	public double getRateTransmit(){
		return rateTransmit; 
	}
	public void setFactorBetaReduceIsolated( double d ){
		factorBetaReduceIsolated = d;
	}
	public double getFactorBetaReduceIsolated(){
		return factorBetaReduceIsolated; 
	}
	public void setStopTime( double d ){
		stopTime = d;
	}
	public double getStopTime(){
		return stopTime; 
	}
	public void setStepSize( double d ){
		stepSize = d;
	}
	public double getStepSize(){
		return stepSize; 
	}
	public void setReportFreq( double d ){
		reportFreq = d;
	}
	public double getReportFreq(){
		return reportFreq; 
	}

	public void setNumStates( int i ){
		numStates = i;
	}
	public int getNumStates(){
		return numStates; 
	}

	public void setNumInitPop( int i ){
		numInitPop = i;
	}
	public int getNumInitPop(){
		return numInitPop; 
	}
	public void setNumInitInfectious( int i ){
		numInitInfectious = i;
	}
	public int getNumInitInfectious(){
		return numInitInfectious; 
	}
	public void setNumInitExposed( int i ){
		numInitExposed = i;
	}
	public int getNumInitExposed(){
		return numInitExposed; 
	}
	public void setNumInitIsolated( int i ){
		numInitIsolated = i;
	}
	public int getNumInitIsolated(){
		return numInitIsolated; 
	}
	
	public void setCumulInc( int i ){
		cumulInc = i;
	}
	public int getCumulInc(){
		return cumulInc; 
	}
	public void setShapeGammaDurationOfIncubation( double d ){
		shapeGammaDurationOfIncubation = d;
	}
	public double getShapeGammaDurationOfIncubation(){
		return shapeGammaDurationOfIncubation; 
	}
	public void setRateGammaDurationOfIncubation( double d ){
		rateGammaDurationOfIncubation = d;
	}
	public double getRateGammaDurationOfIncubation(){
		return rateGammaDurationOfIncubation; 
	}

	public void setShapeGammaDurationOfInfectiousness( double d ){
		shapeGammaDurationOfInfectiousness = d;
	}
	public double getShapeGammaDurationOfInfectiousness(){
		return shapeGammaDurationOfInfectiousness; 
	}
	public void setRateGammaDurationOfInfectiousness( double d ){
		rateGammaDurationOfInfectiousness = d;
	}
	public double getRateGammaDurationOfInfectiousness(){
		return rateGammaDurationOfInfectiousness; 
	}
	
	public void setShapeGammaOffspring( double d ){
		shapeGammaOffspring = d;
	}
	public double getShapeGammaOffspring(){
		return shapeGammaOffspring; 
	}
//	public void setRateGammaOffspring( double d ){
//		rateGammaOffspring = d;
//	}
//	public double getRateGammaOffspring(){
//		return rateGammaOffspring; 
//	}
	
	public void setMaxDurationOfIncubation( double d ){
		maxDurationOfIncubation = d;
	}
	public double getMaxDurationOfIncubation(){
		return maxDurationOfIncubation; 
	}
	public void setMaxDurationOfInfectiousness( double d ){
		maxDurationOfInfectiousness = d;
	}
	public double getMaxDurationOfInfectiousness(){
		return maxDurationOfInfectiousness; 
	}
	public void setMeanDurationOfInfectiousness( double d ){
		meanDurationOfInfectiousness = d;
	}
	public double getMeanDurationOfInfectiousness(){
		return meanDurationOfInfectiousness; 
	}
	public void setVaccEfficacy( double d ){
		vaccEfficacy = d;
	}
	public double getVaccEfficacy(){
		return vaccEfficacy; 
	}
	public void setFracHCW( double d ){
		fracHCW = d;
	}
	public double getFracHCW(){
		return fracHCW; 
	}
	public void setFracPatient( double d ){
		fracPatient = d;
	}
	public double getFracPatient(){
		return fracPatient; 
	}
	public void setFracVisitor( double d ){
		fracVisitor = d;
	}
	public double getFracVisitor(){
		return fracVisitor; 
	}

	public void setTimeImmunityDevelopment( double i ){
		timeImmunityDevelopment = i;
	}
	public double getTimeImmunityDevelopment(){
		return timeImmunityDevelopment;
	}
	public double getTimeNeededForVaccination() {
		return timeNeededForVaccination;
	}
	public void setTimeNeededForVaccination( double timeNeededForVaccination) {
		this.timeNeededForVaccination = timeNeededForVaccination;
	}
	public double getVaccCoverage() {
		return vaccCoverage;
	}
	public void setVaccCoverage( double vaccCoverageContact ) {
		this.vaccCoverage = vaccCoverageContact;
	}
	
//	public double getVaccCoverageContact() {
//		return vaccCoverageContact;
//	}
//	public void setVaccCoverageContact(double vaccCoverageContact) {
//		this.vaccCoverageContact = vaccCoverageContact;
//	}
//	public double getVaccCoverageHCW() {
//		return vaccCoverageHCW;
//	}
//	public void setVaccCoverageHCW(double vaccCoverageHCW) {
//		this.vaccCoverageHCW = vaccCoverageHCW;
//	}
//	public double getVaccCoveragePatient() {
//		return vaccCoveragePatient;
//	}
//	public void setVaccCoveragePatient(double vaccCoveragePatient) {
//		this.vaccCoveragePatient = vaccCoveragePatient;
//	}
	public double getFracReductionLowInfectivity() {
		return fracReductionLowInfectivity;
	}
	public void setFracReductionLowInfectivity(double d) {
		this.fracReductionLowInfectivity = d;
	}
	public double getFracHighInfectivity() {
		return fracHighInfectivity;
	}
	public void setFracHighInfectivity(double d) {
		this.fracHighInfectivity = d;
	}
	public double getDayVaccinationStart() {
		return dayVaccinationStart;
	}
	public void setDayVaccinationStart(double dayVaccinationStart) {
		this.dayVaccinationStart = dayVaccinationStart;
	}

	public double getMinDurationOfIncubation() {
		return minDurationOfIncubation;
	}
	public void setMinDurationOfIncubation(double minDurationOfIncubation) {
		this.minDurationOfIncubation = minDurationOfIncubation;
	}
	public double getFracRingIsolation() {
		return fracRingIsolation;
	}
	public void setFracRingIsolation(double fracRingIsolation) {
		this.fracRingIsolation = fracRingIsolation;
	}

	public double getMaxNumTransmitPerDay() {
		return maxNumTransmitPerDay;
	}
	public void setMaxNumTransmitPerDay( double maxNumTransmitPerDay) {
		this.maxNumTransmitPerDay = maxNumTransmitPerDay;
	}
	public int getNumHospital() {
		return numHospital;
	}
	public void setNumHospital(int numHospital) {
		this.numHospital = numHospital;
	}
//	public int getNumPopInHospital() {
//		return numPopInHospital;
//	}
//	public void setNumPopInHospital(int numPopInHospital) {
//		this.numPopInHospital = numPopInHospital;
//	}
//	public double getRateMoveToOtherHospital() {
//		return rateMoveToOtherHospital;
//	}
//	public void setRateMoveToOtherHospital( double rate ) {
//		this.rateMoveToOtherHospital = rate;
//	}
	public double getHospitalCoverage() {
		return hospitalCoverage;
	}
	public void setHospitalCoverage(double coverage) {
		this.hospitalCoverage = coverage;
	}
	public double getRelativeVaccEfficacyPostExposure() {
		return relativeVaccEfficacyPostExposure;
	}
	public void setRelativeVaccEfficacyPostExposure(double relVaccEfficacyPostExposure) {
		this.relativeVaccEfficacyPostExposure = relVaccEfficacyPostExposure;
	}
	public String getFilePathHospSize() {
		return filePathHospSize;
	}
	public void setFilePathHospSize(String filePath) {
		this.filePathHospSize = filePath;
	}
	public String getFilePathHospLevel() {
		return filePathHospLevel;
	}
	public void setFilePathHospLevel(String filePath) {
		this.filePathHospLevel = filePath;
	}
	public String getFilePathHospRegion() {
		return filePathHospRegion;
	}
	public void setFilePathHospRegion(String filePath) {
		this.filePathHospRegion = filePath;
	}
	public int getNumPopInHospitalFactor() {
		return numPopInHospitalFactor;
	}
	public void setNumPopInHospitalFactor(int numPopInHospitalFactor) {
		this.numPopInHospitalFactor = numPopInHospitalFactor;
	}
	public int getRandomSeed() {
		return randomSeed;
	}
	public void setRandomSeed(int randomSeed) {
		this.randomSeed = randomSeed;
	}
	public int getNumPopAtRiskFirstHospital() {
		return numPopAtRiskFirstHospital;
	}
	public void setNumPopAtRiskFirstHospital(int numPopAtRiskFirstHospital) {
		this.numPopAtRiskFirstHospital = numPopAtRiskFirstHospital;
	}
	public double getTimelagSymptomOnsetToTransmissionFirstCase() {
		return timelagSymptomOnsetToTransmissionFirstCase;
	}
	public void setTimelagSymptomOnsetToTransmissionFirstCase( double timelagSymptomOnsetToTransmissionFirstCase ) {
		this.timelagSymptomOnsetToTransmissionFirstCase = timelagSymptomOnsetToTransmissionFirstCase;
	}
	public double getPropSeekingCareFromOtherHospitals() {
		return propSeekingCareFromOtherHospitals;
	}
	public void setPropSeekingCareFromOtherHospitals( double d) {
		this.propSeekingCareFromOtherHospitals = d;
	}
	public int getDayIntervention() {
		return dayIntervention;
	}
	public void setDayIntervention(int dayIntervention) {
		this.dayIntervention = dayIntervention;
	}
//	public double[] getMeanDelayFromSymptomOnsetToIsolation() {
//		return meanDelayFromSymptomOnsetToIsolation;
//	}
//	public void setMeanDelayFromSymptomOnsetToIsolation(double[] meanDelayFromSymptomOnsetToIsolation) {
//		this.meanDelayFromSymptomOnsetToIsolation = meanDelayFromSymptomOnsetToIsolation;
//	}
//	public double[] getMaxDelayFromSymptomOnsetToIsolation() {
//		return maxDelayFromSymptomOnsetToIsolation;
//	}
//	public void setMaxDelayFromSymptomOnsetToIsolation(double[] maxDelayFromSymptomOnsetToIsolation) {
//		this.maxDelayFromSymptomOnsetToIsolation = maxDelayFromSymptomOnsetToIsolation;
//	}

	public double getFactorHighRiskTransmissibility() {
		return factorHighRiskTransmissibility;
	}
	public void setFactorHighRiskTransmissibility( double factorHighRiskTransmissibility) {
		this.factorHighRiskTransmissibility = factorHighRiskTransmissibility;
	}
	public boolean isAreaTargetedVaccination() {
		return areaTargetedVaccination;
	}
	public void setAreaTargetedVaccination(boolean areaTargetedVaccination) {
		this.areaTargetedVaccination = areaTargetedVaccination;
	}
	
	public double getMeanDelayVaccineInducedImmunity() {
		return meanDelayVaccineInducedImmunity;
	}
	public void setMeanDelayVaccineInducedImmunity(double meanDelayVaccineInducedImmunity) {
		this.meanDelayVaccineInducedImmunity = meanDelayVaccineInducedImmunity;
	}
	public double getFracVaccTargetPopulation() {
		return fracVaccTargetPopulation;
	}
	public void setFracVaccTargetPopulation(double fracVaccTargetPopulation) {
		this.fracVaccTargetPopulation = fracVaccTargetPopulation;
	}
	public int[] getVaccTargetRegionID() {
		return vaccTargetRegionID;
	}
	public void setVaccTargetRegionID( int[] vaccTargetRegionID) {
		this.vaccTargetRegionID = vaccTargetRegionID;
	}
	public boolean isPreEmptiveVaccination() {
		return preEmptiveVaccination;
	}
	public void setPreEmptiveVaccination(boolean preEmptiveVaccination) {
		this.preEmptiveVaccination = preEmptiveVaccination;
	}
	public String getFilePathHospLongitude() {
		return filePathHospLongitude;
	}
	public void setFilePathHospLongitude(String filePathHospLongitude) {
		this.filePathHospLongitude = filePathHospLongitude;
	}
	public String getFilePathHospLatitude() {
		return filePathHospLatitude;
	}
	public void setFilePathHospLatitude(String filePathHospLatitude) {
		this.filePathHospLatitude = filePathHospLatitude;
	}
	public double getRadiusHospitalSearch() {
		return radiusHospitalSearch;
	}
	public void setRadiusHospitalSearch(double radiusHospitalSearch) {
		this.radiusHospitalSearch = radiusHospitalSearch;
	}
	public boolean isHospitalTargetedVaccination() {
		return hospitalTargetedVaccination;
	}
	public void setHospitalTargetedVaccination(boolean hospitalTargetedVaccination) {
		this.hospitalTargetedVaccination = hospitalTargetedVaccination;
	}
	public int getNumHospitalsForTargetedVaccination() {
		return numHospitalsForTargetedVaccination;
	}
	public void setNumHospitalsForTargetedVaccination(int numHospitalsForTargetedVaccination) {
		this.numHospitalsForTargetedVaccination = numHospitalsForTargetedVaccination;
	}
	public int getThresholdNumberCaseForVaccinationInitiation() {
		return thresholdNumberCaseForVaccinationInitiation;
	}
	public void setThresholdNumberCaseForVaccinationInitiation( int vaccTargetCaseNumber ) {
		this.thresholdNumberCaseForVaccinationInitiation = vaccTargetCaseNumber;
	}
	public boolean isDayVaccinationStartAdjusted() {
		return dayVaccinationStartAdjusted;
	}
	public void setDayVaccinationStartAdjusted(boolean dayVaccinationStartAdjusted) {
		this.dayVaccinationStartAdjusted = dayVaccinationStartAdjusted;
	}
	public boolean isUnderVaccinationScenario() {
		return underVaccinationScenario;
	}
	public void setUnderVaccinationScenario(boolean underVaccinationScenario) {
		this.underVaccinationScenario = underVaccinationScenario;
	}
	public double getDayDelayBeforeMovingToAnotherHospital() {
		return dayDelayBeforeMovingToAnotherHospital;
	}
	public void setDayDelayBeforeMovingToAnotherHospital(double dayDelayBeforeMovingToAnotherHospital) {
		this.dayDelayBeforeMovingToAnotherHospital = dayDelayBeforeMovingToAnotherHospital;
	}
	public int getFirstGenerationOffspring2015() {
		return firstGenerationOffspring2015;
	}
	public void setFirstGenerationOffspring2015(int firstGenerationOffspring2015) {
		this.firstGenerationOffspring2015 = firstGenerationOffspring2015;
	}
	public int getMinDelaySymptomOnsetToTransmission2015() {
		return minDelaySymptomOnsetToTransmission2015;
	}

	public void setMinDelaySymptomOnsetToTransmission2015(int minDelaySymptomOnsetToTransmission2015) {
		this.minDelaySymptomOnsetToTransmission2015 = minDelaySymptomOnsetToTransmission2015;
	}
	public int getMaxDelaySymptomOnsetToTransmission2015() {
		return maxDelaySymptomOnsetToTransmission2015;
	}
	public void setMaxDelaySymptomOnsetToTransmission2015(int maxDelaySymptomOnsetToTransmission2015) {
		this.maxDelaySymptomOnsetToTransmission2015 = maxDelaySymptomOnsetToTransmission2015;
	}
	public int getNumberInfectedsImported() {
		return numberInfectedsImported;
	}
	public void setNumberInfectedsImported(int numberInfectedsImported) {
		this.numberInfectedsImported = numberInfectedsImported;
	}
	public int getCutoffHospitalLevel() {
		return cutoffHospitalLevel;
	}
	public void setCutoffHospitalLevel( int cutoffHospitalLevel ) {
		this.cutoffHospitalLevel = cutoffHospitalLevel;
	}
	public double getVaccinationTargetRadius() {
		return vaccinationTargetRadius;
	}
	public void setVaccinationTargetRadius( double radius ) {
		this.vaccinationTargetRadius = radius;
	}
	public int getCumulVaccDose() {
		return cumulVaccDose;
	}
	public void setCumulVaccDose(int cumulVaccDose) {
		this.cumulVaccDose = cumulVaccDose;
	}
	public int getNumSampleRun() {
		return numSampleRun;
	}
	public void setNumSampleRun(int numSampleRun) {
		this.numSampleRun = numSampleRun;
	}
	public int getCumulVaccProtected() {
		return cumulVaccProtected;
	}
	public void setCumulVaccProtected(int cumulVaccProtected) {
		this.cumulVaccProtected = cumulVaccProtected;
	}
	public int getThresholdDayVaccinationInitiation() {
		return thresholdDayVaccinationInitiation;
	}
	public void setThresholdDayVaccinationInitiation(int thresholdDayVaccinationInitiation) {
		this.thresholdDayVaccinationInitiation = thresholdDayVaccinationInitiation;
	}
	
	public int getDebug() {
		return debug;
	}
	public void setDebug( int i ) {
		debug = i;
	}
	public String getVaccinationScenario() {
		return vaccinationScenario;
	}
	public void setVaccinationScenario(String vaccinationScenario) {
		this.vaccinationScenario = vaccinationScenario;
	}
	public String getOutbreakScenario() {
		return outbreakScenario;
	}
	public void setOutbreakScenario(String outbreakScenario) {
		this.outbreakScenario = outbreakScenario;
	}
	public double getTimeIndexCaseConfirmation() {
		return timeIndexCaseConfirmation;
	}
	public void setTimeIndexCaseConfirmation(double timeIndexCaseConfirmation) {
		this.timeIndexCaseConfirmation = timeIndexCaseConfirmation;
	}
	public double getVaccProbPerStepForSusc() {
		return vaccProbPerStepForSusc;
	}

	public void setVaccProbPerStepForSusc(double vaccProbPerStepForSusc) {
		this.vaccProbPerStepForSusc = vaccProbPerStepForSusc;
	}

	public double getVaccProbPerStepForExp() {
		return vaccProbPerStepForExp;
	}

	public void setVaccProbPerStepForExp(double vaccProbPerStepForExp) {
		this.vaccProbPerStepForExp = vaccProbPerStepForExp;
	}

	public double getVaccProbPerStep() {
		return vaccProbPerStep;
	}

	public void setVaccProbPerStep(double vaccProbPerStep) {
		this.vaccProbPerStep = vaccProbPerStep;
	}
	
}



