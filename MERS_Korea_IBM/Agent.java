
import java.util.ArrayList;
import java.util.Date;
import java.util.Iterator;
import java.util.Set;


public class Agent {
		// class variables
		static int                         	       nextID = 0;  // to give each an ID
		private String                             infectionStatus = "S"; // "S" means susceptible
		private boolean                            isolated = false; // can also mean quarantine for susceptibles
		private boolean                            quarantined = false; // can also mean quarantine for susceptibles
		private ArrayList<Hospital> 			   visitedHospitals	= new ArrayList<Hospital> (); //I
		private int 			   				   infectedHospitalID = -999; //
		private int                                ID = 0;
		private boolean                            highInfectivity = false;
		private boolean                            invader = false;

		private int                                generation = -999;
		private int                                infectorID = -999;
		private double                             infectivity = -999.0;
		private int                                numOffspring = 0;
		
		private double                             daySinceInfection = -999.0; // this is reset to zero when infected
		private double                             daySinceSymptomOnset = -999.0;// this is reset to zero when symptoms develop
		private double                             daySinceVaccination = -999.0; // this is reset to zero upon vaccination and 
		private double                             durationOfIncubation = -999.0;
		private double                             durationOfInfectiousness = -999.0;
		private double                             delayVaccineInducedImmunity = -999.0;
		private double                             timeToIsolation = -999.0;
		private boolean                            vaccinated = false;

		private Hospital 						   hospital = null;

		private boolean 						   indexCase = false;





	/////////////////////////////////////////////////////////////////////////////////////////////
	// constructors
	// 
	public Agent(){
		ID = nextID ++;
	}

	
	/////////////////////////////////////////////////////////////////////////////////////////////
	// constructors
	public Agent( String state ) {
		ID = nextID++;
		setInfectionStatus( state );
	}
	
	///////////////////////////////////////////////////////////////////////////
	//
	public void gammaDurationOfIncubation( Parameters pars ) {
		double max = pars.getMaxDurationOfIncubation();
		double min = pars.getMinDurationOfIncubation(); 
		double dur = Model.gammaIncubationPeriod.sample();
		if( dur > max ) {
			dur = max;
		} else if( dur < min ) 
			dur = min;
		else {}
		this.setDurationOfIncubation( dur );
	}
	
	///////////////////////////////////////////////////////////////////////////
	//
	public void gammaDelayVaccineInducedImmunity() {
		double max = 21;
		double min = 3; 
		double dur = Model.gammaDelayVaccineInducedImmunity.sample();
		if( dur > max ) {
			dur = max;
		} else if( dur < min ) 
			dur = min;
		else {}
		this.setDelayVaccineInducedImmunity( dur );
	}
		
		
	///////////////////////////////////////////////////////////////////////////
	//
	public void becomeInfectious( Parameters pars ) {
		String status = this.getInfectionStatus();
		if( ! ( status.equals( "E" ) || status.equals( "QE" ) || status.equals( "VE" ) || status.equals( "QVE" ) ) ){
			System.err.println( "Agent.becomeInfectious: infection status have to be E, QE, VE, or QVE. The current state is " + status );
		}
		pars.setCumulInc( pars.getCumulInc() + 1 ); // increase the cumulative incidence by 1
		Model.dayCaseSymptomOnset.add( Step.currentDay );
		this.setInfectionStatus( "I" ); 
		this.setDurationOfInfectiousness( Model.gammaDurationOfInfectiousness.sample() );
		this.setDaySinceSymptomOnset( 0.0 );

	}
		
	
	///////////////////////////////////////////////////////////////////////////
	//
	public void beIsolated( Parameters pars ) {
		String status = this.getInfectionStatus();
		if( ! status.equals( "I" ) ){
			System.err.println( "Agent.beIsolated: infection status have to be I. The current state is " + status );
		}
		this.setInfectionStatus( "J" );
		this.setIsolated( true );
		this.setInfectivity( this.getInfectivity() * pars.getFactorBetaReduceIsolated() );
	}
	
	
	
	
	///////////////////////////////////////////////////////////////////////////
	//
	public void beQuarantined( Parameters pars ) {
		String status = this.getInfectionStatus();
		if( ! ( status.equals( "S" ) || status.equals( "E" ) || status.equals( "VS" ) || status.equals( "VE" )  || status.equals( "VP" )  || status.equals( "R" ) ) ){
			System.err.println( "Agent.beQuarantined: Check the infection state of the input agent! The current state is " + status );
		}
		if( status.equals( "S" ) ) {
			this.setInfectionStatus( "QS");
		} else if( status.equals( "E" ) ) {
			this.setInfectionStatus( "QE");
		}
		else if( status.equals( "VS" ) ) {
			this.setInfectionStatus( "QVS");
		}
		else if( status.equals( "VE" ) ) {
			this.setInfectionStatus( "QVE");
		}
		else if( status.equals( "VP" ) ) {
			this.setInfectionStatus( "QVP");
		}
		else if( status.equals( "R" ) ) {
			this.setInfectionStatus( "JR");
		}
		this.setQuarantined( true );

	}		
	

	
	///////////////////////////////////////////////////////////////////////////
	//
	public void infect( Parameters pars, ArrayList<Agent> offspring ) {
		//increase the number of offsprings by this agent
		this.setNumOffspring( this.getNumOffspring() + offspring.size() );
		double probHighRisk = pars.getPropSeekingCareFromOtherHospitals();
		for( Agent b: offspring ) {
			b.setInfectorID( this.getID() );
			b.setGeneration( this.getGeneration() + 1 );
			b.setDaySinceInfection( 0.0 );
			b.gammaDurationOfIncubation( pars );
			if( Model.unifFromZeroToOne.sample() < probHighRisk ) {
				b.setHighInfectivity( true );
			}
			String str = b.getInfectionStatus();
			if( str.equals("S") ) {
				b.setInfectionStatus( "E" );
			} else if( str.equals("VS") ) {
				b.setInfectionStatus( "VE" );
			}
		}
	}		
		
		
	
	/////////////////////////////////////////////////////////////////
	// printSelf()
	// print information such as ID, degree, risk group,  
	public void printSelf () {
		System.out.println( 
				"ID = " + getID() + ", infection status = " + getInfectionStatus() );
	}


	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	////////////////////////////////////////////////////////////////////////////////////////////
	// setters and getters()
	//
	// note these are class methods, to set class variables
	public void setInfectionStatus ( String s ) {
		infectionStatus = s;
    }
	public String getInfectionStatus () {
		return infectionStatus;
    }
	public void setID ( int i ) {
		ID = i;
    }
	public int getID () {
		return ID;
	}
	public void setInfectorID ( int i ) {
		infectorID = i;
    }
	public int getInfectorID () {
		return infectorID;
	}
	public void setInfectivity ( double d ) {
		infectivity = d;
    }
	public double getInfectivity() {
		return infectivity;
	}
	public void setNumOffspring ( int i ) {
		numOffspring = i;
    }
	public int getNumOffspring() {
		return numOffspring;
	}
	public void setDurationOfIncubation ( double d ) {
		durationOfIncubation = d;
    }
	public double getDurationOfIncubation() {
		return durationOfIncubation;
	}
	public void setDaySinceInfection ( double d ) {
		daySinceInfection = d;
    }
	public double getDaySinceInfection() {
		return daySinceInfection;
	}
	public void setDurationOfInfectiousness ( double d ) {
		durationOfInfectiousness = d;
    }
	public double getDurationOfInfectiousness() {
		return durationOfInfectiousness;
	}
	public double getTimeToIsolation() {
		return timeToIsolation;
	}
	public boolean isIsolated() {
		return isolated;
	}
	public void setIsolated(boolean isolated) {
		this.isolated = isolated;
	}
	public void setTimeToIsolation( double timeToIsolation ) {
		this.timeToIsolation = timeToIsolation;
	}
	public boolean isHighInfectivity() {
		return highInfectivity;
	}
	public void setHighInfectivity( boolean b ) {
		this.highInfectivity = b;
	}
	public Hospital getHospital() {
		return hospital;
	}
	public void setHospital( Hospital hospital ) {
		this.hospital = hospital;
	}
	public double getDaySinceVaccination() {
		return daySinceVaccination;
	}
	public void setDaySinceVaccination(double daySinceVaccination) {
		this.daySinceVaccination = daySinceVaccination;
	}
	public int getGeneration() {
		return generation;
	}
	public void setGeneration(int generation) {
		this.generation = generation;
	}
	public double getDaySinceSymptomOnset() {
		return daySinceSymptomOnset;
	}
	public void setDaySinceSymptomOnset(double daySinceSymptomOnset) {
		this.daySinceSymptomOnset = daySinceSymptomOnset;
	}
	public boolean isQuarantined() {
		return quarantined;
	}
	public void setQuarantined(boolean quarantined) {
		this.quarantined = quarantined;
	}
	public double getDelayVaccineInducedImmunity() {
		return delayVaccineInducedImmunity;
	}
	public void setDelayVaccineInducedImmunity(double delayVaccineInducedImmunity) {
		this.delayVaccineInducedImmunity = delayVaccineInducedImmunity;
	}
	public ArrayList<Hospital> getVisitedHospitals() {
		return visitedHospitals;
	}
	public void setVisitedHospitals(ArrayList<Hospital> visitedHospitals) {
		this.visitedHospitals = visitedHospitals;
	}
	public int getInfectedHospitalID() {
		return infectedHospitalID;
	}
	public void setInfectedHospitalID(int infectedHospitalID) {
		this.infectedHospitalID = infectedHospitalID;
	}
	public boolean isInvader() {
		return invader;
	}
	public void setInvader(boolean invader) {
		this.invader = invader;
	}
	public boolean isIndexCase() {
		return indexCase;
	}
	public void setIndexCase(boolean indexCase) {
		this.indexCase = indexCase;
	}
	public boolean isVaccinated() {
		return vaccinated;
	}
	public void setVaccinated(boolean vaccinated) {
		this.vaccinated = vaccinated;
	}
}
