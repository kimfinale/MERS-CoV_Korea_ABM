import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Set;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;

public class Hospital {
	static int                     nextID 		= 0;// to give each an ID
	private int                    ID 			= 0;
	private int                    level 		= 0; // 1=primary, ..., 4=higher level of hospitals in Korea
	private double                 longitude 	= 0.0;
	private double                 latitude 	= 0.0;
	private int                    regionID 	= -999;
	private double                 dayVaccinationStarted = 999;
	private int                    vaccDosesGiven = 0;
	private boolean                vaccinationImplemented = false;
	private boolean                infectorInvaded = false; //
	private boolean                indexHosp = false;
	
	
	// 0=Gangwon, 1=Gyeonggi, 2=Gyeongnam,
	// 3=Gyeongbuk, 4=Gwangju, 5=Daegu,
	// 6=Daejeon, 7=Busan, 8=Seoul,
	// 9=Sejong City, 10=Ulsan, 11=Incheon, 12=Jeonnam, 13=Jeonbuk, 
	// 14=Jeju, 15=Chungnam, 16=Chungbuk
//	private int                            popSize 		= 0;
	private ArrayList<Agent> 			   infectious 				= new ArrayList<Agent> (); //I
	private ArrayList<Agent> 			   susceptibles 			= new ArrayList<Agent> (); //S
	private ArrayList<Agent> 			   exposeds 				= new ArrayList<Agent> (); //E
	private ArrayList<Agent> 			   isolateds 				= new ArrayList<Agent> (); //J
	private ArrayList<Agent> 			   removeds 				= new ArrayList<Agent> (); //R
	private ArrayList<Agent> 			   isolatedRemoveds 		= new ArrayList<Agent> (); //JR

	private ArrayList<Agent> 			   vaccinatedSusceptibles 	= new ArrayList<Agent> (); //VS
	private ArrayList<Agent> 			   vaccinatedExposeds 		= new ArrayList<Agent> (); //VE
	private ArrayList<Agent> 			   vaccinatedProtecteds 	= new ArrayList<Agent> (); //VP
	
	private ArrayList<Agent> 			   quarantinedSusceptibles 	= new ArrayList<Agent> (); //QS
	private ArrayList<Agent> 			   quarantinedExposeds 		= new ArrayList<Agent> (); //QE
	
	private ArrayList<Agent> 			   quarantinedVaccinatedSusceptibles 	= new ArrayList<Agent> (); //QVS
	private ArrayList<Agent> 			   quarantinedVaccinatedExposeds 		= new ArrayList<Agent> (); //QVE
	private ArrayList<Agent> 			   quarantinedVaccinatedProtecteds 		= new ArrayList<Agent> (); //QVP
	
	private ArrayList<Agent> 			   vaccineReceived 		= new ArrayList<Agent> (); // to track the total number of vaccine doses
	
	
	//because only effective vaccine recipients are moved to the vacinated list, this can cause 
	//overestimation of the cumul vaccine dose
	//So remove vaccine Recipient from the eligible population 
	static Parameters pars = new Parameters ();
	static Utility util = new Utility ();

	/////////////////////////////////////////////////////////////////////////////////////////////
	// constructors
	// 
	public Hospital(){
		ID = nextID ++;
	}
	public Hospital( int N ){
		ID = nextID ++;
		ArrayList<Agent> list = this.getSusceptibles();
		for( int i = 0; i < N; ++i ) {
			Agent a = new Agent( "S" );
			a.getVisitedHospitals().add( this );
			list.add( a );
			a.setHospital( this );
		}
		
	}
	
	public Hospital( int n, String s ){
		ID = nextID ++;
		ArrayList<Agent> list = new ArrayList<Agent>();
		if( s.equals("S") ) {
			list = this.getSusceptibles();
		}
		else if( s.equals("E") ) {
			list = this.getExposeds();
		}
		else if( s.equals("I") ) {
			list = this.getInfectious();
		}
		else if( s.equals("J") ) {
			list = this.getIsolateds();
		}
		for( int i = 0; i < n; ++i ) {
			Agent a = new Agent();
			a.setInfectionStatus( s );
			list.add( a );
			a.setHospital( this );
			a.getVisitedHospitals().add( this );
		}
	}
				

	

	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// transmission( Parameters pars, double beta, double shape )
	//
	public void transmission( Parameters pars, double beta, double shape ) {
		ArrayList<Agent> SToE = new ArrayList<Agent> ();
		ArrayList<Agent> VSToVE = new ArrayList<Agent> ();
	
		ArrayList<Agent> totalSusc = new ArrayList<Agent>();
		totalSusc.addAll( this.getSusceptibles() );
		totalSusc.addAll( this.getVaccinatedSusceptibles() );
		
		int numSusc = totalSusc.size();
		int popSize = this.getPopulationSize();
		int newInfections = 0;

		double factor = pars.getFactorHighRiskTransmissibility();
		boolean transmissionOccurred = false;
		for( Agent a : this.getInfectious() ) {
			if( pars.getDebug() > 2 )
				System.out.printf( "tick = %.1f, transmission begin .. infectious id = %d\n", Step.currentDay, a.getID() );
			double fracSusc = (double) numSusc / popSize; //updated for each infectious individual
			if( fracSusc > 0 ) {
				PoissonDistribution pois = new PoissonDistribution( Model.RNG, beta * fracSusc, 1E-12, 10000000 );
				if( a.isHighInfectivity() && a.isInvader() ) { // we may just need to check invader status as all invaders have to be high-risk
					double scale = beta * fracSusc * factor / shape; // scale for the gamma distribution is by mean divided by shape
					GammaDistribution gamma = new GammaDistribution( Model.RNG, shape, scale );
					double randomMu = gamma.sample();
					if( randomMu <= 0 ) {
						randomMu = Double.MIN_VALUE;
					}
					pois = new PoissonDistribution( Model.RNG, randomMu, 1E-12, 10000000 );
				}
				newInfections = pois.sample();
				if( newInfections > numSusc ) {
					newInfections = numSusc;
				}
				numSusc = numSusc - newInfections;
				if( newInfections > 0 ) {
					ArrayList<Agent> offspring = new ArrayList<Agent>();
					for( int i = 0; i < newInfections; ++i ) {
						Agent b = totalSusc.get( i );
						offspring.add( b );
						b.setInfectedHospitalID( this.getID() );
						if( b.getInfectionStatus().equals("S") ) {
							SToE.add( b );
						} else if( b.getInfectionStatus().equals("VS") ) {
							VSToVE.add( b );
						}
						
						if( pars.getDebug() > 2 )
							System.out.println( "transmission occurred .. infector id =" + a.getID() + ",  infectee id = " + b.getID() );
					}
					totalSusc.removeAll( offspring );
					a.infect( pars, offspring );
					transmissionOccurred = true;
				}
			}
		}
		if( transmissionOccurred ) {
			if( !Model.hospitalsTransmissionOccurred.contains( this ) ) {
				Model.hospitalsTransmissionOccurred.add( this );
			}
		}
		this.getExposeds().addAll( SToE );
		this.getVaccinatedExposeds().addAll( VSToVE );
		this.getSusceptibles().removeAll( SToE );
		this.getVaccinatedSusceptibles().removeAll( VSToVE );

	}		
	
	
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// transmission( Parameters pars )
	//
	public void becomingInfectious( Parameters pars ) {
		ArrayList<Agent> EToI = new ArrayList<Agent> ();
		for( Agent a : this.getExposeds() ) {
			String status = a.getInfectionStatus();
			if( ! ( status.equals( "E" ) ) ){
				System.err.println( "Hospital.becomingInfectious: infection status have to be E. The current status is " + status );
			}
			else if( a.getDurationOfIncubation() <= a.getDaySinceInfection() ) {
				EToI.add( a );
				a.becomeInfectious( pars );
			}
		}
		this.getExposeds().removeAll( EToI ); 
		this.getInfectious().addAll( EToI );	
		
		ArrayList<Agent> VEToI = new ArrayList<Agent> ();	
		for( Agent a : this.getVaccinatedExposeds() ) {
			String status = a.getInfectionStatus();
			if( ! ( status.equals( "VE" ) ) ){
				System.err.println( "Hospital.infectiousnessDevelopment: infection status have to be VE. The current status is " + status );
			}
			else if( a.getDurationOfIncubation() <= a.getDaySinceInfection() ) {
				VEToI.add( a );
				a.becomeInfectious( pars );
			}
		}
		this.getVaccinatedExposeds().removeAll( VEToI ); 
		this.getInfectious().addAll( VEToI );
		
		ArrayList<Agent> QEToJ = new ArrayList<Agent> ();
		for( Agent a : this.getQuarantinedExposeds() ) {
			String status = a.getInfectionStatus();
			if( ! ( status.equals( "QE" ) ) ){
				System.err.println( "Hospital.infectiousnessDevelopment: infection status have to be QE. The current status is " + status );
			}
			else if( a.getDurationOfIncubation() <= a.getDaySinceInfection() ) {
				QEToJ.add( a );
				a.becomeInfectious( pars ); // needed to get the incidence increased although it immediately isolated 
				a.beIsolated( pars );
			}
		}
		this.getQuarantinedExposeds().removeAll( QEToJ ); 
		this.getIsolateds().addAll( QEToJ );
		
		ArrayList<Agent> QVEToJ = new ArrayList<Agent> ();
		for( Agent a : this.getQuarantinedVaccinatedExposeds() ) {
			String status = a.getInfectionStatus();
			if( ! ( status.equals( "QVE" ) ) ){
				System.err.println( "Hospital.infectiousnessDevelopment: infection status have to be QE. The current status is " + status );
			}
			else if( a.getDurationOfIncubation() <= a.getDaySinceInfection() ) {
				QVEToJ.add( a );
				a.becomeInfectious( pars );// needed to get the incidence increased although it immediately isolated 
				a.beIsolated( pars );
			}
		}
		this.getQuarantinedVaccinatedExposeds().removeAll( QVEToJ ); 
		this.getIsolateds().addAll( QVEToJ );
		
	}		
	
	
		
	
	/////////////////////////////////////////////////////////////////////////////////
	// isolation( Parameter pars )
	// agents get isolated some day after becoming symptomatic: I -> J
	protected void isolation( Parameters pars, double rate, double maxTimeFromSymptomOnsetToIsolation  ){
		ArrayList<Agent> agentsIsolated = new ArrayList<Agent> ();
		for( Agent a : this.getInfectious() ) {
			if( !a.getInfectionStatus().equals( "I" ) ){
				System.err.println( "Hospital.isolation: Infection status have to be I. The current state is " + a.getInfectionStatus() );
			}
			else{
				if( !a.isIsolated() && ( ( Model.unifFromZeroToOne.sample() < rate ) ||
						( a.getDaySinceSymptomOnset() >= maxTimeFromSymptomOnsetToIsolation ) ) ){
					a.beIsolated( pars );
					agentsIsolated.add( a );
				}
			}
		}
		if( agentsIsolated.size() > 0 ) {
			this.getIsolateds().addAll( agentsIsolated );
			this.getInfectious().removeAll( agentsIsolated );
			if( !Model.hospitalsCaseIsolated.contains( this ) ) {
				Model.hospitalsCaseIsolated.add( this );
			}
		}
	}
	

	
	
	
	/////////////////////////////////////////////////////////////////////////////////
	// isolation( Parameter pars )
	// agents get isolated some day after becoming symptomatic: I -> J
	protected void removal( Parameters pars ){
		ArrayList<Agent> IToR = new ArrayList<Agent> ();
		for( Agent a : this.getInfectious() ) {
			String status = a.getInfectionStatus();
			if( ! ( status.equals( "I" ) ) ){
				System.err.println( "Hospital.removal: Infection status have to be I. The current status is " + a.getInfectionStatus() );
			}
			if ( a.getDurationOfInfectiousness() <= a.getDaySinceSymptomOnset() ) {
				a.setInfectionStatus( "R" );
				IToR.add( a );
			}
		}
		this.getInfectious().removeAll( IToR );
		this.getRemoveds().addAll( IToR );
		
		ArrayList<Agent> JToJR = new ArrayList<Agent> ();
		for( Agent a : this.getIsolateds() ) {
			String status = a.getInfectionStatus();
			if( ! ( status.equals( "J" ) ) ){
				System.err.println( "Hospital.removal: Infection status have to be J. The current status is " + a.getInfectionStatus() );
			}
			if ( a.getDurationOfInfectiousness() <= a.getDaySinceSymptomOnset() ) {
				a.setInfectionStatus( "JR" );
				JToJR.add( a );
			}
		}
		this.getIsolateds().removeAll( JToJR );
		this.getIsolatedRemoveds().addAll( JToJR );

	}
	

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// quarantine( Parameters pars, double rate, double maxTimeFromSymptomOnsetToIsolation, Hospital hosp ) {
	// quarantine exposed and susceptible 
	// because quarantine is not applied to the vaccinated susceptibles, who are still susceptible to infection, 
	// vaccination may cause increase in the incidence
	public void quarantine( Parameters pars, double rate, double maxTimeFromSymptomOnsetToIsolation ) {
		ArrayList<Agent> quarantined = new ArrayList<Agent> ();
		for( Agent a : this.getExposeds() ) {
			if( !a.isQuarantined() && ( Model.unifFromZeroToOne.sample() < rate ) ){
				a.beQuarantined( pars );
				quarantined.add( a );
			}
		}
		this.getQuarantinedExposeds().addAll( quarantined );
		this.getExposeds().removeAll( quarantined );
		quarantined.clear();
		
		for( Agent a : this.getSusceptibles() ) {
			if( !a.isQuarantined() && ( Model.unifFromZeroToOne.sample() < rate ) ){
				a.beQuarantined( pars );
				quarantined.add( a );
			}
		}
		this.getQuarantinedSusceptibles().addAll( quarantined );
		this.getSusceptibles().removeAll( quarantined );//removed from infection transmission
		quarantined.clear();
		for( Agent a : this.getVaccinatedSusceptibles() ) {
			if( !a.isQuarantined() && ( Model.unifFromZeroToOne.sample() < rate ) ){
				a.beQuarantined( pars );
				quarantined.add( a );
			}
		}
		this.getQuarantinedVaccinatedSusceptibles().addAll( quarantined );
		this.getVaccinatedSusceptibles().removeAll( quarantined ); //removed from infection transmission
		
		quarantined.clear();
		for( Agent a : this.getVaccinatedExposeds() ) {
			if( !a.isQuarantined() && ( Model.unifFromZeroToOne.sample() < rate ) ){
				a.beQuarantined( pars );
				quarantined.add( a );
			}
		}
		this.getQuarantinedVaccinatedExposeds().addAll( quarantined );
		this.getVaccinatedExposeds().removeAll( quarantined ); //removed from infection transmission
		
		quarantined.clear();
		for( Agent a : this.getVaccinatedProtecteds() ) {
			if( !a.isQuarantined() && ( Model.unifFromZeroToOne.sample() < rate ) ){
				a.beQuarantined( pars );
				quarantined.add( a );
			}
		}
		this.getQuarantinedVaccinatedProtecteds().addAll( quarantined );
		this.getVaccinatedProtecteds().removeAll( quarantined ); 
		
		quarantined.clear();
		for( Agent a : this.getRemoveds() ) {
			if( !a.isQuarantined() && ( Model.unifFromZeroToOne.sample() < rate ) ){
				a.beQuarantined( pars );
				quarantined.add( a );
			}
		}
		this.getIsolatedRemoveds().addAll( quarantined );
		this.getRemoveds().removeAll( quarantined ); 
		
		if( pars.getDebug() > 0 ) {
			System.out.printf( "Day: %.1f, Step.quarantine occurred.\n", Step.currentDay );
			for(Hospital h : Model.hospitals ) {
				h.printSelf();
			}
		}
	}


	
	//////////////////////////////////////////////////////////
	// vaccination(  double vaccProbSusc, double vaccProbExp )
	//
	public void vaccination( Parameters pars, double vaccProbPerStepSusc, double vaccProbPerStepExp, double vaccReceivingProbPerStep ){
		ArrayList<Agent> agentVaccinated = new ArrayList<Agent>();
		ArrayList<Agent> S = this.getSusceptibles();
		ArrayList<Agent> E = this.getExposeds();
		ArrayList<Agent> QS = this.getQuarantinedSusceptibles();
		ArrayList<Agent> QE = this.getQuarantinedExposeds();
		// Ensure that people get vaccinated only once for the correct number of vaccines, i.e., account for those who receive vaccine and don't become immune
		ArrayList<Agent> vaccReceivedAlready = this.getVaccineReceived() ;
		S.removeAll( vaccReceivedAlready ); 
		E.removeAll( vaccReceivedAlready );
		QS.removeAll( vaccReceivedAlready );
		QE.removeAll( vaccReceivedAlready );
		
		for( Agent a : S ) {
			double u = Model.unifFromZeroToOneVacc.sample(); 
			if( u < vaccProbPerStepSusc ) {
				a.setDaySinceVaccination( 0.0 );
				a.setInfectionStatus( "VS" ); // efficacy is interpreted as being the fraction to be 
				agentVaccinated.add( a );
			}
			if( u < vaccReceivingProbPerStep ) {
				vaccReceivedAlready.add( a );
				pars.setCumulVaccDose( pars.getCumulVaccDose() + 1 );
			}
		}
		this.getVaccinatedSusceptibles().addAll( agentVaccinated );
		S.removeAll( agentVaccinated );
		
		agentVaccinated.clear();
		for( Agent a : E ) {
			double u = Model.unifFromZeroToOneVacc.sample();
			if( u < vaccProbPerStepExp ) {
				a.setDaySinceVaccination( 0.0 );
//				a.gammaDelayVaccineInducedImmunity();
				a.setInfectionStatus( "VE" );
				agentVaccinated.add( a );
			}
			if( u < vaccReceivingProbPerStep ) {
				vaccReceivedAlready.add( a );
				pars.setCumulVaccDose( pars.getCumulVaccDose() + 1 );
			}
		}
		this.getVaccinatedExposeds().addAll( agentVaccinated );
		E.removeAll( agentVaccinated );
		
		agentVaccinated.clear();
		for( Agent a : QS ) {
			double u1 = Model.unifFromZeroToOneVacc.sample(); 
			if( u1 < vaccProbPerStepSusc ) {
				a.setDaySinceVaccination( 0.0 );
//				a.gammaDelayVaccineInducedImmunity();
				a.setInfectionStatus( "QVS" ); // efficacy is interpreted as being the fraction to be 
				agentVaccinated.add( a );
				
			}
			if( u1 < vaccReceivingProbPerStep ) {
				vaccReceivedAlready.add( a );
				pars.setCumulVaccDose( pars.getCumulVaccDose() + 1 );
			}
		}
		this.getQuarantinedVaccinatedSusceptibles().addAll( agentVaccinated );
		QS.removeAll( agentVaccinated );
		
		agentVaccinated.clear();
		for( Agent a : QE ) {
			double u2 = Model.unifFromZeroToOneVacc.sample();
			if( u2 < vaccProbPerStepExp ) {
				a.setDaySinceVaccination( 0.0 );
//				a.gammaDelayVaccineInducedImmunity();
				a.setInfectionStatus( "QVE" );
				agentVaccinated.add( a );
			}
			if( u2 < vaccReceivingProbPerStep ) {
				vaccReceivedAlready.add( a );
				pars.setCumulVaccDose( pars.getCumulVaccDose() + 1 );
			}
		}
		this.getQuarantinedVaccinatedExposeds().addAll( agentVaccinated );
		QE.removeAll( agentVaccinated );
	}
	
	
	

	
	////////////////////////////////////////////////////////////////////////////////
	// selectHospital()
	// select one hospital to move arch hospitals within certain distance
	//
	public Hospital selectHospital( Parameters pars ){
		ArrayList<Hospital> targetHosp = new ArrayList<Hospital>();
		ArrayList<Hospital> hospPool = new ArrayList<Hospital>();
		double distCutoff = pars.getRadiusHospitalSearch();
		double myLon = this.getLongitude();
		double myLat = this.getLatitude();
		hospPool.addAll( Model.hospitals );
		hospPool.addAll( Model.uninfectedHospitals );
		hospPool.removeAll( Model.hospitalsCaseIsolated );
		// case can go any hospital from the Level 4 hospitals in Seoul
		if( this.getLevel() == 4 && this.getRegionID() == 8 ) { 
			targetHosp.addAll( Model.uninfectedHospitals );
		}
		else {
			for( Hospital h : hospPool ) {
				double d = util.getDistance( myLat, h.getLatitude(), myLon, h.getLongitude() ); 
				if( d < distCutoff || ( h.getLevel() == 4 && h.getRegionID() == 8 ) ) {// Level 4 Hospitals in Seoul are always accessible 
					targetHosp.add( h );
				}
			}
		}		
		if( targetHosp.size() > 0 ) {// no target hospital then stay where you are			 
			int N = targetHosp.size();
			int[] popCumSum = new int[ N ];
			int popsize = targetHosp.get(0).getPopulationSize();
			popCumSum[ 0 ] = popsize * popsize;
	        for( int i = 1; i < N; i++ ) {
	        	int popsize1 = targetHosp.get(i).getPopulationSize();
	        	popCumSum[ i ] = popCumSum[ i-1  ] +  popsize1 * popsize1;
	        }
	        int index = 0;
	        double r = Model.unifFromZeroToOne.sample();
	        for( int i = 0; i < N; ++ i ) {       	
	        	if( r <= ( (double) popCumSum[i] / popCumSum[N-1] ) ) {
	        		index = i;
	        		break;
	        	}
	        }
		    return targetHosp.get( index );
		} else {
			return null;
		}
	}
	
	
	
/*	
	////////////////////////////////////////////////////////////////////////////////
	// selectHospital()
	// select one hospital to move arch hospitals within certain distance
	//
	public Hospital selectHospital( Parameters pars ){
		ArrayList<Hospital> targetHospList = new ArrayList<Hospital>();
		double distCutoff = pars.getRadiusHospitalSearch();
		double myPop = this.getPopSize();
		double myLon = this.getLongitude();
		double myLat = this.getLatitude();
		int iter = Model.uninfectedHospitals.size(); 
		for( int i = 0; i < iter; ++i ) {
			Hospital h = Model.uninfectedHospitals.get( i );
			double d = util.getDistance( myLat, h.getLatitude(), myLon, h.getLongitude() ); 
//			if( d < distCutoff && myPop < h.getPopSize() ) {
			if( d < distCutoff || ( h.getLevel() == 4 && h.getRegionID() == 8 ) ) {// Level 4 Hospitals in Seoul are always accessible 
				targetHospList.add( h );
			}
		}
				
//		for( int k = 0; k < targetHospList.size(); ++k ) {
//			System.out.println( "time = " + Step.currentDay + ", "
//		
//					+ "level = " +  targetHospList.get(k).getLevel() 
//					+ ", regionID = " +  targetHospList.get(k).getRegionID() 
//					+ ", pop at risk = " +  targetHospList.get(k).getPopSize() );
//		}
		// select hospital in favor of those with larger population (in fact, square of population size)
		if( targetHospList.size() > 0 ) {// no target hospital then stay where you are			 
			int N = targetHospList.size();
			int[] popCumSum = new int[ N ];
			int popsize = targetHospList.get(0).getPopSize();
			popCumSum[ 0 ] = popsize * popsize;
	        for( int i = 1; i < N; i++ ) {
	        	int popsize1 = targetHospList.get(i).getPopSize();
	        	popCumSum[ i ] = popCumSum[ i-1  ] +  popsize1 * popsize1;
	        }
	            
	        int index = 0;
	        double r = Model.unifFromZeroToOne.sample();
	        for( int i = 0; i < N; i++ ) {       	
	        	if( r <= ( (double) popCumSum[i] / popCumSum[N-1] ) ) {
	        		index = i;
	        		break;
	        	}
	        }
//			
//	        System.out.println( "time = " + Step.currentDay + ", selected id = " +  targetHospList.get(index).getID() );
			
		    return targetHospList.get(index);
		} else {
			return null;
		}
	}
	
*/	
	
	

	////////////////////////////////////////////////////////////////////////////////
	// selectHospitalsForVaccinationByDistance()
	// select hospitals within a prespecified distance from the hospitals where cases are detected
	//
	public ArrayList<Hospital> selectHospitalsForVaccinationByDistance( Parameters pars ){
		Set<Hospital> vaccTargetHospList = new HashSet<Hospital>();
		double distCutoff = pars.getVaccinationTargetRadius();
		ArrayList<Hospital> allHospitals = new ArrayList<Hospital>();
		allHospitals.addAll( Model.hospitals );
		allHospitals.addAll( Model.uninfectedHospitals );
		int iter = allHospitals.size(); 
		for( Hospital hosp : Model.hospitals ) {
			// vaccination happens in the hospitals that are within a certain distance from the hospital where cases were isolated			
			if( hosp.getIsolateds().size() > 0 ) { 
				for( int i = 0; i < iter; ++i ) {
					double myLon = hosp.getLongitude();
					double myLat = hosp.getLatitude();
					Hospital h = allHospitals.get(i);
					double d = util.getDistance( myLat, h.getLatitude(), myLon, h.getLongitude() ); 
					if( d < distCutoff ) {
						vaccTargetHospList.add( h );
					}
				}
			}
		}
		return (ArrayList<Hospital>) vaccTargetHospList;
	}

/*		

	////////////////////////////////////////////////////////////////////////////////
	// selectHospitalsForVaccinationByDistance()
	// select hospitals within a prespecified distance from the hospitals where cases are detected
	//
	public ArrayList<Hospital> selectHospitalsForVaccinationByDistance( Parameters pars ){
		ArrayList<Hospital> targetHospList = new ArrayList<Hospital>();
		double distCutoff = pars.getVaccinationTargetRadius();
		double myLon = this.getLongitude();
		double myLat = this.getLatitude();
		ArrayList<Hospital> allHospitals = new ArrayList<Hospital>();
		allHospitals.addAll( Model.hospitals );
		allHospitals.addAll( Model.uninfectedHospitals );
		int iter = allHospitals.size(); 
		for( int i = 0; i < iter; ++i ) {
			Hospital h = allHospitals.get(i);
			double d = util.getDistance( myLat, h.getLatitude(), myLon, h.getLongitude() ); 
			if( d < distCutoff ) {
				targetHospList.add( h );
			}
		}
		return targetHospList;
		
	}

		
*/
	////////////////////////////////////////////////////////////////////////////////
	// electHospitalsForVaccination()
	// select hospitals within a fixed distance from the hospitals where cases are detected
	//
	public ArrayList<Hospital> selectHospitalsForVaccination( Parameters pars ){
		ArrayList<Hospital> targetHospList = new ArrayList<Hospital>();
		double distCutoff = pars.getRadiusHospitalSearch();
		double myLon = this.getLongitude();
		double myLat = this.getLatitude();
		int iter = Model.uninfectedHospitals.size(); 
		for( int i = 0; i < iter; ++i ) {
			Hospital h = Model.uninfectedHospitals.get(i);
			double d = util.getDistance( myLat, h.getLatitude(), myLon, h.getLongitude() ); 
			if( d < distCutoff ) {
				targetHospList.add( h );
			}
		}
		Collections.sort( targetHospList, new SortByPop() );
		if( targetHospList.size() > 0 && targetHospList.size() <= 5 ) {
			return targetHospList;
		}
		else if( targetHospList.size() > 5 ) {
			ArrayList<Hospital> list = new ArrayList<Hospital>();
			for( int r = 0; r < 5; ++r ) {
				list.add( targetHospList.get(r) );
			}
			return list;
		} else {
			return null;
		}
	}
	
	
	////////////////////////////////////////////////////////////////////
	class SortByPop implements Comparator<Hospital> { 
	    public int compare(Hospital a, Hospital b) { 
	        return b.getPopulationSize() - a.getPopulationSize(); 
	    } 
	}
	
	//////////////////////////////////////////////////////////////////////////////
	//
	public ArrayList<Hospital> getTargetHospitalList(Parameters pars){
		ArrayList<Hospital> targetHospList = new ArrayList<Hospital>();
		double distCutoff = pars.getRadiusHospitalSearch();
		double myPop = this.getPopulationSize();
		double myLon = this.getLongitude();
		double myLat = this.getLatitude();
		int iter = Model.uninfectedHospitals.size(); 
		for( int i = 0; i < iter; ++i ) {
			Hospital h = Model.uninfectedHospitals.get(i);
			double d = util.getDistance( myLat, h.getLatitude(), myLon, h.getLongitude() ); 
			if( d < distCutoff  && myPop < h.getPopulationSize() ) {
				targetHospList.add( h );
			}
		}
		return targetHospList;
	}
	
	
	/////////////////////////////////////////////////////////////////////////////////////////////
	// getFracSusceptibles()
	// return the fraction of susceptibles, which is used to proportionally reduce the number of secondary cases
	// quarantined susceptibles are excluded
	public double getFracSusceptibles(){
		int numS = this.getSusceptibles().size();
		int numVS = this.getVaccinatedSusceptibles().size();
		int susc = numS + numVS;
		
		double frac = (double) susc / this.getPopulationSize() ;
		return frac;
	}
	

	

	/////////////////////////////////////////////////////////////////////////////////////////////
	// retrieveExposedAgentsFromHospital()
	// agents exposed to either vaccine or the MERS-CoV
	public ArrayList<Agent> retrieveExposedAgentsFromHospital(){
		ArrayList<Agent> list = new ArrayList<Agent>();
		// exposed to the infected persons
		list.addAll( this.getExposeds() );		
		list.addAll( this.getInfectious() );
		list.addAll( this.getIsolateds() );
		list.addAll( this.getQuarantinedExposeds() );
		list.addAll( this.getQuarantinedVaccinatedExposeds() );
		
		return list;
	}
	

	/////////////////////////////////////////////////////////////////////////////////////////////
	// getPopulationSize()
	// return the fraction of susceptibles, which is used to proportionally reduce the number of secondary cases
	// quarantined susceptibles are excluded
	public int getPopulationSize(){
		int numS = this.getSusceptibles().size();
		int numE = this.getExposeds().size();
		int numI = this.getInfectious().size();
		int numJ = this.getIsolateds().size();
		int numR = this.getRemoveds().size();
		int numJR = this.getIsolatedRemoveds().size();
		
		int numVS = this.getVaccinatedSusceptibles().size();
		int numVE = this.getVaccinatedExposeds().size();
		int numVP = this.getVaccinatedProtecteds().size(); // 
		
		int numQS = this.getQuarantinedSusceptibles().size();
		int numQE = this.getQuarantinedExposeds().size(); // 
		int numQVS = this.getQuarantinedVaccinatedSusceptibles().size();
		int numQVE = this.getQuarantinedVaccinatedExposeds().size();
		int numQVP = this.getQuarantinedVaccinatedProtecteds().size(); // 
		
		int pop = numS + numE + numI + numJ + numR + numJR + numVS + numVE + numVP + numQS + numQE + numQVS + numQVE + numQVP;
				
		return pop;
	}
	
	public void setPopulationDistributionPostVaccination( Parameters pars) {
		// dS/dt = -r*S
		// dVS/dt = +r*S - g*V
		// dVP/dt = g*V
		// S(t) = exp(-r*t)
		// VS(t) = (r*exp(t*(-g-r))*(exp(g*t)-exp(r*t)))/(g-r)
		// VP(t) = (exp(t*(-g-r))*(exp(r*t)*((g-r)*exp(g*t)+r)-g*exp(g*t)))/(g-r)
		double t = Step.currentDay - this.getDayVaccinationStarted();
		double r = pars.getVaccProbPerStepForSusc();
		double g = pars.getStepSize()/pars.getMeanDelayVaccineInducedImmunity();
		double VS = (r*Math.exp(t*(-g-r))*(Math.exp(g*t)-Math.exp(r*t)))/(g-r);
		double VP = (Math.exp(t*(-g-r))*(Math.exp(r*t)*((g-r)*Math.exp(g*t)+r)-g*Math.exp(g*t)))/(g-r);
		double S = Math.exp( - pars.getVaccProbPerStep() * t );
		ArrayList<Agent> susc = getSusceptibles();
		ArrayList<Agent> vaccSusc = getVaccinatedSusceptibles();
		ArrayList<Agent> vaccProtected = getVaccinatedProtecteds();
		ArrayList<Agent> vaccReceived = getVaccineReceived();
		int numVS = (int) ( susc.size() * VS );
		int numVP = (int) ( susc.size() * VP );
		int numVR = (int) ( susc.size() * (1-S) );
//		System.out.println( "Day="+ Step.currentDay  + ", t=" + t + ", id="+ getID() + ", numVS="+ numVS + ", numVP=" + numVP + ", numVR=" + numVR );
		for( int i =0; i < numVR; ++i ) {
			Agent a = susc.get( i ); 
			vaccReceived.add( a );
			if( i < numVS ) {
				a.setInfectionStatus( "VS" );
				vaccSusc.add( a );
			}
			else if( i < (numVS+numVP) ){
				vaccProtected.add( a );
				a.setInfectionStatus( "VP" );
			}
		}
		susc.removeAll( vaccSusc );
		susc.removeAll( vaccProtected );
		
	}
		

		
		
	/////////////////////////////////////////////////////////////////////////////////////////////
	// getAgentsFromHospital()
	// 
	public ArrayList<Agent> retrieveAgentsFromHospital(){
		ArrayList<Agent> list = new ArrayList<Agent>();
		list.addAll( this.getSusceptibles() );
		list.addAll( this.getExposeds() );
		list.addAll( this.getInfectious() );
		list.addAll( this.getIsolateds() );
		list.addAll( this.getRemoveds() );
		list.addAll( this.getIsolatedRemoveds() );
		list.addAll( this.getVaccinatedSusceptibles() );
		list.addAll( this.getVaccinatedExposeds() );
		list.addAll( this.getVaccinatedProtecteds() ); // 
		list.addAll( this.getQuarantinedSusceptibles() );
		list.addAll( this.getQuarantinedExposeds() ); // 
		list.addAll( this.getQuarantinedVaccinatedSusceptibles() );
		list.addAll( this.getQuarantinedVaccinatedExposeds() );
		list.addAll( this.getQuarantinedVaccinatedProtecteds() ); // 
				
		return list;
	}
	
	
	/////////////////////////////////////////////////////////////////
	// printSelf()
	// print information such as ID, degree, risk group,  
	public void printSelf () {
		for( Agent a : getExposeds() ) {
			System.out.println( "Exposed list: " + a.getInfectionStatus() );
		}
		for( Agent a : getInfectious() ) {
			System.out.println( "Infectious list: " + a.getInfectionStatus() );
		}
		
		for( Agent a : getIsolateds() ) {
			System.out.println( "Isolated list: " + a.getInfectionStatus() );
		}
	}


	////////////////////////////////////////////////////////////////////////////////////////////
	// setters and getters()
	// note these are class methods, to set class variables
	public void setID ( int i ) {
		ID = i;
    }
	public int getID () {
		return ID;
	}
	public int getRegionID() {
		return regionID;
	}
	public void setRegionID( int id) {
		this.regionID = id;
	}
	public double getLongitude() {
		return longitude;
	}
	public void setLongitude(double longitude) {
		this.longitude = longitude;
	}
	public double getLatitude() {
		return latitude;
	}
	public void setLatitude(double latitude) {
		this.latitude = latitude;
	}
	public int getLevel() {
		return level;
	}
	public void setLevel(int level) {
		this.level = level;
	}
	
	public ArrayList<Agent> getInfectious() {
		return infectious;
	}
	public void setInfectious(ArrayList<Agent> infectious) {
		this.infectious = infectious;
	}
	public ArrayList<Agent> getSusceptibles() {
		return susceptibles;
	}
	public void setSusceptibles(ArrayList<Agent> susceptibles) {
		this.susceptibles = susceptibles;
	}
	public ArrayList<Agent> getExposeds() {
		return exposeds;
	}
	public void setExposeds(ArrayList<Agent> exposeds) {
		this.exposeds = exposeds;
	}
	public ArrayList<Agent> getIsolateds() {
		return isolateds;
	}
	public void setIsolateds(ArrayList<Agent> isolateds) {
		this.isolateds = isolateds;
	}
	public ArrayList<Agent> getRemoveds() {
		return removeds;
	}
	public void setRemoveds(ArrayList<Agent> removeds) {
		this.removeds = removeds;
	}
	public ArrayList<Agent> getVaccinatedSusceptibles() {
		return vaccinatedSusceptibles;
	}
	public void setVaccinatedSusceptibles(ArrayList<Agent> vaccinatedSusceptibles) {
		this.vaccinatedSusceptibles = vaccinatedSusceptibles;
	}
	public ArrayList<Agent> getVaccinatedExposeds() {
		return vaccinatedExposeds;
	}
	public void setVaccinatedExposeds(ArrayList<Agent> vaccinatedExposeds) {
		this.vaccinatedExposeds = vaccinatedExposeds;
	}
	public ArrayList<Agent> getVaccinatedProtecteds() {
		return vaccinatedProtecteds;
	}
	public void setVaccinatedProtecteds(ArrayList<Agent> vaccinatedProtecteds) {
		this.vaccinatedProtecteds = vaccinatedProtecteds;
	}
	public ArrayList<Agent> getQuarantinedSusceptibles() {
		return quarantinedSusceptibles;
	}
	public void setQuarantinedSusceptibles(ArrayList<Agent> quarantinedSusceptibles) {
		this.quarantinedSusceptibles = quarantinedSusceptibles;
	}
	public ArrayList<Agent> getQuarantinedExposeds() {
		return quarantinedExposeds;
	}
	public void setQuarantinedExposeds(ArrayList<Agent> quarantinedExposeds) {
		this.quarantinedExposeds = quarantinedExposeds;
	}
	public ArrayList<Agent> getVaccineReceived() {
		return vaccineReceived;
	}
	public void setVaccineReceived(ArrayList<Agent> vaccineReceived) {
		this.vaccineReceived = vaccineReceived;
	}
	public double getDayVaccinationStarted() {
		return dayVaccinationStarted;
	}
	public void setDayVaccinationStarted(double dayVaccinationStarted) {
		this.dayVaccinationStarted = dayVaccinationStarted;
	}
	public ArrayList<Agent> getIsolatedRemoveds() {
		return isolatedRemoveds;
	}
	public void setIsolatedRemoveds(ArrayList<Agent> isolatedRemoveds) {
		this.isolatedRemoveds = isolatedRemoveds;
	}
	public ArrayList<Agent> getQuarantinedVaccinatedSusceptibles() {
		return quarantinedVaccinatedSusceptibles;
	}
	public void setQuarantinedVaccinatedSusceptibles(ArrayList<Agent> quarantinedVaccinatedSusceptibles) {
		this.quarantinedVaccinatedSusceptibles = quarantinedVaccinatedSusceptibles;
	}
	public ArrayList<Agent> getQuarantinedVaccinatedExposeds() {
		return quarantinedVaccinatedExposeds;
	}
	public void setQuarantinedVaccinatedExposeds(ArrayList<Agent> quarantinedVaccinatedExposeds) {
		this.quarantinedVaccinatedExposeds = quarantinedVaccinatedExposeds;
	}
	public ArrayList<Agent> getQuarantinedVaccinatedProtecteds() {
		return quarantinedVaccinatedProtecteds;
	}
	public void setQuarantinedVaccinatedProtecteds(ArrayList<Agent> quarantinedVaccinatedProtecteds) {
		this.quarantinedVaccinatedProtecteds = quarantinedVaccinatedProtecteds;
	}
	public boolean isIndexHosp() {
		return indexHosp;
	}
	public void setIndexHosp(boolean indexHosp) {
		this.indexHosp = indexHosp;
	}
	public boolean isVaccinationImplemented() {
		return vaccinationImplemented;
	}
	public void setVaccinationImplemented(boolean vaccinationImplemented) {
		this.vaccinationImplemented = vaccinationImplemented;
	}
	public boolean isInfectorInvaded() {
		return infectorInvaded;
	}
	public void setInfectorInvaded(boolean infectiousInvaded) {
		this.infectorInvaded = infectiousInvaded;
	}
}
