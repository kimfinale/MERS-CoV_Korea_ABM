/**
 * Model.java
 * does all the main actions
 */


import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

import org.apache.commons.math3.random.BitsStreamGenerator;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.distribution.AbstractIntegerDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.distribution.UniformIntegerDistribution;
import org.apache.commons.math3.distribution.UniformRealDistribution;
import org.apache.commons.lang3.ArrayUtils;

public class Step {
	// class variables
	static double currentDay = 0.0;

	//////////////////////////////////////////////////////////////////////////////////////////////
	// constructors
	//
	public Step() {
    }


	////////////////////////////////////////////////////////////////////////////////////
	// MERSTransmission()
	// transmission-related steps
	protected void MERSTransmission ( Parameters pars ){
		stepTransmit( pars );
		stepBecomeInfectious( pars );
		stepIsolate( pars );
		stepHospitalShopping( pars );
		stepRemoval( pars );
		stepAge( pars );
		double timeVaccStart = pars.getDayVaccinationStart();
		if( pars.isUnderVaccinationScenario() && timeVaccStart <= currentDay && 0.001 < pars.getVaccCoverage() ) {
			stepVaccinate( pars );
		}
		stepBecomeImmune( pars );
		// do other stuff
		if( pars.isUnderVaccinationScenario() && !pars.isDayVaccinationStartAdjusted()) {
			adjustVaccinationStartDateByCumulativeIncidence( pars );
		}
		if( pars.getDebug() > 0 )
			System.out.printf( "debug = %d, tick = %.1f, MERSTransmission step method done...\n", pars.getDebug(), Step.currentDay );

	}
	
	
	////////////////////////////////////////////////////////////////////////////////
	// stepTransmit( Parameters pars )
	// susceptibles become infected at each hospital where infectious persons exist.
	public void stepTransmit( Parameters pars ){
		double beta = pars.getRateTransmit() * pars.getStepSize();
		double shape = pars.getShapeGammaOffspring();
		for( Hospital hosp : Model.hospitals ) {
//			System.out.println(   "time = " + Step.currentDay + ", hospital id = " + hosp.getID() + ", inf size = "  + hosp.getInfectious().size() ); 
			if( hosp.getInfectious().size() > 0 ) {
				hosp.transmission( pars, beta, shape );
			}
		}
	}

	
	////////////////////////////////////////////////////////////////////////////////
	// becomeInfectious()
	// exposed become infectious
	public void stepBecomeInfectious( Parameters pars ){
		for( Hospital hosp : Model.hospitals ) {
			hosp.infectiousnessDevelopment( pars );
		}
	}
	

	
	
	////////////////////////////////////////////////////////////////////////////////
	// adjustVaccinationStartDateByCumulativeIncidence()
	// adjust the start date of vaccination so that the vaccination can takes place in response to the number of cases occurred
	// and observed
	public void adjustVaccinationStartDateByCumulativeIncidence( Parameters pars ){
//		boolean adjustedAlready = pars.isDayVaccinationStartAdjusted();
//		if( !adjustedAlready ) {
		int vaccTargetCaseNumber = pars.getThresholdNumberCaseForVaccinationInitiation();
		int vaccTargetDate = pars.getThresholdDayVaccinationInitiation();
		int caseUpToNow = Model.getNumPeople( Model.hospitals , "I" );
		caseUpToNow = caseUpToNow + Model.getNumPeople( Model.hospitals , "J" );
		caseUpToNow = caseUpToNow + Model.getNumPeople( Model.hospitals , "R" );
		caseUpToNow = caseUpToNow + Model.getNumPeople( Model.hospitals , "JR" );
		if( vaccTargetCaseNumber <= caseUpToNow  && Step.currentDay <= vaccTargetDate ) {
			pars.setDayVaccinationStart( Step.currentDay );
			pars.setDayVaccinationStartAdjusted( true );
		}

	}
	


	
	/////////////////////////////////////////////////////////////
	// vaccinate( Parameters pars )
	// vaccinate on a hospital basis
	public void stepVaccinate( Parameters pars ){
		double fracVaccTarget = pars.getFracVaccTargetPopulation(); // <1 in case only HCWs or visitors are vaccinated
		double vaccDur = pars.getTimeNeededForVaccination();
		double vaccCoverage = pars.getVaccCoverage();
		double vaccEff = pars.getVaccEfficacy();
		double relVaccEffPostExp = pars.getRelativeVaccEfficacyPostExposure();
		double actualVaccCovSusc = vaccCoverage * fracVaccTarget * vaccEff ; // assuming an all-or-nothing vaccine
		double actualVaccCovExp = vaccCoverage * fracVaccTarget * vaccEff * relVaccEffPostExp;
		double numberOfStepsForVaccination = vaccDur / pars.getStepSize();
		// 1 - (1-vaccProbPerStep)^numberOfStepsForVaccination = actualVaccCov 
		double vaccProbPerStepSusc = 1 - Math.pow( 1-actualVaccCovSusc, 1/numberOfStepsForVaccination );
		double vaccProbPerStepExp = 1 - Math.pow( 1-actualVaccCovExp, 1/numberOfStepsForVaccination );
		// The number of vaccine doses needs to account for those who receive but don't take the vaccine.
		double vaccReceivingProbPerStep = 1 - Math.pow( (1- (vaccCoverage * fracVaccTarget)), 1/numberOfStepsForVaccination );
		
		Model.updateHospitalsForVaccination( pars );
		
		double durVacc = pars.getTimeNeededForVaccination();
		for( Hospital h : Model.hospitalsVaccinationImplemented ) {
			double dayVaccStart = h.getDayVaccinationStarted();
			if( pars.isUnderVaccinationScenario() && dayVaccStart <= currentDay && currentDay < (dayVaccStart + durVacc) ) {
				h.vaccination( pars, vaccProbPerStepSusc, vaccProbPerStepExp, vaccReceivingProbPerStep );
			}
			if( pars.getDebug() > 0 )
				System.out.printf( "tick = %.1f, hospital id = %d, vaccination done.\n", Step.currentDay, h.getID() );
		}

	}
	

/*
	
	/////////////////////////////////////////////////////////////
	// stepVaccinate( Parameters pars )
	// vaccinate on a hospital basis
	public void stepVaccinate( Parameters pars ){
		double fracVaccTarget = pars.getFracVaccTargetPopulation(); // <1 in case only HCWs or visitors are vaccinated
		double vaccDur = pars.getTimeNeededForVaccination();
		double vaccCoverage = pars.getVaccCoverage();
		double vaccEff = pars.getVaccEfficacy();
		double relVaccEffPostExp = pars.getRelativeVaccEfficacyPostExposure();
		
		double actualVaccFracSusc = vaccCoverage * fracVaccTarget * vaccEff ; // assuming an all-or-nothing vaccine
		double actualVaccFracExp = vaccCoverage * fracVaccTarget * vaccEff * relVaccEffPostExp;
				// (1-vaccProb)^dayNeededForVaccination = 1 - vaccCoverage
		double vaccProbSusc = 0.0;
		double vaccProbExp = 0.0;
		vaccProbSusc = 1 - Math.pow( 1-actualVaccFracSusc, pars.getStepSize()/vaccDur );
		vaccProbExp = 1 - Math.pow( 1-actualVaccFracExp, pars.getStepSize()/vaccDur  );
		// adjust the vaccination start date
		
		double vaccStart = pars.getDayVaccinationStart();
		
		int [] vaccTargetRegionID = pars.getVaccTargetRegionID();
		ArrayList<Hospital> hospVacc = new ArrayList<Hospital>();
		hospVacc.addAll( Model.hospitals ); 
		hospVacc.addAll( Model.uninfectedHospitals );
		if( pars.isHospitalTargetedVaccination() ) { // hopsitals to target are determined
			Hospital indexHospital = ((ArrayList<Hospital>) Model.hospitals).get( 0 );
			ArrayList<Hospital> targetHosp = indexHospital.selectHospitalsForVaccination( pars ); 
			for( Hospital h : targetHosp ) {
				if( (currentDay >= vaccStart) && ((currentDay - vaccStart) < vaccDur) ) { 
					h.vaccination( vaccProbSusc, vaccProbExp );
				}
			}
		} else {
			for( Hospital hosp : hospVacc ) {
				if( (currentDay >= vaccStart) && ((currentDay - vaccStart) < vaccDur) ) { 
					if( ArrayUtils.contains( vaccTargetRegionID, hosp.getRegionID() ) ) { // if the hospital is in the region of interest
						hosp.vaccination( vaccProbSusc, vaccProbExp );
					}
				}
			}
		}
		hospVacc.clear();
	}
	

*/	
	
	////////////////////////////////////////////////////////////////////////////////////////
	// becomeImmune( Parameters pars )
	// vaccine recipients become immune 
	public void stepBecomeImmune( Parameters pars ){
		ArrayList<Agent> agentImmune = new ArrayList<Agent>();
		ArrayList<Hospital> hospList = new ArrayList<Hospital>();
		hospList.addAll( Model.hospitals );
		hospList.addAll( Model.uninfectedHospitals );
		for( Hospital hosp : hospList ) {
			ArrayList<Agent> vaccSusc = hosp.getVaccinatedSusceptibles();
			ArrayList<Agent> vaccProtected = hosp.getVaccinatedProtecteds();
			for( Agent a : vaccSusc ) {
				if( a.getDelayVaccineInducedImmunity() < a.getDaySinceVaccination() ) {
					agentImmune.add( a );
					a.setInfectionStatus( "VP" ); // vaccinated and protected
					pars.setCumulVaccProtected( pars.getCumulVaccProtected() + 1 );
				}
			}
			vaccProtected.addAll( agentImmune );
			vaccSusc.removeAll( agentImmune );
			agentImmune.clear();	

			ArrayList<Agent> vaccExp = hosp.getVaccinatedExposeds();
			for( Agent a : vaccExp ) {
				if( a.getDelayVaccineInducedImmunity() < a.getDaySinceVaccination() ) {
					agentImmune.add( a );
					a.setInfectionStatus( "VP" ); // vaccinated, exposed, and protected
					pars.setCumulVaccProtected( pars.getCumulVaccProtected() + 1 );
				}
			}
			vaccProtected.addAll( agentImmune );
			vaccExp.removeAll( agentImmune );
			agentImmune.clear();
			
			ArrayList<Agent> quarantinedVaccSusc = hosp.getQuarantinedVaccinatedSusceptibles();
			ArrayList<Agent> quarantinedVaccProtected = hosp.getQuarantinedVaccinatedProtecteds();
			for( Agent a : quarantinedVaccSusc ) {
				if( a.getDelayVaccineInducedImmunity() < a.getDaySinceVaccination() ) {
					agentImmune.add( a );
					a.setInfectionStatus( "QVP" ); // vaccinated and protected
					pars.setCumulVaccProtected( pars.getCumulVaccProtected() + 1 );
				}
			}
			quarantinedVaccProtected.addAll( agentImmune );
			quarantinedVaccSusc.removeAll( agentImmune );
			agentImmune.clear();	

			ArrayList<Agent> quarantinedVaccExp = hosp.getVaccinatedExposeds();
			for( Agent a : quarantinedVaccExp ) {
				if( a.getDelayVaccineInducedImmunity() < a.getDaySinceVaccination() ) {
					agentImmune.add( a );
					a.setInfectionStatus( "QVP" ); // vaccinated, exposed, and protected
					pars.setCumulVaccProtected( pars.getCumulVaccProtected() + 1 );
				}
			}
			quarantinedVaccProtected.addAll( agentImmune );
			quarantinedVaccExp.removeAll( agentImmune );
			agentImmune.clear();
		}
	}
	
	
	
    ///////////////////////////////////////////////////////////////////////
	public int getNumTransmissions( Agent a, int n, Parameters pars ) {
		double maxNumTransmit = pars.getMaxNumTransmitPerDay();
		if( n > maxNumTransmit ) {
			n = (int) maxNumTransmit;
		}
		return n;
	}	
	
	

	
	////////////////////////////////////////////////////////////////////////////////
	// stepHospitalShopping()
	// infected person move from one hospital to another
	//
	public void stepHospitalShopping( Parameters pars ){
		ArrayList<Hospital> newHospWithInfectious = new ArrayList<>(); // 
		double rateMove = ( 1 / pars.getDayDelayBeforeMovingToAnotherHospital() ) * pars.getStepSize(); 
		ArrayList<ArrayList<Agent>> agentsToMove = new ArrayList< ArrayList<Agent> >();
		ArrayList<ArrayList<Hospital>> hospitalsToMoveTo = new ArrayList<ArrayList<Hospital>>();
		ArrayList<Hospital> hospitalsToMoveFrom = new ArrayList<Hospital>();  
		for( Hospital h : Model.hospitals ) {
			ArrayList<Agent> inf = h.getInfectious();
			ArrayList<Agent> agentsTemp = new ArrayList<Agent>(); // temporary container for agents to move from each hospital
			ArrayList<Hospital> hospitalsTemp = new ArrayList<Hospital>();// temporary container for hospitals to move to from each hospital
			for( Agent a : inf ) {
				if( a.isHighInfectivity() && Model.unifFromZeroToOne.sample() < rateMove ) {
					Hospital h2 = h.selectHospital( pars );
					if( h2 != null ) {
						agentsTemp.add( a );
						hospitalsTemp.add( h2 );
						newHospWithInfectious.add( h2 );
					} 
					else{ 
//						System.out.println( "tick = " + currentDay + ", null returned!" ); 
					}
				}
			}
			if( agentsTemp.size() != hospitalsTemp.size() ) {
				System.err.println( "agents and the number of hospitals to "
						+ "which they are transferred must be the same. Now agents are " 
						+ agentsTemp.size() + " and hospitals are " + hospitalsTemp.size() );
			}
			if( agentsTemp.size() > 0 ) {
				hospitalsToMoveFrom.add( h );
				agentsToMove.add( agentsTemp );
				hospitalsToMoveTo.add( hospitalsTemp );
			}
		}
		// now do the transfer between hospitals
		int size = agentsToMove.size();
		for( int i = 0; i < size; ++ i ) {
			ArrayList<Agent> agentList = agentsToMove.get( i );
			Hospital from = hospitalsToMoveFrom.get( i );
			ArrayList<Hospital> toHospList = hospitalsToMoveTo.get( i );
			int len = agentList.size();
//			System.out.println( "time = " + currentDay + ", i = " + i + ", from id = " + from.getID() + ", from inf before = " + from.getInfectious().size() + ", agents size = " + len );
			for( int j = 0; j < len; ++ j ) {
				Agent a = agentList.get( j ); 
				from.getInfectious().remove( a );
//				System.out.println( "i = " + i + ", j = " + j + ", to id = " + to.get(j).getID() + ", to inf before = " + to.get(j).getInfectious().size() );
				toHospList.get( j ).getInfectious().add( a );
				a.setInvader( true );
//				System.out.println( "i = " + i + ", j = " + j + ", to id = " + to.get(j).getID() + ", to inf = " + to.get(j).getInfectious().size() );
//				System.out.println( "i = " + i + ", j = " + j + ", from id = " + from.getID() + ", from inf = " + from.getInfectious().size() );
			}
		}
		//to remove duplicate using Set
		Set<Hospital> hospSet = new LinkedHashSet<Hospital>( newHospWithInfectious );
		newHospWithInfectious.clear();
		newHospWithInfectious.addAll( hospSet );
		Model.hospitals.addAll( newHospWithInfectious );
		Model.uninfectedHospitals.removeAll( newHospWithInfectious );
		
	}

	
	
/*		
	////////////////////////////////////////////////////////////////////////////////
	// stepHospitalShopping()
	// infected person move from one hospital to another
	//
	public void stepHospitalShopping( Parameters pars ){
		Set<Hospital> newHospitals = new HashSet<>(); // set doesn't allow duplicates, which makes it easy when transferring hospitals form uninfected hospital list to hospital list 
		int numHosp = Model.hospitals.size();
		double prop = pars.getPropSeekingCareFromOtherHospitals() * pars.getStepSize();
		for( int i = 0; i < numHosp; ++i ) {
			Hospital h1 = Model.hospitals.get( i ); 
			ArrayList<Agent> inf1 = h1.getInfectious();
			for (Iterator<Agent> it = inf1.iterator(); it.hasNext(); ) {
			    Agent a = it.next();
				if( Model.unifFromZeroToOne.sample() < prop ) {
					Hospital h2 = h1.selectHospital( pars );
					if( h2 != null ) {
						newHospitals.add( h2 );
						a.setHighInfectivity( true ); //those who move to new hospitals are more infectious 
						h2.getInfectious().add( a );
						it.remove();
					} 
					else{ 
//						System.out.println( "tick = " + currentDay + ", null returned!" ); 
					}
				}
			}
		}
		
		// in case the vaccination started already, the newly selected hospitals need vaccination as well, which was not implemented to increase efficiency 
		for( Hospital h : newHospitals ) {
			int regionID = h.getRegionID();
			int[] targetRegionIDs = pars.getVaccTargetRegionID();
			if( ArrayUtils.contains( targetRegionIDs, regionID ) && currentDay >= pars.getDayVaccinationStart() ) {
				h.delayedVaccination( pars );
			}
			Model.hospitals.add( h );
			Model.uninfectedHospitals.remove( h  );
		}
		
	}

*/	
	
//	////////////////////////////////////////////////////////////////////////////////
//	// move()
//	// infected person move from one hospital to another
//	//
//	public void stepHospitalShopping( Parameters pars ){
//		ArrayList<Agent> agentToMove = new ArrayList<Agent>();
//		ArrayList<Hospital> newHospitals = new ArrayList<Hospital>();
//		int numHosp = Model.hospitals.size();
//		double prop = pars.getPropSeekingCareFromOtherHospitals() * pars.getStepSize();
//		for( int i = 0; i < numHosp; ++i ) {
//			Hospital h1 = Model.hospitals.get( i ); 
//			ArrayList<Agent> inf1 = h1.getInfectious();
//			for( Agent a : inf1 ) {
//				if( Model.unifFromZeroToOne.sample() < prop ) {
//					Hospital h2 = h1.selectHospital( pars );
//					if( h2 != null ) {
//						newHospitals.add( h2 );
//						a.setHighInfectivity( true ); //those who move to new hospitals are more infectious 
//						agentToMove.add( a );
//						int regionID = h2.getRegionID();
//						int[] targetRegionIDs = pars.getVaccTargetRegionID();
//						if( ArrayUtils.contains( targetRegionIDs, regionID ) && currentDay >= pars.getDayVaccinationStart() ) {
//							h2.delayedVaccination( pars );
//						}
//						
//					} 
//					else{ 
////						System.out.println( "tick = " + currentDay + ", null returned!" ); 
//					}
//				}
//			}
//		}
//		
//		
//		// move individuals between hospitals		
//		for( int i = 0; i < agentToMove.size(); ++ i ) {
//			Agent a = agentToMove.get( i );
//			Hospital newHosp = newHospitals.get( i );
//			Hospital oldHosp = a.getHospital();
//			a.setHospital( newHosp );
//			newHosp.getInfectious().add( a );
//			oldHosp.getInfectious().remove( a );
//		}
//		// move hospitals between lists
//		int numNewHosp = newHospitals.size();
//		for( int i = 0; i < numNewHosp; ++ i ) {
//			Hospital h = newHospitals.get(i);
//			Model.hospitals.add( h );
//			Model.uninfectedHospitals.remove( h  );
//		}
//	}

	
	
	////////////////////////////////////////////////////////////////////////////////
	// becomeInfectious()
	// exposed become infectious
	public void stepRemoval( Parameters pars ){
		for( Hospital hosp : Model.hospitals ) {
			hosp.removal( pars );
		}
	}
	
	
	////////////////////////////////////////////////////////////////////////////////
	// stepAge()
	// exposed or infectious people increase their infection days
	// this is used to transfer to the next compartment.
	public void stepAge( Parameters pars ){
		double delta = pars.getStepSize(); 
		ArrayList<Agent> agents = new ArrayList<Agent>();;
		for( Hospital hosp : Model.hospitals ) {
			agents.addAll( hosp.retrieveExposedAgentsFromHospital() );
		}
		for( Hospital hosp : Model.hospitalsVaccinationImplemented ) {
			agents.addAll( hosp.retrieveExposedAgentsFromHospital() );
		}
		for( Agent a : agents ) {
			double daySinceInfect = a.getDaySinceInfection();
			double daySinceVacc = a.getDaySinceVaccination();
			if( daySinceInfect > -100 )
				a.setDaySinceInfection( a.getDaySinceInfection() + delta );
			if( a.getDurationOfIncubation() > -1 )
				a.setDaySinceSymptomOnset( a.getDaySinceSymptomOnset() + delta );
			if( daySinceVacc > -1 )
				a.setDaySinceVaccination( a.getDaySinceVaccination() + delta );
		}
	}

	
	/////////////////////////////////////////////////////////////////////////////////
	// stepIsolate( Parameter pars )
	// agents get isolated some day after becoming symptomatic: I -> J
	protected void stepIsolate( Parameters pars  ){
		double[] meanDelay = pars.getMeanTimeToIsolation();
		double[] maxDelay = pars.getMaxTimeToIsolation();
		int[] dayCutoff = pars.getPeriodCutOff();
		
		double mean = meanDelay[ 0 ];
		double max = maxDelay[ 0 ];
		if( dayCutoff[0] <= currentDay && currentDay < dayCutoff[1] ) {
			mean = meanDelay[ 1 ];
			max = maxDelay[ 1 ];
		} 
		else if( dayCutoff[1] <= currentDay ){
			mean = meanDelay[ 2 ];
			max = maxDelay[ 2 ];
		}
		double rate = (1 / mean) * pars.getStepSize();
		for( Hospital hosp : Model.hospitals ) {
			int numAlreadyIsolated = hosp.getIsolateds().size(); // record the number of people who were already isolated
			hosp.isolation( pars, rate, max );
			if( (hosp.getIsolateds().size() - numAlreadyIsolated) > 0 ) { // if anybody were newly isolated
				hosp.quarantine( pars, rate, max );
			}
		}
	}
	
}
