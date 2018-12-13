/**
 * Model.java
 * does all the main actions
 */

import java.util.*;


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
		if( pars.getOutbreakScenario().equalsIgnoreCase("2015")  ) {
			adjustIndexHospIndexCase( pars );
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
		ArrayList<Hospital> targetHospitals = Model.hospitals;
//		targetHospitals.removeAll( Model.hospitalsCaseIsolated );  // transmissions do not occur once cases are isolated
		for( Hospital hosp : targetHospitals ) {
			if( hosp.getInfectious().size() > 0 ) {
				hosp.transmission( pars, beta, shape );
			}
		}
		if( pars.getDebug() > 0 ) {
			System.out.printf( "Day: %.1f, Step.stepTransmit occurred.\n", Step.currentDay );
			for(Hospital h : Model.hospitals ) {
				h.printSelf();
			}
		}
	
	}

	
	////////////////////////////////////////////////////////////////////////////////
	// becomeInfectious()
	// exposed become infectious
	public void stepBecomeInfectious( Parameters pars ){
		for( Hospital hosp : Model.hospitals ) {
			hosp.becomingInfectious( pars );
		}
		if( pars.getDebug() > 0 ) {
			System.out.printf( "Day: %.1f, Step.stepBecomeInfectious occurred.\n", Step.currentDay );
			for(Hospital h : Model.hospitals ) {
				h.printSelf();
			}
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
		int isolatedCaseUpToNow = Model.getNumPeople( Model.hospitals , "J" ) + Model.getNumPeople( Model.hospitals , "JR" );	
		if( vaccTargetCaseNumber <= isolatedCaseUpToNow  && Step.currentDay <= vaccTargetDate ) {
			pars.setDayVaccinationStart( Step.currentDay );
			pars.setDayVaccinationStartAdjusted( true );
		}

	}
	
	
	
	////////////////////////////////////////////////////////////////////////////////
	// adjustIndexHospIndexCase()
	// adjust the start date of vaccination so that the vaccination can takes place in response to the number of cases occurred
	// and observed
	public void adjustIndexHospIndexCase( Parameters pars ){
		double observedTimeToIndexCase = pars.getTimeIndexCaseConfirmation();
		if( currentDay <= observedTimeToIndexCase && 
				observedTimeToIndexCase < currentDay + pars.getStepSize() ) {
			for( Hospital h : Model.hospitals ) {
				if( h.isIndexHosp()) {
					if( ! Model.hospitalsCaseIsolated.contains( h ) ) {
						Model.hospitalsCaseIsolated.add( h );
					}
					ArrayList<Agent> list = new ArrayList<Agent>();
					for( Agent a : h.getRemoveds() ) {
						if( a.isIndexCase() ) {
							a.setInfectionStatus( "JR" );
							list.add( a );
						}
						break;
					}
					h.getRemoveds().removeAll( list );
					h.getIsolatedRemoveds().addAll( list );
					break;	
				}
			}
			if( pars.getDebug() > 3 )
				System.out.println( "index case method done.. tick = " + Step.currentDay );
		}
		
	}
	
	
	/////////////////////////////////////////////////////////////
	// vaccinate( Parameters pars )
	// vaccinate on a hospital basis
	public void stepVaccinate( Parameters pars ){
		Model.updateHospitalsForVaccination( pars ); // update the list of hospitals 
		
		double vaccDur = pars.getTimeNeededForVaccination();
		double vaccCoverage = pars.getVaccCoverage();
		double vaccEff = pars.getVaccEfficacy();
		double relVaccEffPostExp = pars.getRelativeVaccEfficacyPostExposure();
		double actualVaccCovSusc = vaccCoverage  * vaccEff ; // assuming an all-or-nothing vaccine
		double actualVaccCovExp = vaccCoverage * vaccEff * relVaccEffPostExp;
		double numberOfStepsForVaccination = vaccDur / pars.getStepSize();
		// 1 - (1-vaccProbPerStep)^numberOfStepsForVaccination = actualVaccCov 
		double vaccProbPerStepSusc = 1 - Math.pow( 1-actualVaccCovSusc, 1/numberOfStepsForVaccination );
		double vaccProbPerStepExp = 1 - Math.pow( 1-actualVaccCovExp, 1/numberOfStepsForVaccination );
		// The number of vaccine doses needs to account for those who receive but don't take the vaccine.
		double vaccReceivingProbPerStep = 1 - Math.pow( (1- (vaccCoverage )), 1/numberOfStepsForVaccination );
		
		double durVacc = pars.getTimeNeededForVaccination();
		for( Hospital h : Model.hospitalsVaccinationImplemented ) {
			double dayVaccStart = h.getDayVaccinationStarted();
			if( dayVaccStart <= currentDay && currentDay < (dayVaccStart + durVacc) ) {
				h.vaccination( pars, vaccProbPerStepSusc, vaccProbPerStepExp, vaccReceivingProbPerStep );
			}
			if( pars.getDebug() > 0 )
				System.out.printf( "tick = %.1f, hospital id = %d, vaccination done.\n", Step.currentDay, h.getID() );
		}

		System.out.printf( "Day: %.1f, Step.stepVaccinate occurred.\n", Step.currentDay );
		for(Hospital h : Model.hospitals ) {
			h.printSelf();
		}
	}
	
	
	////////////////////////////////////////////////////////////////////////////////////////
	// becomeImmune( Parameters pars )
	// vaccine recipients become immune 
	public void stepBecomeImmune( Parameters pars ){
		ArrayList<Agent> agentImmune = new ArrayList<Agent>();
		ArrayList<Hospital> hospList = new ArrayList<Hospital>();
		hospList.addAll( Model.hospitals );
		hospList.addAll( Model.uninfectedHospitals );
		for( Hospital hosp : hospList ) {
			ArrayList<Agent> VS = hosp.getVaccinatedSusceptibles();
			for( Agent a : VS ) {
				if( a.getDelayVaccineInducedImmunity() < a.getDaySinceVaccination() ) {
					agentImmune.add( a );
					a.setInfectionStatus( "VP" ); // vaccinated and protected
					pars.setCumulVaccProtected( pars.getCumulVaccProtected() + 1 );
				}
			}
			hosp.getVaccinatedProtecteds().addAll( agentImmune );
			VS.removeAll( agentImmune );
			agentImmune.clear();	

			ArrayList<Agent> VE = hosp.getVaccinatedExposeds();
			for( Agent a : VE ) {
				if( a.getDelayVaccineInducedImmunity() < a.getDaySinceVaccination() ) {
					agentImmune.add( a );
					a.setInfectionStatus( "VP" ); // vaccinated, exposed, and protected
					pars.setCumulVaccProtected( pars.getCumulVaccProtected() + 1 );
				}
			}
			hosp.getVaccinatedProtecteds().addAll( agentImmune );
			VE.removeAll( agentImmune );
			agentImmune.clear();
			
			ArrayList<Agent> QVS = hosp.getQuarantinedVaccinatedSusceptibles();
			for( Agent a : QVS ) {
				if( a.getDelayVaccineInducedImmunity() < a.getDaySinceVaccination() ) {
					agentImmune.add( a );
					a.setInfectionStatus( "QVP" ); // vaccinated and protected
					pars.setCumulVaccProtected( pars.getCumulVaccProtected() + 1 );
				}
			}
			hosp.getQuarantinedVaccinatedProtecteds().addAll( agentImmune );
			QVS.removeAll( agentImmune );
			agentImmune.clear();	

			ArrayList<Agent> QVE = hosp.getQuarantinedVaccinatedExposeds();
			for( Agent a : QVE ) {
				if( a.getDelayVaccineInducedImmunity() < a.getDaySinceVaccination() ) {
					agentImmune.add( a );
					a.setInfectionStatus( "QVP" ); // vaccinated, exposed, and protected
					pars.setCumulVaccProtected( pars.getCumulVaccProtected() + 1 );
				}
			}
			hosp.getQuarantinedVaccinatedProtecteds().addAll( agentImmune );
			QVE.removeAll( agentImmune );
		}
		if( pars.getDebug() > 0 ) {
			System.out.printf( "Day: %.1f, Step.stepBecomeImmune occurred.\n", Step.currentDay );
			for( Hospital h : Model.hospitals ) {
				h.printSelf();
			}
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
				System.err.println( "Step.stepHospitalShopping: agents and the number of hospitals to "
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
		for( Hospital h : newHospWithInfectious ) {
			if( !Model.hospitals.contains( h ) ) {
				Model.hospitals.add( h );
			}
		}
		Model.uninfectedHospitals.removeAll( newHospWithInfectious );
		if( pars.getDebug() > 0 ) {
			System.out.printf( "Day: %.1f, Step.stepHospitalShopping occurred.\n", Step.currentDay );
			for(Hospital h : Model.hospitals ) {
				h.printSelf();
			}
		}
	}


	
	////////////////////////////////////////////////////////////////////////////////
	// becomeInfectious()
	// exposed become infectious
	public void stepRemoval( Parameters pars ){
		for( Hospital hosp : Model.hospitals ) {
			hosp.removal( pars );
		}
		if( pars.getDebug() > 0 ) {
			System.out.printf( "Day: %.1f, Step.stepRemoval occurred.\n", Step.currentDay );
			for(Hospital h : Model.hospitals ) {
				h.printSelf();
			}
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
		if( pars.getDebug() > 0 ) {
			System.out.printf( "Day: %.1f, Step.stepAge occurred.\n", Step.currentDay );
			for(Hospital h : Model.hospitals ) {
				h.printSelf();
			}
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
		if( pars.getDebug() > 0 ) {
			System.out.printf( "Day: %.1f, Step.stepIsolate occurred.\n", Step.currentDay );
			for(Hospital h : Model.hospitals ) {
				h.printSelf();
			}
		}
	}
	
}
