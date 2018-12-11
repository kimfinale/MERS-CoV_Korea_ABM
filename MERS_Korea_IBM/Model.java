
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.Writer;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.Random;
import java.util.Scanner;
import java.util.Set;

//import com.sun.xml.internal.ws.policy.privateutil.PolicyUtils.Collections;

//import org.apache.commons.math3.util.FastMath;

//import org.apache.commons.math3.random.MersenneTwister;
//import org.apache.commons.math3.random.RandomDataGenerator;
//import org.apache.commons.math3.util.FastMath;

//import jonghoonkim.methods.ode.RungeKutta;
//import cern.jet.random.*; 
//import cern.jet.random.tdouble.*;
//import cern.jet.random.tdouble.engine.DoubleMersenneTwister;
//import cern.jet.random.tdouble.engine.DoubleRandomEngine;

import org.apache.commons.math3.random.BitsStreamGenerator;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.distribution.AbstractIntegerDistribution;
import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.distribution.UniformIntegerDistribution;
import org.apache.commons.math3.distribution.UniformRealDistribution;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;



public class Model {

	static ArrayList<Hospital> hospitals = new ArrayList<Hospital>();
	static ArrayList<Hospital> uninfectedHospitals = new ArrayList<Hospital> ();
	// the above two hospitals give the total hospitals if summed
	static Set<Hospital> hospitalsTransmissionOccurred = new HashSet<Hospital>();
	// to help implement vaccination
	static ArrayList<Hospital> hospitalsCaseIsolated = new ArrayList<Hospital>();
	static ArrayList<Hospital> hospitalsVaccinationImplemented = new ArrayList<Hospital>();

	static Parameters pars = new Parameters ();
	static Hospital hosp = new Hospital ();
	static Utility util = new Utility ();
	
	static ArrayList<Integer>  	hospLevel 		= pars.readFile( pars.getFilePathHospLevel() );
	static ArrayList<Integer>  	hospRegionID 	= pars.readFile( pars.getFilePathHospRegion() );
	static ArrayList<Double>  	hospPopAtRisk 	= pars.readFileDouble( pars.getFilePathHospSize() );
	static ArrayList<Double>  	hospLongitude 	= pars.readFileDouble( pars.getFilePathHospLongitude() );
	static ArrayList<Double>  	hospLatitude 	= pars.readFileDouble( pars.getFilePathHospLatitude() );
	// random number generators
	static BitsStreamGenerator RNG = new MersenneTwister();
	static AbstractRealDistribution unifFromZeroToOne = new UniformRealDistribution( RNG, 0, 1 );
	static AbstractRealDistribution gammaIncubationPeriod = 
			new GammaDistribution( RNG, pars.getShapeGammaDurationOfIncubation(), 1/pars.getRateGammaDurationOfIncubation() );
	static AbstractRealDistribution gammaDurationOfInfectiousness = 
			new GammaDistribution( RNG, pars.getShapeGammaDurationOfInfectiousness(), 1/pars.getRateGammaDurationOfInfectiousness() );
	static AbstractRealDistribution gammaDelayVaccineInducedImmunity = 
			new GammaDistribution( RNG, pars.getShapeGammaDurationOfIncubation(), pars.getMeanDelayVaccineInducedImmunity() / pars.getShapeGammaDurationOfIncubation() );
	static AbstractIntegerDistribution unifHosp = new UniformIntegerDistribution( RNG, 0, hospPopAtRisk.size() - 1 );
	static AbstractIntegerDistribution uniformDelay = 
			new UniformIntegerDistribution( RNG, pars.getMinDelaySymptomOnsetToTransmission2015(), pars.getMaxDelaySymptomOnsetToTransmission2015() );
	static AbstractRealDistribution unifFromZeroToOneVacc = new UniformRealDistribution( RNG, 0, 1 );
	
	
	public static void main( String[] args ) throws IOException {
		
		DescriptiveStatistics cumCase = new DescriptiveStatistics();
		DescriptiveStatistics offspringVariance = new DescriptiveStatistics();
		DescriptiveStatistics hospitalsWithInfectious = new DescriptiveStatistics();
		DescriptiveStatistics hospitalsTransmissionOccurred = new DescriptiveStatistics();
		DescriptiveStatistics cumulVaccDose = new DescriptiveStatistics();
		DescriptiveStatistics cumulVaccDoseSet = new DescriptiveStatistics();
		

		for( int i = 0; i < pars.getNumSampleRun(); ++i ){
			System.out.println( "i = "  + i );			
			int numSteps= (int) ( pars.getStopTime() / pars.getReportFreq() );
//			Parameters pars = new Parameters ();
//			pars.setRandomSeed( i );

			double [][] out = runModel( pars );
			
			
			cumCase.addValue( out[numSteps-1][6] );
			offspringVariance.addValue( out[numSteps-1][7] );
			hospitalsWithInfectious.addValue( Model.hospitals.size() );
			hospitalsTransmissionOccurred.addValue( out[numSteps-1][8]  );
			cumulVaccDose.addValue( pars.getCumulVaccDose()  );
			cumulVaccDoseSet.addValue( getTotalVaccineDosesUsed() );

		}
		double mean = cumCase.getMean();
		double sd = cumCase.getStandardDeviation();
		double lower = cumCase.getPercentile( 2.5 );
		double med = cumCase.getPercentile( 50 );
		double upper = cumCase.getPercentile( 97.5 );
		System.out.println( "cumul case mean = " + mean + ", sd = " + sd  + ", lower = " + lower  + ", med = " + med  + ", upper = " + upper );
		double meanOffspringVariance = offspringVariance.getMean();
		double sdOffspringVariance = offspringVariance.getStandardDeviation();
		double lowerOffspringVariance = offspringVariance.getPercentile( 2.5 );
		double medOffspringVariance = offspringVariance.getPercentile( 50 );
		double upperOffspringVariance = offspringVariance.getPercentile( 97.5 );
		System.out.println( "offspring var mean = " + meanOffspringVariance + ", sd = " + sdOffspringVariance  + ", lower = " + lowerOffspringVariance  
				+ ", med = " + medOffspringVariance  + ", upper = " + upperOffspringVariance );
		System.out.println( "hospitalsVisited = " + hospitalsWithInfectious.getMean() + ", sd = " + + hospitalsWithInfectious.getStandardDeviation() );
		System.out.println( "hospitalsTransmissionOccurred = " + hospitalsTransmissionOccurred.getMean() + ", sd= " + + hospitalsTransmissionOccurred.getStandardDeviation() );
		System.out.println( "cumulVaccDose = " + cumulVaccDose.getMean() + ", sd = " + + cumulVaccDose.getStandardDeviation() );
//		System.out.println( "cumulVaccProtected = " + pars.getCumulVaccProtected() );
		System.out.println( "cumulVaccDoseSet = " + cumulVaccDoseSet.getMean() + ", sd = " + + cumulVaccDoseSet.getStandardDeviation() );
	
	}



	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// runModel()
	// run the model 
	public static double [][] runModel( Parameters pars ){
		renewRandomNumberGenerators( pars );
		modelSetup( pars );
		double delta = pars.getStepSize();
		int nIter = (int) ( pars.getStopTime() / delta );
		Step step = new Step();
		int totalStepsReported = (int) ( pars.getStopTime() / pars.getReportFreq() );
		double [][] out = new double[ totalStepsReported ][ 9 ];// day, s, e, i, j, r, ci, varToMeanRatio, numHospitalVisited 				
		for( int i = 0; i < totalStepsReported; ++ i ) {
			out[ i ][ 0 ] = i; 
		}
		int numSusc = getNumPeople( hospitals, "S" );
		int numExp = getNumPeople( hospitals, "E" );
		int numInf = getNumPeople( hospitals, "I" );
		int numIso = getNumPeople( hospitals, "J" );
		int numRem = getNumPeople( hospitals, "R" );
		
		out[ 0 ][ 1 ] = numSusc ;//pars.getNumInitPop() - exp - inf - iso; //s, susceptible
		out[ 0 ][ 2 ] = numExp; //e, exposed
		out[ 0 ][ 3 ] = numInf; //i, infectious
		out[ 0 ][ 4 ] = numIso; //j, isolated
		out[ 0 ][ 5 ] = numRem; //r, removed
		out[ 0 ][ 6 ] = pars.getCumulInc(); //cumulative number of infections
		out[ 0 ][ 7 ] = 0.0; //variance to mean ratio in the number of offspring
		out[ 0 ][ 8 ] = 1.0; //number of hospitals infected
		int index = 1; // initial values of out at index=0 was defined previously
		int everyNthStep = (int) ( pars.getReportFreq() / delta );
		for( int i = 0; i < nIter; ++ i ){
			step.MERSTransmission( pars );
			Step.currentDay = Step.currentDay + delta;
			if( i > 0 && i % everyNthStep == 0 ){
				out[ index ][ 1 ] = getNumPeople( hospitals, "S" );
				out[ index ][ 2 ] = getNumPeople( hospitals, "E" );
				out[ index ][ 3 ] = getNumPeople( hospitals, "I" );
				out[ index ][ 4 ] = getNumPeople( hospitals, "J" ); 
				out[ index ][ 5 ] = getNumPeople( hospitals, "R" ); // protected by vaccine
				out[ index ][ 6 ] = pars.getCumulInc();
				out[ index ][ 7 ] = offspringVarianceToMeanRatio( pars );
//				out[ index ][ 8 ] = getNumberHospitalInfected( pars );
//				out[ index ][ 8 ] = new HashSet<Hospital>( hospitalsTransmissionOccurred ).size(); // to remove duplicates
				out[ index ][ 8 ] = hospitalsTransmissionOccurred.size(); // to remove duplicates
				index ++;
				if( pars.getDebug() > 2 )
					System.out.println( "total pop = " + getTotalPopulationSize() +  ", cumul inc = " + pars.getCumulInc() );
			}

		}
		return out;
	}
	
	
	
	
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// modelSetup()
	// populate agents according to the number of initial population size and initially infecteds
	// assumes that the index case is high-risk, which is the case for 2015 outbreak
	public static void modelSetup( Parameters pars ){
		// without setting the arraylists that were created to null, the updated distribution at the end of the previous simulation
		// become the initial distribution when running from rJava.
		hospitals = new ArrayList<Hospital>();
		uninfectedHospitals = new ArrayList<Hospital> ();
		hospitalsTransmissionOccurred = new HashSet<Hospital>();
		hospitalsCaseIsolated = new ArrayList<Hospital>();
		hospitalsVaccinationImplemented = new ArrayList<Hospital>();
		// reset variables for vaccination to occur properly
		pars.setDayVaccinationStartAdjusted( false ); // for the vaccination start date to be adjusted 
		pars.setDayVaccinationStart( pars.getStopTime() ); // reset to the stop time
		pars.setCumulVaccDose( 0 );
		pars.setCumulVaccProtected( 0 );
		
		Agent.nextID = 0;
		Hospital.nextID = 0;
		Step.currentDay = 0;
		pars.setCumulInc( 0 );
		
		// generate hospitals 
		int numHosp = hospLatitude.size();	
		for( int i = 0; i < numHosp; ++ i ) {
			if( hospLevel.get( i ) > pars.getCutoffHospitalLevel() ) {
				int pop = (int) ((double) hospPopAtRisk.get(i) );
				Hospital h = new Hospital( pop );
				h.setLatitude( hospLatitude.get( i ) );
				h.setLongitude( hospLongitude.get( i ) );
				h.setRegionID( hospRegionID.get( i ) );
				h.setLevel( hospLevel.get( i ) );
				uninfectedHospitals.add( h );
			}
		}
		
		String scenario = pars.getOutbreakScenario();
		
		if( scenario.equalsIgnoreCase("2015") ) { // Pyeongtaek St. Mary's hospital, index case being high-risk, successful transmission 4-6 days after symptom onset
//			double lon = 127.074340007;
//			double lat = 37.0084094413;
			int index = 0; 
			for( int i = 0; i < uninfectedHospitals.size(); ++ i  ) {
				double hospLon = uninfectedHospitals.get(i).getLongitude();
				if( 127.07434 < hospLon && hospLon < 127.074341 ) {
					index = i;
					break;
				}
			}
			
			Hospital indexHosp = uninfectedHospitals.get( index  ); // 
			uninfectedHospitals.remove( indexHosp );
			hospitals.add( indexHosp );	// hospitals where infected peoples exist are separately tracked in the list hospitals
			hospitalsTransmissionOccurred.add( indexHosp );
			
			int offspringSize = pars.getFirstGenerationOffspring2015();
			
			ArrayList<Agent> susc = indexHosp.getSusceptibles();
			ArrayList<Agent> isolatedRemoveds = indexHosp.getIsolatedRemoveds();
			
			Agent indexCase = susc.get( 0 );
			isolatedRemoveds.add( indexCase );
			indexCase.setInfectionStatus( "R" ); // I
			indexCase.setNumOffspring( offspringSize );
			indexCase.setHospital( indexHosp );
			indexCase.setGeneration( 0 );
			susc.remove( indexCase );
			pars.setCumulInc( 1 );
			
			ArrayList<Agent> newlyExposed = new ArrayList<Agent>();
			
			for( int i = 0; i < offspringSize; ++i ){
				Agent a = susc.get( i );
				a.setInfectorID( indexCase.getID() );
				a.setGeneration( indexCase.getGeneration() + 1 );
				int delay = uniformDelay.sample();
				a.setDaySinceInfection( - 1 * delay ); //
				a.gammaDurationOfIncubation( pars );
				a.setDurationOfIncubation( a.getDurationOfIncubation() + delay );
				a.setInfectionStatus( "E" );
				if( unifFromZeroToOne.sample() < pars.getPropSeekingCareFromOtherHospitals() ) {
					a.setHighInfectivity( true );
				}
				newlyExposed.add( a );
			}
			susc.removeAll( newlyExposed );
			indexHosp.getExposeds().addAll( newlyExposed );
		}
		else if( scenario.equalsIgnoreCase( "Importation") ) {
			// should I include all hospitals instead of Level 3 or 4 hospitals?
			Hospital indexHosp = uninfectedHospitals.get( unifHosp.sample() ); // the index case appears in a random hospital
			uninfectedHospitals.remove( indexHosp );
			hospitals.add( indexHosp );	// hospitals where infected peoples exist are separately tracked in the list hospitals
			
			ArrayList<Agent> susc = indexHosp.getSusceptibles();
			
			int initInf = pars.getNumberInfectedsImported();
			ArrayList<Agent> agentsInfectious = new ArrayList<Agent>();
			for( int i = 0; i < initInf; ++i ){
				Agent a = susc.get( i );
				a.becomeInfectious( pars );
				double r = unifFromZeroToOne.sample();
				if( r < pars.getPropSeekingCareFromOtherHospitals() ) {
					a.setHighInfectivity( true ); //the index case is a super-spreader
				}
				agentsInfectious.add( a );
			}
			susc.removeAll( agentsInfectious );
			indexHosp.getInfectious().addAll( agentsInfectious );
		}
		else {
			AbstractIntegerDistribution unifHospIndexCase = new UniformIntegerDistribution( RNG, 0, uninfectedHospitals.size() - 1 );
			Hospital indexHosp = uninfectedHospitals.get( unifHospIndexCase.sample() );
			uninfectedHospitals.remove( indexHosp );
			hospitals.add( indexHosp );	// hospitals where infected peoples exist are separately tracked in the list hospitals
			int initNumInf = pars.getNumInitInfectious();
			ArrayList<Agent> susc = indexHosp.getSusceptibles();
			ArrayList<Agent> agentsInfectious = new ArrayList<Agent>();
			for( int i = 0; i < initNumInf; ++i ){
				Agent a = susc.get( i );
				a.becomeInfectious( pars );
				a.setHighInfectivity( true ); //the index case is a super-spreader
				a.setInvader( true );
				agentsInfectious.add( a );
			}
			susc.removeAll( agentsInfectious );
			indexHosp.getInfectious().addAll( agentsInfectious );
		}	
		
	}	

	
	

	// These methods are created to return one-dimensional longitude and latitude array separately, 
	// to make it easier to call and treat in R using rJava
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// getAffectedHospitalCoordinateArray()
	// retrieve the array of longitude of the hospitals where transmission took place
	public static double [] getAffectedHospitalLongitudeArray(){
		int n = hospitalsTransmissionOccurred.size();
		double[] lonArray = new double[ n ];
		int index = 0;
		for( Hospital h : hospitalsTransmissionOccurred ){
			lonArray[ index ] = h.getLongitude();
			index ++;
		}
		return lonArray;
	}
	
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// getAffectedHospitalLatitudeArray()
	// retrieve the array of latitude of the hospitals where transmission took place
	public static double [] getAffectedHospitalLatitudeArray(){
		int n = hospitalsTransmissionOccurred.size();
		double[] latArray = new double[ n ];
		int index = 0;
		for( Hospital h : hospitalsTransmissionOccurred ){
			latArray[ index ] = h.getLatitude();
			index ++;
		}
		return latArray;
	}


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// getAffectedHospitalArray()
	// retrieve the array of IDs of the hospitals where transmission took place
	public static int [] getAffectedHospitalIDArray(){
		int n = hospitalsTransmissionOccurred.size();
		int[] IDArray = new int[ n ];
		int index = 0;
		for( Hospital h : hospitalsTransmissionOccurred ){
			IDArray[ index ] = h.getID();
		}
		return IDArray;
	}
		
	
	
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// getTransmissionGenerationArray()
	// return a list of infectious (I or J) individuals that have generated secondary transmissions
	public static int[] getTransmissionGenerationArray(){
		ArrayList<Agent> potentialInfectors = getPotentialInfectors ( Model.hospitals );
		ArrayList<Integer> offspring = new ArrayList<Integer>();
		for( int i = 0; i < potentialInfectors.size(); ++ i ) {
			offspring.add( potentialInfectors.get( i ).getNumOffspring()  );
		}
		int[] arr = offspring.stream().mapToInt(i -> i).toArray();
		return arr;	
	}
	
	
	
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// getNumOffspringArray()
	// return a list of infectious (I or J) individuals that have generated secondary transmissions
	public static int[] getNumOffspringArray(){
		ArrayList<Agent> potentialInfectors = getPotentialInfectors ( Model.hospitals );
		ArrayList<Integer> offspring = new ArrayList<Integer>();
		for( int i = 0; i < potentialInfectors.size(); ++ i ) {
			offspring.add( potentialInfectors.get(i).getNumOffspring()  );
		}
		int[] arr = offspring.stream().mapToInt(i -> i).toArray();
		return arr;	
	}
	

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// numOffspringList( Parameters pars )
	// return a list of infectious (I or J) individuals that could have generated secondary transmissions
	public static ArrayList<Integer> numOffspringList( Parameters pars ){
		ArrayList<Agent> potentialInfectors = getPotentialInfectors ( Model.hospitals );
		ArrayList<Integer> offspring = new ArrayList<Integer>();
		for( Agent a : potentialInfectors ) {
			offspring.add( a.getNumOffspring() );
		}
		if( pars.getDebug() > 3 )
			System.out.printf( "tick = %.1f, num potential infector = %d\n", Step.currentDay, potentialInfectors.size() );
		
		return offspring;	
	}
	
	
	
	

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// getPotentialInfectors( ArrayList<Hospital> list)
	// return a list of agents who are potential infectors (i.e., I, J, R, JR )
	public static ArrayList<Agent> getPotentialInfectors( ArrayList<Hospital> list ){
		ArrayList<Agent> potentialInfectors = new ArrayList<Agent>();
		for( Hospital hosp : list ) {
			potentialInfectors.addAll( hosp.getInfectious() );
			potentialInfectors.addAll( hosp.getIsolateds() );
			potentialInfectors.addAll( hosp.getRemoveds() );
			potentialInfectors.addAll( hosp.getIsolatedRemoveds() );
//			potentialInfectors.addAll( hosp.getExposeds() );
//			potentialInfectors.addAll( hosp.getQuarantinedExposeds() );
//			potentialInfectors.addAll( hosp.getVaccinatedExposeds() );
		}
		
		return potentialInfectors;	
	}
	

	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// offspringVarianceToMeanRatio()
	// variance to mean regarding the offspring size
	public static double offspringVarianceToMeanRatio( Parameters pars ){
		ArrayList<Integer> offspringDist = numOffspringList( pars );
		double varToMeanRatio = 0;
		double len = offspringDist.size();
		double sum = 0;
		double mean = 0;
		double var = 0;
		double ssq = 0;
		if( len > 0 ) {  // 
			for( int i = 0; i < len; ++ i ) {
				sum += offspringDist.get(i);
			}
			mean = sum / len;
			for( int i = 0; i < len; ++ i ) {
				ssq += (offspringDist.get(i) - mean ) * (offspringDist.get(i) - mean ) ;
			}
			var = ssq / (len-1);
		}
		if( mean > 0 ) 
			varToMeanRatio = var / mean; 
		if( pars.getDebug() > 2 )
			System.out.printf( "tick = %.1f, offspring size mean = %.2f, variance = %.2f\n", Step.currentDay, mean, var );
		return varToMeanRatio;
	}


	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// renewRandomNumberGenerators()
	// renew the random number generators such that we get the same result in case we don't change the random seed in R (via rJava) and Java
	//
	public static void renewRandomNumberGenerators( Parameters pars ){
		RNG = new MersenneTwister( pars.getRandomSeed() );
		unifFromZeroToOne = new UniformRealDistribution( RNG, 0, 1 );
		gammaIncubationPeriod = 
				new GammaDistribution( RNG, pars.getShapeGammaDurationOfIncubation(), 1/pars.getRateGammaDurationOfIncubation() );
		gammaDurationOfInfectiousness = 
				new GammaDistribution( RNG, pars.getShapeGammaDurationOfInfectiousness(), 1/pars.getRateGammaDurationOfInfectiousness() );
		gammaDelayVaccineInducedImmunity = 
				new GammaDistribution( RNG, pars.getShapeGammaDurationOfIncubation(), pars.getMeanDelayVaccineInducedImmunity() / pars.getShapeGammaDurationOfIncubation() );
		unifHosp = new UniformIntegerDistribution( RNG, 0, hospPopAtRisk.size() - 1 );
		uniformDelay = 
				new UniformIntegerDistribution( RNG, pars.getMinDelaySymptomOnsetToTransmission2015(), pars.getMaxDelaySymptomOnsetToTransmission2015() );
		unifFromZeroToOneVacc = new UniformRealDistribution( RNG, 0, 1 );
		
	}

	
	/////////////////////////////////////////////////////////////////////////////////////////////
	// getNumPeople( ArrayList<Hospital> list, String str )
	// return the number of people from a list of hospitals in a state specified by the string 
	public static int getNumPeople( ArrayList<Hospital> list, String str ){
		int num = 0; 
		if( str.equals( "S" ) ) {
			for( Hospital h : list ) {
				num += h.getSusceptibles().size();
				num += h.getVaccinatedSusceptibles().size();
				num += h.getQuarantinedSusceptibles().size();
			}
		} else if( str.equals( "E" ) ) {
			for( Hospital h : list ) {
				num += h.getExposeds().size();
				num += h.getVaccinatedExposeds().size();
				num += h.getQuarantinedExposeds().size();
			}
		} else if( str.equals( "I" ) ) {
			for( Hospital h : list ) {
				num += h.getInfectious().size();
			}
		} else if( str.equals( "J" ) ) {
			for( Hospital h : list ) {
				num += h.getIsolateds().size();
			}
		}
		else if( str.equals( "R" ) ) {
			for( Hospital h : list ) {
				num += h.getRemoveds().size();
				num += h.getIsolatedRemoveds().size();
				num += h.getVaccinatedProtecteds().size();
			}
		}

		return num;
	}

	
	

	/////////////////////////////////////////////////////////////////////////////////////////////
	// getNumPeople( ArrayList<Hospital> list, String str )
	// return the number of people from a list of hospitals in a state specified by the string 
	public static int getTotalPopulationSize(){
		int totalPop = 0; 
		for ( Hospital h : Model.hospitals ) {
			totalPop += h.getPopulationSize();
		}
		for ( Hospital h : Model.uninfectedHospitals ) {
			totalPop += h.getPopulationSize();
		}
		
		return totalPop;
	}

	
	
	
	/////////////////////////////////////////////////////////////////////////////////////////////
	// getDayVaccinationStartAdjusted( Parameters pars )
	// return the number of people from a list of hospitals in a state specified by the string 
	public static int getDayVaccinationStartAdjusted( Parameters pars ){
		int status = 0; 
		if( pars.isDayVaccinationStartAdjusted() ) {
			status = 1;
		}
		return status;
	}
	
	
	/////////////////////////////////////////////////////////////////////////////////////////////
	// getNumPeople( ArrayList<Hospital> list, String str )
	// return the number of people from a list of hospitals in a state specified by the string 
	public static int getTotalVaccineDosesUsed(){
		int num = 0; 
		for( Hospital h : hospitalsVaccinationImplemented ) {
			num += h.getVaccineReceived().size();
		}
		return num;
	}
	
	
	
	/////////////////////////////////////////////////////////////////////////////////////////////
	// getTotalPopulationSize( ArrayList<Hospital> list )
	// return the number of people from hospital and sum them up 
	public static int getTotalPopulationSize( ArrayList<Hospital> list ){
		int num = 0; 
		for( Hospital h : list ) {
			num += h.getPopulationSize();					
		}
		return num;
	}

	
	
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// testRNG( int iter )
	// Call the method from R and confirm the right distributions are generated
	// This becomes important and different programming languages adopt different defaults 
	public static double[] testRNG( Parameters pars, int iter ) {
		double [] out = new double[ iter ];
		double shape = pars.getShapeGammaOffspring();
		double beta = pars.getRateTransmit();
		double scale = beta / shape;
		GammaDistribution Ga = new GammaDistribution( RNG, shape, scale );
		
		for( int i = 0; i < iter; ++i ) {
			PoissonDistribution Pois = new PoissonDistribution( Model.RNG, Ga.sample(), 1E-12, 10000000 );
//			out[ i ] = rngPois.nextInt() ;
//			out[ i ] = rngGamma.nextDouble();
//			out[ i ] = Ga.sample();
			out[ i ] = Pois.sample();
		}
		return out ;
	}
	

	////////////////////////////////////////////////////////////////////////////////////////////////
	// updateHospitalsForVaccination
	// 
	public static void updateHospitalsForVaccination( Parameters pars ){
		String scenario = pars.getVaccinationScenario();
		if( scenario.equalsIgnoreCase( "Distance") ) {
			updateHospitalsForVaccinationByDistance( pars );
		}
		else if( scenario.equalsIgnoreCase( "Region") ) {
			updateHospitalsForVaccinationByRegion( pars );
		}
		else if( scenario.equalsIgnoreCase( "Hospital") ) {
			updateHospitalsForVaccinationByHospital( pars );
		}
	}
	
	
	////////////////////////////////////////////////////////////////////////////////
	// updateHospitalsForVaccinationByDistance()
	// update the list of hospitals for vaccination based on the distance from the hospitals where cases have been isolated
	//
	public static void updateHospitalsForVaccinationByDistance( Parameters pars ){
		// hospitals where isolation occurred but vaccination didn't start
		ArrayList<Hospital> hospWithIsolationNoVaccination = new ArrayList<Hospital> ();
		hospWithIsolationNoVaccination.addAll( Model.hospitalsCaseIsolated );
		hospWithIsolationNoVaccination.removeAll( Model.hospitalsVaccinationImplemented );
		for( Hospital h : hospWithIsolationNoVaccination ) {
			hospitalsVaccinationImplemented.add( h );
			h.setDayVaccinationStarted( Step.currentDay );
		}
		// decide target pool
		//vaccines are given to new hospitals when cases are  this is updated in the isolate() method in the Step class
		double distCutoff = pars.getVaccinationTargetRadius();
		ArrayList<Hospital> hospitalsSearched = new ArrayList<Hospital>();
		hospitalsSearched.addAll( Model.hospitals );
		hospitalsSearched.addAll( Model.uninfectedHospitals );
		hospitalsSearched.removeAll( Model.hospitalsVaccinationImplemented );
		for( Hospital hosp : hospWithIsolationNoVaccination ) { // target hospitals are within a certain distance from the hospital where cases are isolated
			double myLon = hosp.getLongitude();
			double myLat = hosp.getLatitude();
			for( Hospital h : hospitalsSearched ) {
				double d = util.getDistance( myLat, h.getLatitude(), myLon, h.getLongitude() ); 
				if( d < distCutoff ) {
					hospitalsVaccinationImplemented.add( h );
					h.setDayVaccinationStarted( Step.currentDay );
				}
			}
		}
		if( pars.getDebug() > 4  ) {
			for( Hospital hosp : hospitalsVaccinationImplemented ) {
				System.out.printf( "day=%.1f, day vacc started =%.1f, vaccinated hospital ID = %d\n", Step.currentDay, 
				hosp.getDayVaccinationStarted(), hosp.getID() );
			}
		}
	}
	

	////////////////////////////////////////////////////////////////////////////////
	// updateHospitalsForVaccinationByArea()
	// update the list of hospitals for vaccination based on the region they belong to
	//
	public static void updateHospitalsForVaccinationByRegion( Parameters pars ){
		
		// hospitals where isolation occurred but vaccination didn't start
		ArrayList<Hospital> hospWithIsolationNoVaccination = new ArrayList<Hospital> ();
		hospWithIsolationNoVaccination.addAll( Model.hospitalsCaseIsolated );
		hospWithIsolationNoVaccination.removeAll( Model.hospitalsVaccinationImplemented );
		for( Hospital h : hospWithIsolationNoVaccination ) {
			hospitalsVaccinationImplemented.add( h );
			h.setDayVaccinationStarted( Step.currentDay );
		}
		// decide target pool
		//vaccines are given to new hospitals when cases are  this is updated in the isolate() method in the Step class
		ArrayList<Hospital> hospitalsSearched = new ArrayList<Hospital>();
		hospitalsSearched.addAll( Model.hospitals );
		hospitalsSearched.addAll( Model.uninfectedHospitals );
		hospitalsSearched.removeAll( Model.hospitalsVaccinationImplemented );
		for( Hospital hosp : hospWithIsolationNoVaccination ) { // target hospitals are within a certain distance from the hospital where cases are isolated
			for( Hospital h : hospitalsSearched ) {
				if( hosp.getRegionID() == h.getRegionID() || h.getRegionID() == 8 ) { //same region or Seoul
					hospitalsVaccinationImplemented.add( h );
					h.setDayVaccinationStarted( Step.currentDay );
				}
			}
		}
	}
	
	

	////////////////////////////////////////////////////////////////////////////////
	// updateHospitalsForVaccinationByArea()
	// update the list of hospitals for vaccination based on the region they belong to
	//
	public static void updateHospitalsForVaccinationByHospital( Parameters pars ){
//		Seoul Asan 127.109525587	37.5251582003
//		Samsung Seoul 127.089591626	37.4903457681
//		Yonsei Severance 126.940823731	37.5614071269
//		Seoul National University Hospital 127.00039849	37.5779137772
//		Catholic University Seoul St. Mary's Hospital 127.005862753	37.5023918137
		ArrayList<Hospital> hospitalsSearched = new ArrayList<Hospital>();
		hospitalsSearched.addAll( Model.hospitals );
		hospitalsSearched.addAll( Model.uninfectedHospitals );
		hospitalsSearched.removeAll( Model.hospitalsVaccinationImplemented );
		
		for( Hospital h: hospitalsSearched  ) {
			double lon = h.getLongitude();
			if( 127.10952558 < lon && lon < 127.10952559 ) {
				Model.hospitalsVaccinationImplemented.add( h );
			}
			if( 127.08959162 < lon && lon < 127.08959163 ) {
				Model.hospitalsVaccinationImplemented.add( h );
			}
			if( 126.94082373 < lon && lon < 126.94082374 ) {
				Model.hospitalsVaccinationImplemented.add( h );
			}
			if( 127.0003984 < lon && lon < 127.0003985 ) {
				Model.hospitalsVaccinationImplemented.add( h );
			}
			if( 127.00586275 < lon && lon < 127.00586276 ) {
				Model.hospitalsVaccinationImplemented.add( h );
			}
		}
		if( pars.getDebug() > 3 ) {
			for( Hospital h: hospitalsVaccinationImplemented  ) {
				System.out.printf( "day=%.1f, hospial pop =%d\n", Step.currentDay, h.getPopulationSize() );
			}
		}
	}
	
	
	
}
