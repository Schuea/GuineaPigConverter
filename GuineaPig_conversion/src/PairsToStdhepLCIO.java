import java.io.*;
import java.util.*;

import static java.lang.Math.sqrt;
import static java.lang.Math.toRadians;
import hep.io.stdhep.StdhepWriter;
import hep.io.stdhep.StdhepEvent;
import hep.io.stdhep.StdhepBeginRun;
import hep.io.stdhep.StdhepEndRun;

import hep.lcio.event.*;
import hep.lcio.io.*;
import hep.lcio.implementation.event.*;
import hep.lcio.implementation.io.*;


/**
 * Converter for a GuineaPig pairs.dat file to stdhep or LCIO
 *
 */

public class PairsToStdhepLCIO {

	public static void main(String[] args) throws IOException {
		Start();
		if (args.length < 1) usage();
		if (args.length %2 != 0){
			System.out.println("Please check your arguments!\n"
					+ "I guess you forgot to set a flag... Type -h / --help for the USAGE.");
			System.exit(1);
		}

		//Input and Output files:
		boolean inputfile_set = false;
		boolean outputfile_set = false;
		String input_filename = null;
		String output_filename = null;
		int total_number_of_particles = 0;

		for (int i = 0; i < args.length; i++){
			if ( args[i].equals("-h") || args[i].equals("--help")) usage();
			if ( args[i].equals("-i")){
				input_filename = args[i+1];
				inputfile_set = true;
			}
			if ( args[i].equals("-o")){
				output_filename = args[i+1];
				outputfile_set = true;
			}
		}
		
		if (!inputfile_set || !outputfile_set){
			System.out.println("You didn't give an inputfile/outputfile. Please try again!\n");
			usage();
		}
		
		File inputfile = new File(input_filename);
		if (!inputfile.exists()) {
			System.out.println("Input file " + input_filename + " does not exist!");
			System.exit(1);
		}
		total_number_of_particles = countLines(input_filename);
		
		FileInputStream pairs_file = new FileInputStream(inputfile);
		BufferedReader pairs = new BufferedReader(new InputStreamReader(
				pairs_file));
	
		File outputfile = new File(output_filename);
		if (outputfile.exists()) {
			System.out.println("Output file " + output_filename + " already exists!\n"
					+ "Please pick another output filename!");
			System.exit(1);
		}
		int dot = output_filename.lastIndexOf(".");
		String file_format = output_filename.substring(dot+1);
		String output_name = output_filename.substring(0,dot);

		//Cuts:	
		double pT_cut_low = 0.0D;
		double pT_cut_high = 999.9D;
		double Theta_cut_low = 0.0D;
		double Theta_cut_high = 2.0D*Math.PI;
		int nmax = total_number_of_particles;

		for (int i = 0; i < args.length; i++){
			if ( args[i].equals("-n")){
				nmax = Integer.parseInt(args[i+1]);
			}
			if ( args[i].equals("-pl") || args[i].equals("--ptcut_low")){
				pT_cut_low = Double.parseDouble(args[i+1]);
			}
			if ( args[i].equals("-ph") || args[i].equals("--ptcut_high")){
				pT_cut_high = Double.parseDouble(args[i+1]);
			}
			if ( args[i].equals("-tl") || args[i].equals("--thetacut_low")){
				Theta_cut_low = Math.toRadians(Double.parseDouble(args[i+1]));
			}
			if ( args[i].equals("-th") || args[i].equals("--thetacut_high")){
				Theta_cut_high = Math.toRadians(Double.parseDouble(args[i+1]));
			}
		}
		
		if (pT_cut_low < 0    || pT_cut_high < 0    ||
		    Theta_cut_low < 0 || Theta_cut_high < 0 ||
		    nmax < 0) {
			System.out.println("Please give positive values for the pT cuts, Theta cuts and the maximum number of particles!");
			System.exit(1);
		}
		

		boolean MoreOutputfiles = false;
		
		if (file_format.equals("stdhep")) {
			if (nmax < total_number_of_particles) {
				MoreOutputfiles = YesNo_more_outputfiles();
			}
			ToStdhep(output_name, pairs, pT_cut_low, pT_cut_high, Theta_cut_low, Theta_cut_high, nmax, MoreOutputfiles);
			pairs_file.close();
		}
		else if (file_format.equals("slcio")){
			if (nmax < total_number_of_particles) {
				MoreOutputfiles = YesNo_more_outputfiles();
			}
		 	ToLCIO(output_name, pairs, pT_cut_low, pT_cut_high, Theta_cut_low, Theta_cut_high, nmax, MoreOutputfiles);
		 	pairs_file.close(); 
		}
		else {
			System.out.println("Unknown file format! Please type output.slcio or output.stdhep!");
			System.exit(1);
		}

	}//end main

	private static void ToStdhep(String outputFilename, BufferedReader pairs,
			double pT_cut_low, double pT_cut_high, double Theta_cut_low, double Theta_cut_high, int _nmax, boolean More_outputfiles) {
		
		double[] values = new double[7];

		// Constants and variables
		double mass = 0.000510998928;
		int pdg = 11;
		double energy = -999;
		double[] mom = new double[3];
		double[] beta = new double[3];
		double[] pos = new double[3];

		int _n = 0;
		int _eventnum = 0;
		int _i = 1;
		int[] _fst = new int[_nmax];
		int[] _id = new int[_nmax];
		int[] _jmo = new int[2 * _nmax];
		int[] _jda = new int[2 * _nmax];
		double[] _p = new double[5 * _nmax];
		double[] _v = new double[4 * _nmax];
		
		// Dummy values...
		int nevtreq = 1;
		int nevtgen = 1;
		int nevtwrt = 1;
		float stdecom = 2.F;
		float stdxsec = 2.F;
		double stdseed1 = 137.;
		double stdseed2 = 138.;

		StdhepWriter w = null;
		String New_outputFilename = outputFilename;
		
		try {
			w = new StdhepWriter(outputFilename+".stdhep", "Stdhep events",
					"converted from a GuineaPig pairs.dat file", 10);
			w.setCompatibilityMode(false);
			w.writeRecord(new StdhepBeginRun(nevtreq, nevtgen, nevtwrt,
					stdecom, stdxsec, stdseed1, stdseed2));

		} catch (java.io.IOException e) {
			System.err.println("Error opening file: " + outputFilename);
			e.printStackTrace();
			System.exit(1);
		}

		try {
			_n = 0;
			_eventnum = 0;
			
			String line;
			while ((line = pairs.readLine()) != null) {
												
				if (_eventnum >=_nmax && More_outputfiles){
					New_outputFilename = outputFilename + "_" + Integer.toString(_i);
					File NEWoutputfile = new File(New_outputFilename+".stdhep");
					if (NEWoutputfile.exists()) {
						System.out.println("\nAs I wanted to create several output files with each "+ _nmax 
							+ " MCParticles, I was about to create an output file with the filename " 
							+ New_outputFilename + ".stdhep.\n"
							+ "But such a file already exists.\n"
							+ "Please pick another output filename or move the existing file!");
						System.exit(1);
					}
					try {
						w = new StdhepWriter(New_outputFilename+".stdhep", "Stdhep events",
								"converted from a GuineaPig pairs.dat file", 10);
						w.setCompatibilityMode(false);
						w.writeRecord(new StdhepBeginRun(nevtreq, nevtgen, nevtwrt,
								stdecom, stdxsec, stdseed1, stdseed2));
					} catch (java.io.IOException e) {
						System.err.println("Error opening file: " + outputFilename);
						e.printStackTrace();
						System.exit(1);
					}
					_eventnum = 0;
				}
						
				int j = 0;
				StringTokenizer st = new java.util.StringTokenizer(line, " ");
				while (st.hasMoreElements()) {
					values[j++] = Double.valueOf(st.nextToken()).doubleValue();
				}

				energy = values[0];
				beta[0] = values[1];
				beta[1] = values[2];
				beta[2] = values[3];

				if (energy < 0) {
					pdg *= -1;
					energy *= -1.D;
				}
				mom[0] = beta[0] * energy;
				mom[1] = beta[1] * energy;
				mom[2] = beta[2] * energy;

				double pT = 0.0D;
				double theta = 0.0D;
				pT = Math.sqrt(mom[0] * mom[0] + mom[1] * mom[1]);
				theta = Math.atan(pT / Math.abs(mom[2]));
				if (theta < 0)	theta += Math.PI;

				if ((theta < Theta_cut_low || theta > Theta_cut_high) 
				       || (pT < pT_cut_low || pT > pT_cut_high)){
					continue;
				}

				_fst[_n] = 1; // final state particle
				_id[_n] = (int) pdg;
				_p[0 + 5 * _n] = mom[0]; // px
				_p[1 + 5 * _n] = mom[1]; // py
				_p[2 + 5 * _n] = mom[2]; // pz
				_p[3 + 5 * _n] = energy; // E
				_p[4 + 5 * _n] = mass; // mass
				_v[0 + 4 * _n] = pos[0]; // x
				_v[1 + 4 * _n] = pos[1]; // y
				_v[2 + 4 * _n] = pos[2]; // z
				// increment the number of particles in this event
				_n++;

				StdhepEvent event = new StdhepEvent(_eventnum++, _n, _fst, _id,
						_jmo, _jda, _p, _v);
				w.writeRecord(event);
				// reset the particle per event count
				_n = 0;
		
				if (_eventnum >=_nmax && More_outputfiles){
														
					w.writeRecord(new StdhepEndRun(nevtreq, nevtgen, nevtwrt, stdecom,
							stdxsec, stdseed1, stdseed2));
					try{
						System.out.println("\n DONE! Closing file "+New_outputFilename+".stdhep with "+_eventnum+" MCParticles.");
						w.close();
					}
					catch(java.io.IOException ex){
						System.err.println("Error closing file: " + New_outputFilename + ".stdhep");
		                ex.printStackTrace();
		                System.exit(1);
					}
					_i++;
				}
				if (_eventnum >=_nmax && !More_outputfiles) break;
			}
			
			try{
				
				w.writeRecord(new StdhepEndRun(nevtreq, nevtgen, nevtwrt, stdecom,
					stdxsec, stdseed1, stdseed2));
			
				System.out.println("\n DONE! Closing file "+New_outputFilename+".stdhep with "+_eventnum+" MCParticles.");
				w.close();
			}
			catch(java.io.IOException ex){
				System.err.println("Error opening/closing file: " + New_outputFilename + ".stdhep");
                ex.printStackTrace();
                System.exit(1);
			}
			
		} catch (java.io.IOException ex) {
			ex.printStackTrace();
		}

	}//end ToStdhep()

	private static void ToLCIO(String outputFilename, BufferedReader pairs, double pT_cut_low, double pT_cut_high, double Theta_cut_low, double Theta_cut_high, int _nmax, boolean More_outputfiles) { 
		
		String New_outputFilename=outputFilename;
		
		double[] values = new double[7];

		// Constants and variables
		double mass = 0.000510998928;
		int pdg = 11;
		float charge = -999.0F;
		double energy = -999.0D;
		double[] mom = new double[3];
		double[] beta = new double[3];

		int _n = 0;
		int _pairsnum = 0;
		int _i = 1;
		
		LCWriter lcWriter = null ;
		ILCEvent event = null;	
		ILCCollection GP_pairs = null;
		IMCParticle pair = null;
		
		try{
			event = new ILCEvent();	
			event.setEventNumber(_i);
			event.setRunNumber(1);
			event.setDetectorName("UNKNOWN");
			GP_pairs = new ILCCollection(LCIO.MCPARTICLE);
						
			String line; 
			while ((line = pairs.readLine()) != null){
				int j = 0 ; StringTokenizer st = new java.util.StringTokenizer(line, " "); 
			
				while(st.hasMoreElements()){ 
					values[j++] = Double.valueOf(st.nextToken()).doubleValue(); 
				}
				
				energy = values[0];
				beta[0] = values[1];
				beta[1] = values[2];
				beta[2] = values[3];

				charge = -1.0F;
				if (energy < 0) {
					pdg *= -1;
					energy *= -1.0D;
					charge *= -1.0F;
				}
				mom[0] = beta[0] * energy;
				mom[1] = beta[1] * energy;
				mom[2] = beta[2] * energy;

				double pT = 0;
				double theta = 0;
				pT = Math.sqrt(mom[0] * mom[0] + mom[1] * mom[1]);
				theta = Math.atan(pT / Math.abs(mom[2]));
				if (theta < 0) 	theta += Math.PI;

				if ((theta < Theta_cut_low || theta > Theta_cut_high) 
				       || (pT < pT_cut_low || pT > pT_cut_high)){
					continue;
				}
				
				pair = new IMCParticle();
				
				pair.setPDG(pdg);
				pair.setMass((float) mass);
				pair.setCharge(charge);
				pair.setMomentum(mom[0],mom[1],mom[2]);
				
				GP_pairs.add(pair);
				
				_n++;
				_pairsnum++;
				
				if(_n >= _nmax && !More_outputfiles){
					break;
				}
				if(_n >= _nmax && More_outputfiles) { 
					// open and write new outputfile 
					New_outputFilename = outputFilename + "_" + Integer.toString(_i);
					File NEWoutputfile = new File(New_outputFilename+".slcio");
					if (NEWoutputfile.exists()) {
						System.out.println("\nAs I wanted to create several output files with each "+ _nmax 
							+ " MCParticles, I was about to create an output file with the filename " 
							+ New_outputFilename + ".slcio.\n"
							+ "But such a file already exists.\n"
							+ "Please pick another output filename or move the existing file!");
						System.exit(1);
					}
					try{
						lcWriter = LCFactory.getInstance().createLCWriter() ;
						lcWriter.open(New_outputFilename);
						event = new ILCEvent();	
						event.setDetectorName("UNKNOWN");
						event.setRunNumber(1);
						event.setEventNumber(_i);
						event.addCollection(GP_pairs, "MCParticle");
						lcWriter.writeEvent(event);
						System.out.println("\n DONE! Closing file "+ New_outputFilename +" with "+_n+" MCParticles.");
						lcWriter.close();
					}
					catch (java.io.IOException e) {
						System.err.println("Error opening file: " + New_outputFilename);
						e.printStackTrace();
						System.exit(1);
					}
					
					++_i;
					_n = 0;
				}
												
			} 
			if( _i > 1){
				New_outputFilename = outputFilename + "_" + Integer.toString(_i);
			}
			else if (_i == 1){
				New_outputFilename = outputFilename;
			}
			else {
				System.out.println("\n Something went wrong with looping over the input file.");
                		System.exit(1);
			}
			try{
				File NEWoutputfile = new File(New_outputFilename+".slcio");
					if (NEWoutputfile.exists()) {
						System.out.println("\nI was about to create an output file with the filename " 
							+ New_outputFilename + ".slcio.\n"
							+ "But such a file already exists.\n"
							+ "Please pick another output filename or move the existing file!");
						System.exit(1);
					}
				lcWriter = LCFactory.getInstance().createLCWriter() ;
				lcWriter.open(New_outputFilename);
				event.addCollection(GP_pairs, "MCParticle");
				lcWriter.writeEvent(event);
				lcWriter.close();
				System.out.println("\n DONE! Closing file "+ New_outputFilename +".slcio with "+_n+" MCParticles.\n"
					+"In total, "+ _pairsnum +" MCParticles have been processed.");
			}
			catch(java.io.IOException ex){
				System.err.println("Error with opening/closing file: " + New_outputFilename);
				ex.printStackTrace();
                		System.exit(1);
			}
		} catch(java.io.IOException ex){
			ex.printStackTrace(); 
		} 
	}//end of ToLCIO()
	
	private static boolean YesNo_more_outputfiles() {
		boolean more_outputfiles = false;
		Scanner keyboard = null;

		try {
			System.out
					.println("\n The maximum number of events you have given is smaller than the number of events in the input file. \n"
							+ "Do you want me to create several output files (with each nmax events) until the end of the input file is reached? \n"
							+ "(y/n) \n");
			keyboard = new Scanner(System.in);
			String yes_no = keyboard.nextLine();

			if (yes_no.equals("n")) {
				System.out
						.println("\n Okay, I just convert nmax events and don't care about the rest.");
			} 
			else if (yes_no.equals("y")){
				more_outputfiles = true;
				System.out
						.println("\n Okay, I will create several output files to convert all events of the input file.");
			}
			else {
				System.out.println("\n Please, try again.\n");
				YesNo_more_outputfiles();
			}
			return more_outputfiles;
		} finally {
			if (keyboard != null)
				keyboard.close();
		}

	}

	private static int countLines(String filename) throws IOException {
		BufferedReader reader = new BufferedReader(new FileReader(filename));
		boolean empty = true;
		try {
			int lines = 0;
			while (reader.readLine() != null) {
				empty = false;
				lines++;
			}
			//System.out.println("The input file contains "+lines+" MCParticles."); 
			return (lines == 0 && !empty) ? 1 : lines;
		} finally {
			if (empty) {
				System.out.println("Input file is empty!");
				System.exit(1);
			}
			reader.close();
		}
	}

	private static void Start() {
		System.out.print(String.format("%65s", " ").replace(' ', '*'));
		System.out.printf("%n%-4s%20s%17s%20s%4s%n","****"," ","PairsToStdhepLCIO"," ","****");
		System.out.printf("%-4s%57s%4s%n","****"," Converting GuineaPig pairs.dat files to Stdhep or LCIO ","****");
		System.out.printf("%-4s%18s%21s%18s%4s%n","****"," ","Author:  Anne Schuetz"," ","****");
		System.out.print(String.format("%65s", " ").replace(' ', '*'));
		System.out.println("\n");
	}
	
	private static void usage() {
		System.out.println("\nPairsToStdhepLCIO: \n"
			+ "Application to convert pairs.dat output files from GuineaPig to stdhep format or slcio format.\n"
			+ "Cuts on pT and Theta can be applied in the following way: pTcut_low < pT [GeV] < pTcut_high, and thetacut_low < theta [degrees] < thetacut_high. \n"
			+ "With giving an integer number, the number of particles that are to be converted can be defined.");
		System.out.println("\nUSAGE: \n"
			+ ">> java -cp bin:lib/*  PairsToStdhepLCIO -i PATH/TO/input.dat -o output<.stdhep / .slcio> <more options> \n");
		System.out.println("\nRequired Arguments:\n");
		System.out.printf("%-25s%s%n","-i:","<GuineaPig input dat file>");
		System.out.printf("%-25s%s%n","-o:","<output filename.stdhep / .slcio>");
		System.out.println("\nOPTIONS:\n");
		System.out.printf("%-25s%s%n","-h / --help:","Usage");
		System.out.printf("%-25s%s%n","-n:","<number of particles>");
		System.out.printf("%-25s%s%n","-pl / --ptcut_low:","<lower limit for pT in GeV>");
		System.out.printf("%-25s%s%n","-ph / --ptcut_high:","<higher limit for pT in GeV>");
		System.out.printf("%-25s%s%n","-tl / --thetacut_low:","<lower limit for theta in GeV>");
		System.out.printf("%-25s%s%n","-th / --thetacut_high:","<higher limit for theta in GeV>");
		System.out.println("\n For example: \n"
			+ ">> java -cp bin:lib/* PairsToStdhepLCIO -i pairs.dat -o pairs.slcio -n 3000 -pl 0.01 -ph 1 -tl 0.2 -th 30");
		System.exit(0);
	}//end usage()
}//end PairsToSthepLCIO class

class Particle{
	public Particle(double[] qualities){
		
	}
	public double[] getMomentum(){
		return mom;
	}
	public double getEnergy(){
		return energy;
	}
	public int getPDG(){
		return pdg;
	}
	
	private double[] mom = {0.D};
	private double[] beta = {0.D};
	private double energy = -999.0D;
	private int pdg = 0;
}
