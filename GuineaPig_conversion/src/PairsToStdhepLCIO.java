import hep.io.stdhep.StdhepBeginRun;
import hep.io.stdhep.StdhepEndRun;
import hep.io.stdhep.StdhepEvent;
import hep.io.stdhep.StdhepWriter;

import hep.lcio.event.LCIO;
import hep.lcio.implementation.event.ILCCollection;
import hep.lcio.implementation.event.ILCEvent;
import hep.lcio.implementation.event.IMCParticle;
import hep.lcio.implementation.io.LCFactory;
import hep.lcio.io.LCWriter;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Scanner;
import java.util.StringTokenizer;


/**
 * Converter for a GuineaPig pairs.dat file to stdhep or slcio
 * Author: Anne Schuetz
 */

public class PairsToStdhepLCIO {

	/**
	 * Calls functions ToLCIO() or ToStdhep() depending on command line input.
	 * @param args String array of command line input
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		Start();
		if (args.length < 1) Usage();
		if (args.length %2 != 0 && !(args[0].equals("-h") || args[0].equals("--help"))){
			System.out.println("Please check your arguments!\n"
					+ "I guess you forgot to set a flag... Type for the USAGE:\n"
					+ ">> java -cp bin:lib/* PairsToStdhepLCIO -h / --help");
			System.exit(1);
		}

		//Input and Output files:
		boolean inputfile_set = false;
		boolean outputfile_set = false;
		String input_filename = null;
		String output_filename = null;
		int total_number_of_particles = 0;

		//Find the input and output file:
		for (int i = 0; i < args.length; i++){
			if ( args[i].equals("-h") || args[i].equals("--help")) Usage();
			if ( args[i].equals("-i")){
				input_filename = args[i+1];
				inputfile_set = true;
			}
			if ( args[i].equals("-o")){
				output_filename = args[i+1];
				outputfile_set = true;
			}
		}
		//Check if input and output filenames are set by command line input:
		if (!inputfile_set || !outputfile_set){
			System.out.println("You didn't give an inputfile/outputfile. Please try again!\n");
			Usage();
		}
		//Check if input file with given name exists:
		File inputfile = new File(input_filename);
		if (!inputfile.exists()) {
			System.out.println("Input file " + input_filename + " does not exist!");
			System.exit(1);
		}
		//CountLines: counts lines in file, will exit the program if the file is empty:
		total_number_of_particles = CountLines(input_filename);
		
		//Create the followings for reading in text out of input file:
		FileInputStream pairs_file = new FileInputStream(inputfile);
		BufferedReader pairs = new BufferedReader(new InputStreamReader(
				pairs_file));
	
		//Check if output file already exists -> if yes, program is exited:
		File outputfile = new File(output_filename);
		if (outputfile.exists()) {
			System.out.println("Output file " + output_filename + " already exists!\n"
					+ "Please pick another output filename!");
			System.exit(1);
		}
		//Split output filename into name and format ending:
		int dot = output_filename.lastIndexOf(".");
		String file_format = output_filename.substring(dot+1);
		String output_name = output_filename.substring(0,dot);

		//Default values for the optional arguments:	
		int runnum = 1;
		int nmax = total_number_of_particles;
		double pT_cut_low = 0.0D;
		double pT_cut_high = 999.9D;
		double Theta_cut_low = 0.0D;
		double Theta_cut_high = 2.0D*Math.PI;

		//Save values for the optional input values:
		for (int i = 0; i < args.length; i++){
			if ( args[i].equals("-r")){
				runnum = Integer.parseInt(args[i+1]);
			}
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
		//Program exited, if values for p_T or Theta are negative:
		if (pT_cut_low < 0    || pT_cut_high < 0    ||
		    Theta_cut_low < 0 || Theta_cut_high < 0 ||
		    nmax < 0) {
			System.out.println("Please give positive values for the pT cuts, Theta cuts and the maximum number of particles!");
			System.exit(1);
		}
		

		boolean MoreOutputfiles = false;
		
		if (file_format.equals("stdhep")) {
			/* If given number of particles (nmax) is below total number of particles in input file,
			 * the choice is given to create several output files with each a new set of nmax particles
			 */
			if (nmax < total_number_of_particles) {
				MoreOutputfiles = YesNo_MoreOutputfiles();
			}
			//Call function to convert into stdhep:
			ToStdhep(output_name, pairs, total_number_of_particles, nmax, MoreOutputfiles, pT_cut_low, pT_cut_high, Theta_cut_low, Theta_cut_high);
			pairs_file.close();
		}
		else if (file_format.equals("slcio")){
			/* If given number of particles (nmax) is below total number of particles in input file,
			 * the choice is given to create several output files with each a new set of nmax particles
			 */
			if (nmax < total_number_of_particles) {
				MoreOutputfiles = YesNo_MoreOutputfiles();
			}
			//Call function to convert into slcio:
		 	ToLCIO(output_name, pairs, runnum, nmax, MoreOutputfiles, pT_cut_low, pT_cut_high, Theta_cut_low, Theta_cut_high);
		 	pairs_file.close(); 
		}
		else {
			System.out.println("Unknown file format! Please type output.slcio or output.stdhep!");
			System.exit(1);
		}

	}//end main

	/**
	 * @param outputFilename Name of output file
	 * @param pairs BufferedReader for input pairs.dat file
	 * @param tot_num Total number of particles in the input file
	 * @param _nmax Maximum number of particles per output file
	 * @param More_outputfiles Boolean choice whether several output files shall be created with each _nmax particles
	 * @param pT_cut_low Value for lower limit of p_T
	 * @param pT_cut_high Value for higher limit of p_T
	 * @param Theta_cut_low Value for lower limit of theta
	 * @param Theta_cut_high Value for higher limit of theta
	 */
	public static void ToStdhep(String outputFilename, BufferedReader pairs, int tot_num, int _nmax, boolean More_outputfiles,
			double pT_cut_low, double pT_cut_high, double Theta_cut_low, double Theta_cut_high) {
		
		//Double array of values per line in the input file:
		double[] values = new double[7];

		int _n = 0; //Number of particles
		int _eventnum = 1; //For counting up output files
		
		//Arrays for StdhepEvent values:
		int[] _fst = new int[_nmax];
		int[] _id = new int[_nmax];
		int[] _jmo = new int[2 * _nmax];
		int[] _jda = new int[2 * _nmax];
		double[] _p = new double[5 * _nmax];
		double[] _v = new double[4 * _nmax];
				
		//Dummy values:
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
			//Create a new StdhepWriter with the given output filename:
			w = new StdhepWriter(outputFilename+".stdhep", "Stdhep events",
					"converted from a GuineaPig pairs.dat file", 10);
			w.setCompatibilityMode(false);
			w.writeRecord(new StdhepBeginRun(nevtreq, nevtgen, nevtwrt,
					stdecom, stdxsec, stdseed1, stdseed2));
		}	
		//Catch the exception if opening the output file didn't work:
		catch (java.io.IOException e) {
			System.err.println("Error opening file: " + outputFilename + ".stdhep");
			e.printStackTrace();
			System.exit(1);
		}

		try {
			_n = 0;
			_eventnum = 1;
			
			//Loop over lines in input file:
			String line;
			while ((line = pairs.readLine()) != null) {	
				
				//If more output files are wished, open new output file and create new event:
				if (_n >= _nmax && More_outputfiles){
					
					//The new output filename has an number extension (_eventnum is counted up for every new output file):
					New_outputFilename = outputFilename + "_" + Integer.toString(_eventnum);
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
						//Create a new StdhepWriter with the given output filename:
						w = new StdhepWriter(New_outputFilename+".stdhep", "Stdhep events",
								"converted from a GuineaPig pairs.dat file", 10);
						w.setCompatibilityMode(false);
						w.writeRecord(new StdhepBeginRun(nevtreq, nevtgen, nevtwrt,
								stdecom, stdxsec, stdseed1, stdseed2));
					}
					//Catch the exception if opening the StdhepWriter for the output file didn't work:
					catch (java.io.IOException e) {
						System.err.println("Error opening file: " + New_outputFilename + ".stdhep");
						e.printStackTrace();
						System.exit(1);
					}
					_n = 0;
					
					/* New arrays for StdhepEvent values for each new output file:
					 * (In this way the arrays always have the right size.)
					 */
					int eff_num = _nmax;
					if (tot_num - (_eventnum-1)*_nmax < _nmax){
						eff_num = tot_num - (_eventnum-1)*_nmax; //The remaining particles
					}
					System.out.println("eff_num = " + eff_num);
					_fst = new int[eff_num];
					_id = new int[eff_num];
					_jmo = new int[2 * eff_num];
					_jda = new int[2 * eff_num];
					_p = new double[5 * eff_num];
					_v = new double[4 * eff_num];
				}
						
				int j = 0;
				StringTokenizer st = new java.util.StringTokenizer(line, " ");
				
				//Store values in the current line into array:
				while (st.hasMoreElements()) {
					values[j++] = Double.valueOf(st.nextToken()).doubleValue();
				}
				//Initialise a MCParticle with these values:
				Particle p = new Particle(values);
				p.initialise();
				
				//Skip this particle, if theta or p_T are outside the given limits:
				if ((p.getTheta() < Theta_cut_low || p.getTheta() > Theta_cut_high) 
				       || (p.getPT() < pT_cut_low || p.getPT() > pT_cut_high)){
					continue;
				}
				
				//Store the qualities of the particle into arrays for the stdhep file:
				_fst[_n] = 1; // final state particle
				_id[_n] = (int) p.getPDG();
				_p[0 + 5 * _n] = p.getMomentum()[0]; // px
				_p[1 + 5 * _n] = p.getMomentum()[1]; // py
				_p[2 + 5 * _n] = p.getMomentum()[2]; // pz
				_p[3 + 5 * _n] = p.getEnergy(); // E
				_p[4 + 5 * _n] = p.getMass(); // mass
				_v[0 + 4 * _n] = p.getPosition()[0]; // x
				_v[1 + 4 * _n] = p.getPosition()[1]; // y
				_v[2 + 4 * _n] = p.getPosition()[2]; // z
				
				//Increment the number of particles in this event
				_n++;
				
				 
				//If more output files are wished, open new output file and create new event:
				if (_n >= _nmax && More_outputfiles){
					try{
						System.out.println("\n _eventnum = " + _eventnum + ", _n = " + _n + "\n");
						//Write out arrays to the stdhep event:
						StdhepEvent event = new StdhepEvent(_eventnum, _n, _fst, _id,
							_jmo, _jda, _p, _v);
						//Write out stdhep event into the stdhep file:
						w.writeRecord(event);									
						w.writeRecord(new StdhepEndRun(nevtreq, nevtgen, nevtwrt, stdecom,
							stdxsec, stdseed1, stdseed2));
					
						//Write and close the stdhep file:
						System.out.println("\n DONE! Closing file " + New_outputFilename + ".stdhep with " + _n + " MCParticles.");
						w.close();
					}
					//Catch the exception if writing and closing the stdhep file didn't work:
					catch(java.io.IOException ex){
						System.err.println("Error closing file: " + New_outputFilename + ".stdhep");
		                ex.printStackTrace();
		                System.exit(1);
					}
					//_eventnum is counted up for every new output file:
					_eventnum++;
				}
				//Break if maximum number of particles is reached and no more output files are wished:
				if (_n >=_nmax && !More_outputfiles) break;
			}
			
			try{
				System.out.println("\n _eventnum = " + _eventnum + ", _n = " + _n + "\n");
				//Write out the (last) arrays to the (last) stdhep event:
				StdhepEvent event = new StdhepEvent(_eventnum, _n, _fst, _id,
						_jmo, _jda, _p, _v);
				//Write out the (last) stdhep event into the (last) stdhep file:
				w.writeRecord(event);
				w.writeRecord(new StdhepEndRun(nevtreq, nevtgen, nevtwrt, stdecom,
					stdxsec, stdseed1, stdseed2));
				
				//Write and close the (last) stdhep file:
				System.out.println("\n DONE! Closing file " + New_outputFilename + ".stdhep with " + _n + " MCParticles.");
				w.close();
			}
			//Catch the exception if writing and closing the stdhep file didn't work:
			catch(java.io.IOException ex){
				System.err.println("Error opening/closing file: " + New_outputFilename + ".stdhep");
                ex.printStackTrace();
                System.exit(1);
			}
			
		} catch (java.io.IOException ex) {
			ex.printStackTrace();
		}

	}//end ToStdhep()

	/**
	 * @param outputFilename Name of output file
	 * @param pairs BufferedReader for input pairs.dat file
	 * @param run_num Number of run
	 * @param _nmax Maximum number of particles per output file
	 * @param More_outputfiles Boolean choice whether several output files shall be created with each _nmax particles
	 * @param pT_cut_low Value for lower limit of p_T
	 * @param pT_cut_high Value for higher limit of p_T
	 * @param Theta_cut_low Value for lower limit of theta
	 * @param Theta_cut_high Value for higher limit of theta
	 */
	public static void ToLCIO(String outputFilename, BufferedReader pairs, int run_num, int _nmax, boolean More_outputfiles, double pT_cut_low, double pT_cut_high, double Theta_cut_low, double Theta_cut_high) { 
		
		String New_outputFilename=outputFilename;
		
		//Double array of values per line in the input file:
		double[] values = new double[7];

		// Constants and variables:
		int _n = 0;
		int _pairsnum = 0;
		int _i = 1;
		
		LCWriter lcWriter = null ;
		ILCEvent event = null;	
		ILCCollection GP_pairs = null;
		IMCParticle pair = null;
		
		try{
			//Create an event with a MCParticle Collection with MCParticles called "pair":
			event = new ILCEvent();	
			event.setEventNumber(_i);
			event.setRunNumber(run_num);
			event.setDetectorName("UNKNOWN");
			GP_pairs = new ILCCollection(LCIO.MCPARTICLE);
						
			//Loop over lines in input file:
			String line; 
			while ((line = pairs.readLine()) != null){
				int j = 0 ; StringTokenizer st = new java.util.StringTokenizer(line, " "); 
			
				//Store values in the current line into array:
				while(st.hasMoreElements()){ 
					values[j++] = Double.valueOf(st.nextToken()).doubleValue(); 
				}
				//Initialise a MCParticle with these values:
				Particle p = new Particle(values);
				p.initialise();
				
				//Skip this particle, if theta or p_T are outside the given limits:
				if ((p.getTheta() < Theta_cut_low || p.getTheta() > Theta_cut_high) 
				       || (p.getPT() < pT_cut_low || p.getPT() > pT_cut_high)){
					continue;
				}
				//Make the MCParticle a part of the MCParticle Collection:
				pair = new IMCParticle();
				pair.setPDG(p.getPDG());
				pair.setMass((float) p.getMass());
				pair.setCharge(p.getCharge());
				pair.setMomentum(p.getMomentum());
				pair.setVertex(p.getPosition());
				pair.setGeneratorStatus(p.getGeneratorStatus());
				
				GP_pairs.add(pair);
				
				//Count the numbers of particles and already converted particles up:
				_n++;
				_pairsnum++;
				
				//Break if maximum number of particles is reached and no more output files are wished:
				if(_n >= _nmax && !More_outputfiles){
					break;
				}
				//If more output files are wished, open new output file and create new event:
				if(_n >= _nmax && More_outputfiles) { 
					// open and write new output file 
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
						//Write out this stored event. Then create new event and collection:
						++_i;
						lcWriter = LCFactory.getInstance().createLCWriter() ;
						lcWriter.open(New_outputFilename);
						event.addCollection(GP_pairs, "MCParticle");
						lcWriter.writeEvent(event);
						System.out.println("\n DONE! Closing file "+ New_outputFilename +".slcio with "+_n+" MCParticles.");
						lcWriter.close();
						
						event = new ILCEvent();	
						event.setDetectorName("UNKNOWN");
						event.setRunNumber(run_num);
						event.setEventNumber(_i);
						GP_pairs = new ILCCollection(LCIO.MCPARTICLE);
						
					}
					//Catch exception if writing to file didn't work:
					catch (java.io.IOException e) {
						System.err.println("Error writing to file: " + New_outputFilename + ".slcio");
						e.printStackTrace();
						System.exit(1);
					}
					_n = 0;
				}
												
			} 
			/* Check if several output files were written.
			 * If not, just use the original given output filename without extension:
			 */
			if( _i > 1){
				New_outputFilename = outputFilename + "_" + Integer.toString(_i);
			}
			else if (_i == 1){
				New_outputFilename = outputFilename;
			}
			else {
				//If _i = 0:
				System.out.println("\n Something went wrong with looping over the input file.");
                System.exit(1);
			}
			try{
				//Writing out of (last) set of MCParticles:
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
			//Catch exception if writing to file didn't work:
			catch(java.io.IOException ex){
				System.err.println("Error with opening/closing file: " + New_outputFilename + ".slcio");
				ex.printStackTrace();
                		System.exit(1);
			}
		} catch(java.io.IOException ex){
			ex.printStackTrace(); 
		} 
	}//end of ToLCIO()
	
	/**
	 * @return Returns boolean for choice, whether more output files should be created or not
	 */
	private static boolean YesNo_MoreOutputfiles() {
		boolean more_outputfiles = false;
		Scanner keyboard = null;

		try {
			System.out
					.println("\n The maximum number of events you have given is smaller than the number of events in the input file. \n"
							+ "Do you want me to create several output files (with each nmax events) until the end of the input file is reached? \n"
							+ "(y/n) \n");
			//Keyboard input is stored in String yes_no:
			keyboard = new Scanner(System.in);
			String yes_no = keyboard.nextLine();

			//Change value of boolean more_outputfiles depending on keyboard input:
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
				YesNo_MoreOutputfiles();
			}
			return more_outputfiles;
		} finally {
			//End use of keyboard:
			if (keyboard != null) keyboard.close();
		}

	}

	/**
	 * @param filename Name of file for which lines are to be counted
	 * @return Returns int number of lines
	 * @throws IOException
	 */
	private static int CountLines(String filename) throws IOException {
		//Initialise reading file content:
		BufferedReader reader = new BufferedReader(new FileReader(filename));
		boolean empty = true;
		try {
			int lines = 0;
			//Count up lines till the end of the file:
			while (reader.readLine() != null) {
				empty = false;
				lines++;
			}
			//If file not empty, return the number of lines:
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
	
	private static void Usage() {
		System.out.println("\nPairsToStdhepLCIO: \n"
			+ "Application to convert pairs.dat output files from GuineaPig to stdhep format or slcio format.\n"
			+ "\nCuts on pT and Theta can be applied in the following way: pTcut_low < pT [GeV] < pTcut_high, and thetacut_low < theta [degrees] < thetacut_high. \n"
			+ "With giving an integer number, the number of particles that are to be converted can be defined.\n"
			+ "Passing a run number as an argument will allow to distinguish single simulation files after merging.\n"
			+ "If no run number will be passed the default value 1 will be set as the run number for this file.");
		System.out.println("\nUSAGE: \n"
			+ ">> java -cp bin:lib/*  PairsToStdhepLCIO -i PATH/TO/input.dat -o output<.stdhep / .slcio> <more options> \n");
		System.out.println("\nRequired Arguments:\n");
		System.out.printf("%-25s%s%n","-i:","<GuineaPig input dat file>");
		System.out.printf("%-25s%s%n","-o:","<output filename.stdhep / .slcio>");
		System.out.println("\nOPTIONS:\n");
		System.out.printf("%-25s%s%n","-h / --help:","Usage");
		System.out.printf("%-25s%s%n","-r:","<run number>");
		System.out.printf("%-25s%s%n","-n:","<maximum number of particles that are to be converted>");
		System.out.printf("%-25s%s%n","-pl / --ptcut_low:","<lower limit for pT in GeV>");
		System.out.printf("%-25s%s%n","-ph / --ptcut_high:","<higher limit for pT in GeV>");
		System.out.printf("%-25s%s%n","-tl / --thetacut_low:","<lower limit for theta in degree>");
		System.out.printf("%-25s%s%n","-th / --thetacut_high:","<higher limit for theta in degree>");
		System.out.println("\n For example: \n"
			+ ">> java -cp bin:lib/* PairsToStdhepLCIO -i pairs.dat -o pairs.slcio -r 2 -n 3000 -pl 0.01 -ph 1 -tl 0.2 -th 30");
		System.exit(0);
	}//end Usage()
}//end PairsToSthepLCIO class


class Particle{
	/**
	 * @param qualities Array of values for the MCParticle energy, beta and vertex
	 */
	public Particle(double[] qualities){
		//Take array of particle qualities and store them in class member vectors:
		energy = qualities[0];
		beta = new double[3];
		beta[0] = qualities[1];
		beta[1] = qualities[2];
		beta[2] = qualities[3];
		pos = new double[3];
		pos[0] = qualities[4];
		pos[1] = qualities[5];
		pos[2] = qualities[6];

		//Default values for particle mass, charge and state of simulation:
		GeneratorStatus = 1; //These particles are all generated.
		mass = 0.000510998928; //electron mass
		charge = -1; //electron charge sign
		
	}
	/* Initialise the particle by calculating the momentum, p_T and theta.
	 * Additionally, check whether particle is an electron or positron.
	 */
	public void initialise(){
		if (energy < 0) {
			ChangeToPositron();
		}

		mom = new double[3];
		mom[0] = beta[0] * energy;
		mom[1] = beta[1] * energy;
		mom[2] = beta[2] * energy;

		pT = calculatePT(mom);
		theta = calculateTheta(pT, mom);
	}
	private double calculatePT(double[] momentum){
		double p_T = Math.sqrt(momentum[0] * momentum[0] + momentum[1] * momentum[1]);
		return p_T;
	}
	private double calculateTheta(double pT, double[] momentum){
		double Theta = Math.atan(pT / Math.abs(momentum[2]));
		if (Theta < 0) 	Theta += Math.PI;
		return Theta;
	}
	//Switch PDG number, energy and charge value to correspondent values for positrons
	public void ChangeToPositron(){
		pdg *= -1;
		energy *= -1.0D;
		charge *= -1;
	}
	public final double[] getMomentum(){
		return mom;
	}
	public final double[] getPosition(){
		return pos;
	}
	public final double getEnergy(){
		return energy;
	}
	public final int getPDG(){
		return pdg;
	}
	public final int getCharge(){
		return charge;
	}
	public final double getMass(){
		return mass;
	}
	public final double getPT(){
		return pT;
	}
	public final double getTheta(){
		return theta;
	}
	public final int getGeneratorStatus(){
		return GeneratorStatus;
	}
	
	private double[] mom;
	private double[] beta;
	private double[] pos;
	private double energy;
	private int pdg = 11; //PDG number of an electron
	private int charge = -1; //electron charge
	private double mass = 0.000510998928; //electron mass
	private double pT;
	private double theta;
	private int GeneratorStatus;
}
