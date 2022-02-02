/*   TIDD: A Java(TM) based peptide feature extraction software.
 *
 *   Copyright (C) 2022 BIS Labs, Hanyang Univ. Korea
 *
 *   TIDD is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   TIDD is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *   
 *   You should have received a copy of the GNU General Public License
 *   along with TIDD.  If not, see <http://www.gnu.org/licenses/>.
 */

package TIDD;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.StringTokenizer;
import java.util.Vector;

import dst.AnnoPeakT;
import dst.PeptideT;
import dst.SearchSummaryT;
import dst.SpectrumT;

public class ReadResultFiles {
	public SearchSummaryT ss;
	public Spectrum spectrum;
	public SpectrumT spectrumt;
	
	public boolean isHighResolution;
	public HashMap<String, Double> fixedmod;
	public String fMGFpath;
	private long pre_time;
	public static int specID;
	public static BufferedReader mgfbr;
	
	public boolean itraq;
	
	private int[] pIndex;
	private double[] inten;
	
	public static final double PROTON = 1.007276;
	public static final double H_2_O = 1.007825035*2 + 15.99491463;
	public static final double N_H_3 = 1.007825035*3 + 14.003074;
	
	public static final double aminoacid_masss[]={
		71.03711,	0,			103.00919, 115.02694,	129.04259, 
		147.06841,	57.02146,	137.05891, 113.08406,	0, 
		128.09496,	113.08406,	131.04049, 114.04293,	0, 
		97.05276,	128.05858,	156.10111, 87.03203,	101.04768,
		0,			99.06841,	186.07931, 0,			163.06333,	0};
	
	public ReadResultFiles(long pre_time) {
		// TODO Auto-generated constructor stub
		ss = null;
		specID = 0;
		itraq = false;
		mgfbr = null;
		this.pre_time = pre_time;
		fixedmod = new HashMap<String, Double>();
		spectrum = new Spectrum();
	}
	
//	public List<PeptideSpecMatchT> readMODa_fdr(String filepath){
//		BufferedReader br = null;
//		List<PeptideSpecMatchT> psmlist = new ArrayList<PeptideSpecMatchT>();
//		
//		ProxDB pdb= null;
//		TreeMap<String, Integer> proxMap = null;
//		try {
//			System.out.print((System.currentTimeMillis() - pre_time)/1000 + " seconds\t");
//			System.out.println("Read protein database.");
//			
//			pdb= new ProxDB();
//			pdb.readFasta(filepath.substring(0,filepath.lastIndexOf('\\')+1) + ss.database);
//			proxMap = new TreeMap<String, Integer>();
//			for(int i=0; i<pdb.size(); i++){
//				proxMap.put(pdb.get(i).getAccession(), i+1);
//			}
//		} catch (Exception e1) {
//			// TODO Auto-generated catch block
////			e1.printStackTrace();
//		}
//		
//		System.out.print((System.currentTimeMillis() - pre_time)/1000 + " seconds\t");
//		System.out.println("Start reading a moda result file.");
//		
//		try{
//			BufferedReader preview = new BufferedReader( new FileReader(filepath) );
//			HashSet<String> ptmList = new HashSet<>();
//	//		preview.readLine();
//			String line;
//			
//			int filesize = 0;
//			
//			while( (line = preview.readLine()) != null )
//			{
//				filesize++;
////				System.out.println(line);
//				StringTokenizer tok= new StringTokenizer(line, "\t"); 
//				if( tok.countTokens() < 9 ) continue;		
//				String sptname=tok.nextToken(), sptindex=tok.nextToken(), obsvmw=tok.nextToken(), charge=tok.nextToken(); 
//				String calcmw=tok.nextToken(), delta=tok.nextToken(), strscore=tok.nextToken(), prob=tok.nextToken();
//				String pept=tok.nextToken(), prot=tok.nextToken(), pepstart=tok.nextToken();
//	
//				int pep_len= pept.length();
//				for(int i=2; i<pep_len-2; i++){
//					if( !Character.isLetter( pept.charAt(i) ) ){											
//						StringBuffer mod= new StringBuffer();
//						mod.append( pept.charAt(i-1) );
//						mod.append( pept.charAt(i) );
//						int k;
////						for(k=i+1; !Character.isAlphabetic( pept.charAt(k) ); k++){
//						for(k=i+1; k < pept.length() ; k++){
//							if( Character.isAlphabetic( pept.charAt(k) ) )
//								break;
//							mod.append( pept.charAt(k) );
//						}
//						ptmList.add(mod.toString());
//						i = k-1;
//					}
//				}		
//			}				
//			preview.close();
//			
//			ArrayList<String> ptmMat = new ArrayList(ptmList);
//			Collections.sort( ptmMat );
//			
//			File fMGF = new File(fMGFpath);
//			
//			br = new BufferedReader(new FileReader(filepath));
//			
//			int filelineNum = 0;
//			ScanContainer scanDB = null;
//			
////			int cumul_scan = 0;
//			
//			if(fMGF.isDirectory() ){
//				File[] fractions = fMGF.listFiles();
//				int mgfid = -1;
//				String prevFrac = "";
//				while( (line=br.readLine()) != null ){
//					
//					filelineNum++;
//					if( filelineNum%100 == 0 ){
//						System.out.print((System.currentTimeMillis() - pre_time)/1000 + " seconds\t");
//						System.out.println("Processed PSM count : " + filelineNum);
//						
//					}
//					
//					if( line.startsWith("SpectrumFile") )
//						continue;
//					
//					StringTokenizer tokens = new StringTokenizer(line, "\t");
//					String mgf_frac = tokens.nextToken();
//					
//					
////					int specID = Integer.parseInt(tokens.nextToken())-1 - cumul_scan;
//					int specID = Integer.parseInt(tokens.nextToken())-1;
//					double observedMW = Double.parseDouble(tokens.nextToken());
//					int charge = Integer.parseInt(tokens.nextToken());
//					double calcMW = Double.parseDouble(tokens.nextToken());
//					double deltaM = Double.parseDouble(tokens.nextToken());
//					int score = Integer.parseInt(tokens.nextToken());
//					double prob = Double.parseDouble(tokens.nextToken());
//					
//					String tmp = tokens.nextToken();
//					char prevA = tmp.charAt(0);
//					char postA = tmp.charAt(tmp.length()-1);
//					String peptide = tmp.substring(2, tmp.length()-2);
//					String prot = tokens.nextToken();
//					StringTokenizer protst = new StringTokenizer(prot, ";");
//					String temp_prot = "";
//					while( protst.hasMoreTokens() ){
//						String temp = protst.nextToken();
//						if( temp.startsWith("ENST") )
//							temp_prot += temp; 
//						else{
//							if( temp.contains("_U") || temp.contains("_B") )
//								temp_prot = temp + temp_prot;
//							else
//								temp_prot += temp;
//						}
//					}
//					prot = temp_prot;
//					
//					String peppos = tokens.nextToken();
//					
//					PeptideSpecMatchT psm = new PeptideSpecMatchT();
//					psm.observed_mass = observedMW;
//					psm.cs = charge;
//					
//					if( mgf_frac.compareTo(prevFrac) != 0 ){
//						System.out.print((System.currentTimeMillis() - pre_time)/1000 + " seconds\t");
//						System.out.println("Read specturm file : " + mgf_frac);
//						
//						prevFrac = mgf_frac;
//						mgfid++;
//						scanDB = new ScanContainer(charge);
//						if( fractions[mgfid].getName().toLowerCase().endsWith(".mgf") ){
//							scanDB.parseFromMGF(fractions[mgfid].getAbsolutePath());
//						}
//						System.out.print((System.currentTimeMillis() - pre_time)/1000 + " seconds\t");
//						System.out.println("End reading");
//					}
//					
//					
//					PeptideT pept = new PeptideT();
//					int pep_len= peptide.length();
//					
//					StringBuffer seq= new StringBuffer();
//					for(int i=0; i<pep_len; i++){
//						if( Character.isLetter( peptide.charAt(i) ) ){
//							seq.append(peptide.charAt(i));
//						}
//					}
//					pept.sequence = seq.toString();
//					
////					ArrayList<String> modifications = new ArrayList<String>();
//					
//					pept.modified_chk = new boolean[pept.sequence.length()];
//					pept.modified_position = new double[pept.sequence.length()];
//					for(int k=0; k<pept.sequence.length(); k++){
//						pept.modified_position[k] = 0;
//						pept.modified_chk[k] = false;
//					}
//					pept.modset = new ArrayList<ModificationT>();
//					
//					int site= 0;
//					for(int i=0; i<pep_len; i++){
//						if( Character.isLetter( peptide.charAt(i) ) ){
//							seq.append(peptide.charAt(i));
//							site++;
//						}else{
//							StringBuffer diff= new StringBuffer();
////							diff.append( peptide.charAt(i-1) );
//							diff.append( peptide.charAt(i) );
//							int k = 0;
//							
//							
//							
////							try{
//							for(k=i+1; k < peptide.length() ; k++){
//								if( Character.isAlphabetic( peptide.charAt(k) ) )
//									break;
//								diff.append( peptide.charAt(k) );
//							}
////							}catch(Exception e){
////								System.out.println("sequence : " + peptide);
////								System.out.println("index : " + k);
////								System.exit(0);
////							}
//							
////							modifications.add(mod.toString());
//							ModificationT mod = new ModificationT();
//							mod.diff = Double.parseDouble(diff.toString());
//							mod.site = peptide.charAt(i-1) + "";
//							mod.name = mod.site + diff.toString();
//							mod.position = "ANYWHERE";
//							mod.varID = Collections.binarySearch(ptmMat, mod.name);
//							if( !ss.ptms.containsKey(mod.varID) ){
//								ss.ptms.put(mod.varID, mod);
//							}
////							try{
//							pept.modified_chk[site-1] = true;
////							}catch(Exception e){
////								System.out.println(pept.sequence);
////								System.out.println(pept.modified_chk.length);
////								System.exit(0);
////							}
//							pept.modified_position[site-1] += mod.diff;
//							pept.modset.add(mod);
//							i = k-1;
//						}
//					}	
//										
//					pept.mw = ""+calcMW;
//					pept.deltaM = "" + deltaM;
//					pept.score = prob + "";
//					pept.p_n_res = prevA+""+postA;
//					
//					if( pdb != null ){
//						pept.protein = prot;
//						try{
//							pept.protein_id = proxMap.get(prot);
//						}catch(Exception e ){
//							pept.protein_id = -1;
//						}
//						pept.pep_start = Integer.parseInt(peppos.split("~")[0]);
//					}else{
//						pept.protein = null;
//					}
//					
//					
//					
////					try{
//					ss.psmOffset.add(scanDB.get(specID).getOffset());
////					}catch(Exception e){
//////						System.out.print((System.currentTimeMillis() - pre_time)/1000 + " seconds\t");
//////						System.out.println("Error!");
//////						
//////						System.out.println(line);
//////						System.out.println(specID);
//////						System.out.println(fractions[mgfid].getName());
////////						System.err.println(e.getMessage());
//////						BufferedReader readmgf = new BufferedReader(new FileReader(fractions[mgfid]));
//////						String templine;
//////						int scancnt = 0;
//////						while( ( templine = readmgf.readLine() ) != null ){
//////							if( templine.startsWith("BEGIN") )
//////								scancnt++;
//////						}
//////						System.out.println("Mgf id : " + mgfid);
//////						System.out.println("Scan cnt : " + scancnt);
//////						System.out.println("ScanDB size : " + scanDB.size());
//////						System.exit(1);
////						cumul_scan += scanDB.size();
////						mgfid++;
////						scanDB = new ScanContainer(charge);
////						if( fractions[mgfid].getName().toLowerCase().endsWith(".mgf") ){
////							scanDB.parseFromMGF(fractions[mgfid].getAbsolutePath());
////						}
////						
////						specID -= cumul_scan;
////						System.out.print((System.currentTimeMillis() - pre_time)/1000 + " seconds\t");
////						System.out.println("End reading");
////						ss.psmOffset.add(scanDB.get(specID).getOffset());
////					}
//					double frag_tol = 0;
//					if( ss.fragTolUnit.compareTo("ppm") == 0 ){
//						frag_tol = ss.fragTol * calcMW / 1000000;
//					}
//					psm.spectrum = spectrum.getProcessedSpectrum(fractions[mgfid], (int)scanDB.get(specID).getOffset(), (observedMW - 1.007276)*charge, observedMW, frag_tol);
//					
//					
//					
//					psm.peptides = new PeptideT[1];
//					psm.peptides[0] = pept;
//					psmlist.add(psm);
//				}
//			}
//			br.close();
//		}catch(IOException e){
//			System.err.println(e.getMessage());
//		}
//		return psmlist;
//	}
//	
//	public SearchSummaryT getSummary(){
//		return ss;
//	}
	
	public void readParam(String filepath) throws IOException{
		String line = "";
		BufferedReader br = new BufferedReader(new FileReader(filepath));
		ss = new SearchSummaryT();
		
		ss.aaMass = aminoacid_masss.clone();
		
		while( (line=br.readLine()) != null ){
			if( line.startsWith("#") || line.isEmpty() )
				continue;
			StringBuilder sb = new StringBuilder(line);
			
			
			String name = sb.substring(0, sb.indexOf("="));
//			System.out.println(sb.toString());
			String value = sb.substring(sb.indexOf("=")+1);
			
			
			switch (name) {
			case "Spectra":
//				spectra = value;
				ss.msmsdata = value;
				break;
			case "Instrument":
//				instrument = value;
				break;
			case "Fasta":
//				fasta = value;
				ss.database = value;
				break;
			case "PPMTolerance":
//				ppm = value;
				ss.peptTol = Double.parseDouble(value);
				ss.peptTolUnit = "ppm";
				break;
			case "PeptTolerance":
//				pepmasstol = value;
				ss.peptTol = Double.parseDouble(value);
				ss.peptTolUnit = "Da";
				break;
			case "FragTolerance":
//				fragtol = value;
				ss.fragTol = Double.parseDouble(value);
				ss.fragTolUnit = "Da";
				break;
			case "BlindMode":
//				blind = value;
				break;
			case "MinModSize":
//				minsize = value;
				break;
			case "MaxModSize":
//				maxsize = value;
				break;
			case "Enzyme":
//				enzyme = value;
				break;
			case "enzyme_constraint min_number_termini":
//				ntt = value;
				break;
			case "MissedCleavage":
//				miscleavage = value;
				break;
			case "ADD":
				StringTokenizer st = new StringTokenizer(value, ",");
				String site = st.nextToken().trim();
				String delta = st.nextToken().trim();
 				fixedmod.put(site,Double.parseDouble(delta));
				break;
			case "iTRAQSearch":
				if( value.startsWith("ON") )
					itraq = true; 
				else
					itraq = false;
				break;
			case "HighResolution":
				if( value.startsWith("ON") )
					isHighResolution = true;
				else
					isHighResolution = false;
				break;
			default:
				break;
			}
		}
		br.close();
	}
	
//	public static Map<String, Double> mass_map;
//	private static void initMap(){ 
////		C[Num]H[Num]N[Num]O[Num]S[Num]P[Num]Br[Num]Cl[Num]Fe[Num]
//		mass_map = new HashMap<String, Double>();
//		mass_map.put(new String("C"), 12.0107);
//		mass_map.put(new String("H"), 1.00794);
//		mass_map.put(new String("N"), 14.0067);
//		mass_map.put(new String("O"), 15.9994);
//		mass_map.put(new String("S"), 32.065);
//		mass_map.put(new String("Br"), 79.904);
//		mass_map.put(new String("Cl"), 35.453);
//		mass_map.put(new String("Fe"), 55.845);
//	}
//	
//	public void readMSGFParam(String filepath) throws IOException{
//		initMap();
//		
//		ss.aaMass = aminoacid_masss.clone();
//		
//		String line = "";
//		BufferedReader br = new BufferedReader(new FileReader("D:/占쌜억옙 占쏙옙占쏙옙/160621_MSGFtoMODp/Mods.txt"));
//		
//		while( (line=br.readLine()) != null ){
//			line = line.trim();
//			if( line.startsWith("#") || line.isEmpty() )
//				continue;
////			StringBuilder sb = new StringBuilder(line);
//			
//			if( line.startsWith("NumMods=") ){
////				String name = sb.substring(0, sb.indexOf("="));
//////				System.out.println(sb.toString());
////				String value = sb.substring(sb.indexOf("=")+1);	
//			}else{
//				StringTokenizer st = new StringTokenizer(line,"");
//				String first_tk = st.nextToken();
//				
//				String[] stok_mod = first_tk.split(",");
//				
//				double mass = 0;
//				try{
//					mass = Double.parseDouble(stok_mod[0]);
//				}catch(Exception e){
//					StringBuffer sb_tmp = new StringBuffer();
//					StringBuffer sb_dit = null;
//					String str = stok_mod[0];
//					String atom = null;
//					
//					for( int i=0;i < str.length(); i++ ){
//						if( Character.isAlphabetic(str.charAt(i)) ){
//							if( atom != null && sb_dit != null ){
////								System.out.println(mass_map.get(atom) + "\t" + atom + "\t" + 2);
//								mass += mass_map.get(atom)*Double.parseDouble(sb_dit.toString());
//								atom = null;
//								sb_dit = null;
//							}else if( mass_map.containsKey(sb_tmp.toString()) ){
////								System.out.println(mass_map.get(sb_tmp.toString()) + "\t" + sb_tmp.toString() + "\t" + 3);
//								mass += mass_map.get(sb_tmp.toString())*1;
//								atom = null;
//								sb_dit = null;
//								sb_tmp = new StringBuffer();
//							}
//							
//							sb_tmp.append(str.charAt(i));
//						}else{
//							if( atom == null ){
//								atom = sb_tmp.toString();
//								sb_tmp = new StringBuffer();
//								sb_dit = new StringBuffer();
//							}
//							
//							sb_dit.append(str.charAt(i));
//						}
//						
//					}
//					if( atom != null && sb_dit != null ){
////						System.out.println(mass_map.get(atom) + "\t" + atom + "\t" + 1);
//						mass += mass_map.get(atom)*Double.parseDouble(sb_dit.toString());
//					}
//				}
////				System.out.println(mass);
//				
//				String residue = stok_mod[1];
//				
//				if( stok_mod[2].contains("fix") ){
//	 				fixedmod.put(residue,mass);
//				}
//			}
//			
//		}
//		br.close();
//	}
	
	public void setMGFpath(String mgf){
		fMGFpath = mgf;
	}
	
	public Vector<AnnoPeakT> getAnnotatedPeaks(SpectrumT spectrumT, PeptideT peptide, SearchSummaryT ss) {
		int i, matIndex = -1;
		
		Vector<AnnoPeakT> annoPeaks = new Vector<AnnoPeakT>();
		
		//	占싱뤄옙占쏙옙占쏙옙占쏙옙 占쏙옙占쏙옙 y-ion占쏙옙 占쏙옙占쏙옙占쏙옙占쏙옙 찾占승댐옙.
		for( i = 0; i < peptide.theo_y_ions.length; i++ ){
			
			//	y-ion占쏙옙 m/z占쏙옙占쏙옙 占쏙옙치(tolerance 占쏙옙 占싱삼옙 占쏙옙占싱놂옙占쏙옙 占십댐옙)占싹댐옙 peak占쏙옙 찾占쏙옙 占쏙옙占쏙옙占싼댐옙.
			matIndex = getMatchedPeak(spectrumT, peptide.theo_y_ions[i], 0.00);
			
			if( matIndex > -1 ){
				AnnoPeakT annoPeak = new AnnoPeakT();
				annoPeak.type = "ION|Y1";
				annoPeak.mass = spectrumT.peakMassList.get(matIndex);
				annoPeak.inten = spectrumT.peakIntensityList.get(matIndex);
				annoPeak.priority = 1;
				annoPeak.anno = "y"+ (peptide.sequence.length() - i);
				annoPeaks.add(annoPeak);
				
			}
			
			//	charge state占쏙옙 2 占싱삼옙占쏙옙 占쏙옙占쏙옙占� m/z占쏙옙占쏙옙 占쏙옙占쏙옙臼占� match 占쏙옙占쏙옙占쏙옙 확占쏙옙占싹울옙 占쏙옙치占싹댐옙 占쏙옙占� 占쏙옙占쏙옙占싼댐옙.
			for( int tcs = 2; tcs <= spectrumT.cs; tcs++ ){
				double cutOff = ( tcs == spectrumT.cs )? 0.5 : 0.000;
				matIndex = getMatchedPeak( spectrumT, (peptide.theo_y_ions[i]+PROTON*(tcs-1))/tcs, cutOff);
				if( matIndex > -1 ){
					AnnoPeakT annoPeak = new AnnoPeakT();
					annoPeak.type = "ION|Y"+tcs;
					annoPeak.mass = spectrumT.peakMassList.get(matIndex);
					annoPeak.inten = spectrumT.peakIntensityList.get(matIndex);
					annoPeak.priority = 3;
					if( tcs == 2 ) annoPeak.anno = "y" + (peptide.sequence.length() - i) + "++";
					//	annotation 표占썩에 charge state占쏙옙 표占쏙옙占싼댐옙.
					else annoPeak.anno = "y" + (peptide.sequence.length() - i) + "(" + tcs + "+)";
					annoPeaks.add(annoPeak);
				}
			}
			
			//	neutral loss占쏙옙 占싹어났占쏙옙 占쏙옙占쏙옙占� m/z占쏙옙占쏙옙 占쏙옙占쏙옙臼占� match 占쏙옙占쏙옙占쏙옙 확占쏙옙占싼댐옙.
			AnnoPeakT LOSS, NH3;
			LOSS = new AnnoPeakT();
			NH3 = new AnnoPeakT();
			LOSS.inten = 0;
			NH3.inten = 0;
			matIndex = getMatchedPeak(spectrumT, peptide.theo_y_ions[i]-H_2_O, .1);
			//System.out.println(matIndex + "  " + (peptide.theo_y_ions[i]-H_2_O));
			if( matIndex > -1 ){
				LOSS.type = "LOSS|Y";
				LOSS.mass = spectrumT.peakMassList.get(matIndex);
				LOSS.inten = spectrumT.peakIntensityList.get(matIndex);
				LOSS.anno = "y" + (peptide.sequence.length() - i) + "-H2O";
			}
			matIndex = getMatchedPeak(spectrumT, peptide.theo_y_ions[i]-N_H_3, .1);
			if( matIndex > -1 ){
				NH3.type = "LOSS|Y";
				NH3.mass = spectrumT.peakMassList.get(matIndex);
				NH3.inten = spectrumT.peakIntensityList.get(matIndex);
				NH3.anno = "y" + (peptide.sequence.length() - i) +"-NH3";
			}
			if( LOSS.inten < NH3.inten ){
				LOSS.type = NH3.type;
				LOSS.mass = NH3.mass;
				LOSS.inten = NH3.inten;
				LOSS.anno = NH3.anno;
			}
			if( LOSS.inten != 0 ){
				AnnoPeakT annoPeak = new AnnoPeakT();
				annoPeak.type = LOSS.type;
				annoPeak.mass = LOSS.mass;
				annoPeak.inten = LOSS.inten;
				annoPeak.priority = 6;
				annoPeak.anno = LOSS.anno;
				annoPeaks.add(annoPeak);
			}
			
		}
		
		//System.out.println();
		
		//	占싱뤄옙占쏙옙占쏙옙占쏙옙 占쏙옙占쏙옙 b-ion占쏙옙 占쏙옙占쏙옙占쏙옙占쏙옙 찾占승댐옙.	
		for( i = 0 ; i < peptide.theo_b_ions.length; i++ ){
			
			//	b-ion占쏙옙 m/z占쏙옙占쏙옙 占쏙옙치(tolerance 占쏙옙 占싱삼옙 占쏙옙占싱놂옙占쏙옙 占십댐옙)占싹댐옙 peak占쏙옙 찾占쏙옙 占쏙옙占쏙옙占싼댐옙.
			matIndex = getMatchedPeak(spectrumT, peptide.theo_b_ions[i], 0.00);
			if( matIndex > -1 ){
				AnnoPeakT annoPeak = new AnnoPeakT();
				annoPeak.type = "ION|B1";
				annoPeak.mass = spectrumT.peakMassList.get(matIndex);
				annoPeak.inten = spectrumT.peakIntensityList.get(matIndex);
				annoPeak.priority = 2;
				annoPeak.anno = "b" + (i+1);
				annoPeaks.add(annoPeak);
			}
			
			//	charge state占쏙옙 2 占싱삼옙占쏙옙 占쏙옙占쏙옙占� m/z占쏙옙占쏙옙 占쏙옙占쏙옙臼占� match 占쏙옙占쏙옙占쏙옙 확占쏙옙占싹울옙 占쏙옙치占싹댐옙 占쏙옙占� 占쏙옙占쏙옙占싼댐옙.
			for( int tcs = 2; tcs <= spectrumT.cs; tcs++ ){
				double cutOff = ( tcs == spectrumT.cs )? .5 : .000;
				
				matIndex = getMatchedPeak(spectrumT, (peptide.theo_b_ions[i]+PROTON*(tcs-1))/tcs , cutOff);
				if( matIndex > -1 ){
					AnnoPeakT annoPeak = new AnnoPeakT();
					annoPeak.type = "ION|B"+ tcs;
					annoPeak.mass = spectrumT.peakMassList.get(matIndex);
					annoPeak.inten = spectrumT.peakIntensityList.get(matIndex);
					annoPeak.priority = 3;
					if( tcs == 2 )	annoPeak.anno = "b" + (i+1) + "++";
					else	annoPeak.anno = "b" + (i+1) + "(" + tcs + "+)";
					annoPeaks.add(annoPeak);
				}
			}
			//	neutral loss占쏙옙 占싹어났占쏙옙 占쏙옙占쏙옙占� m/z占쏙옙占쏙옙 占쏙옙占쏙옙臼占� match 占쏙옙占쏙옙占쏙옙 확占쏙옙占싼댐옙.
			AnnoPeakT LOSS, NH3;
			LOSS = new AnnoPeakT();
			NH3 = new AnnoPeakT();
			LOSS.inten = 0;
			NH3.inten = 0;
			matIndex = getMatchedPeak(spectrumT, peptide.theo_b_ions[i]-H_2_O, .1);
			if( matIndex > -1 ){
				LOSS.type = "LOSS|B";
				LOSS.mass = spectrumT.peakMassList.get(matIndex);
				LOSS.inten = spectrumT.peakIntensityList.get(matIndex);
				LOSS.anno = "b" + (i+1) + "-H2O";
			}
			matIndex = getMatchedPeak(spectrumT, peptide.theo_b_ions[i]-N_H_3, .1);
			if( matIndex > -1 ){
				NH3.type = "LOSS|B";
				NH3.mass = spectrumT.peakMassList.get(matIndex);
				NH3.inten = spectrumT.peakIntensityList.get(matIndex);
				NH3.anno = "b" + (i+1) +"-NH3";
			}
			if( LOSS.inten < NH3.inten ){
				LOSS.type = NH3.type;
				LOSS.mass = NH3.mass;
				LOSS.inten = NH3.inten;
				LOSS.anno = NH3.anno;
			}
			if( LOSS.inten != 0 ){
				AnnoPeakT annoPeak = new AnnoPeakT();
				annoPeak.type = LOSS.type;
				annoPeak.mass = LOSS.mass;
				annoPeak.inten = LOSS.inten;
				annoPeak.priority = 6;
				annoPeak.anno = LOSS.anno;
				annoPeaks.add(annoPeak);
			}
		}
		
		//	a-ion占쏙옙 占쏙옙占쏙옙占쏙옙 확占쏙옙占싼댐옙.
		for( i = 1; i < 4 && i < peptide.theo_b_ions.length; i++ ){
			matIndex = getMatchedPeak(spectrumT, peptide.theo_b_ions[i]  - 27.9949, .1);
			if( matIndex > -1 ){
				AnnoPeakT annoPeak = new AnnoPeakT();
				annoPeak.type = "ION|A";
				annoPeak.mass = spectrumT.peakMassList.get(matIndex);
				annoPeak.inten = spectrumT.peakIntensityList.get(matIndex);
				annoPeak.priority = 4;
				annoPeak.anno = "a" + (i+1);
				annoPeaks.add(annoPeak);
			}
		}
		
		//	Immonium ion占쏙옙 占쏙옙占쏙옙占쏙옙 확占쏙옙		
		HashSet<Character> imm = new HashSet<Character>();
		for( i = 0; i < peptide.sequence.length(); i++ ){
			if( imm.contains(peptide.sequence.charAt(i)) )
				continue;
			matIndex = getMatchedPeak(spectrumT,
					ss.getAAMass(peptide.sequence.charAt(i)) + peptide.modified_position[i] - 27, .1);
			if( matIndex > -1 ){
				imm.add(peptide.sequence.charAt(i));
				AnnoPeakT annoPeak = new AnnoPeakT();
				annoPeak.type = "IMM|B";
				annoPeak.mass = spectrumT.peakMassList.get(matIndex);
				annoPeak.inten = spectrumT.peakIntensityList.get(matIndex);
				annoPeak.priority = 6;
				if( peptide.modified_position[i] == 0 )
					annoPeak.anno = peptide.sequence.charAt(i) + "*";
				else if( peptide.modified_position[i] > 0 )
					annoPeak.anno = peptide.sequence.charAt(i) + "+" + Math.round(peptide.modified_position[i]) + "*";
				else
					annoPeak.anno = "" + peptide.sequence.charAt(i) + Math.round(peptide.modified_position[i]) + "*";
				annoPeaks.add(annoPeak);
			}
		}
		if( annoPeaks.size() > 0 )
			annoPeaks = qSort(annoPeaks);
		return annoPeaks;
	}
	
	//	Vector 占쏙옙占쏙옙 占쏙옙 AnnoPeak占쌘료구占쏙옙占쏙옙 占쏙옙치占쏙옙 swap占싼댐옙.
	private Vector<AnnoPeakT> swap(Vector<AnnoPeakT> vector, int i, int j){
		vector.set(j, vector.set(i, vector.get(j)));
		return vector;
	}
	//	Annotation占쏙옙 占썼열占쏙옙 mass占쏙옙 priority占쏙옙占쏙옙占쏙옙 sorting占싼댐옙.
	public Vector<AnnoPeakT> qSort(Vector<AnnoPeakT> annoPeaks){
		if( annoPeaks.size() == 0 )
			return new Vector<AnnoPeakT>();
		else if( annoPeaks.size() == 1 )
			return annoPeaks;
		else if( annoPeaks.size() == 2 ){
			if( compareAnnoPeaks(annoPeaks.get(0), annoPeaks.get(1)) > 0 )
				annoPeaks = swap(annoPeaks, 0, 1);
			return annoPeaks;
		}
		AnnoPeakT pivot = annoPeaks.lastElement();
		Vector<AnnoPeakT> lefts = new Vector<AnnoPeakT>();
		Vector<AnnoPeakT> rights = new Vector<AnnoPeakT>();
		
		for(int i = 0; i < annoPeaks.size()-1; i++){
			if( compareAnnoPeaks(annoPeaks.get(i), pivot) >= 0 )
				rights.add(annoPeaks.get(i));
			else
				lefts.add(annoPeaks.get(i));
		}
		lefts = qSort(lefts);
		rights = qSort(rights);
		lefts.add(pivot);
		lefts.addAll(rights);
		return lefts;
	}
	
	//	특占쏙옙 m/z占쏙옙占쏙옙 占쏙옙占쏙옙 peak占쏙옙 spectrum占쏙옙占쏙옙 찾占쏙옙 index占쏙옙 占싼곤옙占쌔댐옙.
	public int getMatchedPeak(SpectrumT spectrumT, double d, double cutOff) {
		int x =spectrum.getMatchedPeak(spectrumT, d, cutOff);
		return x;
	}
	
	private int compareAnnoPeaks(AnnoPeakT a, AnnoPeakT b){
		double x = a.mass, y = b.mass;
		int xp = a.priority, yp = b.priority;
		
		if( x > y )
			return 1;
		else if( x == y ){
			if( xp > yp )
				return 1;
			else if ( xp == yp )
				return 0;
			else 
				return -1;
		}	
		else 
			return -1;
	}
	
	public SpectrumT getProcessedSpectrum(String title, double MW, double pmz, double frag_tol){
		String line = "";
		
		double baseInt = .0, baseMZ = .0, TIC = .0, SME = .0;
		double tarMass = .0, tarInten = .0, tolerance = frag_tol;
//		System.out.println(title);
//		System.out.println("OK");
		
		spectrumt = new SpectrumT();
		//	mgf file 占쏙옙占쏙옙 占쏙옙치占쏙옙 찾占싣곤옙 mass占쏙옙 intensity占쏙옙占쏙옙 占싻억옙占쏙옙灌占�.
		try {
			//long os = offset;
			
//			spectrum.maxME = -99999999.0;
//			spectrum.minME = 99999999.0;
			
			spectrumt.peakMassList = new Vector<Double>();
			spectrumt.peakIntensityList = new Vector<Double>();
			
			spectrumt.originalPeakMassList = new Vector<Double>();
			spectrumt.originalPeakIntensityList = new Vector<Double>();
			
			int first_scan = 0;
			
			//rf.seek(os);
			while( (line = mgfbr.readLine()) != null ){
//				System.out.println(line);
//				if( line.startsWith("BEGIN") ){
					specID++;
					if( line.startsWith("TITLE")  ){
//						if( first_scan < 8 )
//						if( line.endsWith("6380") ){
//							System.out.println("\tFirst Scan after reading PSM: " + line);
//							System.exit(0);
//						}
						
						first_scan++;
//						System.out.println(line);
//						System.out.println(line.endsWith(title));
//						System.out.println(line.contains(title));
//						System.out.println(title);
						if( line.endsWith(title)){
							while( (line=mgfbr.readLine()) != null ){
								line.trim();
								String[] splt = line.split("\\s+");
								if(splt.length == 2){
									try{
										Double.parseDouble(splt[0]);
										break;
									}catch(Exception e){
										
									}
								}
								
							}
							break;
						}
//						if( line.contains("FN01_N181T182_180min_10ug_C2_033114.10054") ){
//							System.out.println("End:END");
//							System.exit(0);
//						}
					}
//					if( specid == specID ){
//						while( (line=mgfbr.readLine()) != null ){
//							line.trim();
//							String[] splt = line.split("\\s+");
//							if(splt.length == 2){
//								try{
//									Double.parseDouble(splt[0]);
//									break;
//								}catch(Exception e){11
//									
//								}
//							}
//							
//						}
//						break;
//					}
//				}
			}
//			System.out.println();
//			System.out.println("Start: "+line);
			
			while(!(line).contains("END") && line!=null){
				double mass = .0, intensity = .0;
				line.trim();
				if( line.isEmpty() ){
					line = mgfbr.readLine();
					continue;
				}
				
				//	mass占쏙옙 intensity占쏙옙 占쏙옙占쏙옙占싹곤옙 占쌍댐옙 占쏙옙占쏙옙占싶울옙占쏙옙 占쏙옙占쏙옙 占쌨아울옙 (占쏙옙占쏙옙占싶울옙 占쏙옙占쏙옙 占쌘울옙 占쌩곤옙占쏙옙 占쏙옙占쏙옙 占쌕댐옙 占쏙옙理� 占쏙옙占쏙옙占쏙옙 占쌕뤄옙占쏙옙 占쏙옙占쏙옙)
				StringTokenizer tokenizer = new StringTokenizer(line);
				mass = Double.parseDouble(tokenizer.nextToken());
//				System.out.println(line);
				intensity = Double.parseDouble(tokenizer.nextToken());
				
				/*
				 * 
				 * 2014-01-15占쏙옙占쏙옙 占쌘듸옙
				 * //	m/z 占쏙옙占쏙옙 intensity占쏙옙 占쏙옙占쏙옙占쏙옙 split 占쏙옙占쌘몌옙 占싱울옙占싹울옙 占쏙옙 占싻몌옙
				String[] temp = line.split("\\t");
				if( temp.length != 2 )	temp = line.split(" ");
				if( temp.length != 2 )	break;
				mass = Double.parseDouble(temp[0]);
				intensity = Double.parseDouble(temp[1]);
				*/
				
				/*
				 * 	2014-05-29占쌘듸옙 占쌩곤옙
				 * 
				 *  占쏙옙占쏙옙 占쏙옙占쏙옙占쌍댐옙 spectrum占쏙옙 占쏙옙占쏙옙歐占� 占쏙옙占쏙옙 占쏙옙占쏙옙 占쏙옙占쏙옙占쏙옙
				 */
				spectrumt.originalPeakMassList.add(mass);
				spectrumt.originalPeakIntensityList.add(intensity);
				
				if( mass > MW )	break;
				
				if( mass < 1 || Math.abs(mass-pmz) < 3 ){
					line = mgfbr.readLine();
					continue;
				}
				
				
				if( intensity != 0 ){
					if( (mass-tarMass) < tolerance ){
						double sum = tarInten + intensity;
						tarMass = tarMass * (tarInten / sum) + mass * (intensity / sum);
						tarInten += intensity;
						spectrumt.peakMassList.add(tarMass);
						spectrumt.peakIntensityList.add(tarInten);
						SME += mass-tarMass;
//						if( spectrum.minME > mass-tarMass ){
//							spectrum.minME = mass-tarMass;
//						}if( spectrum.maxME < mass-tarMass ){
//							spectrum.maxME = mass-tarMass;
//						}
					}else{
						tarMass = mass;
						tarInten = intensity;
						spectrumt.peakMassList.add(tarMass);
						spectrumt.peakIntensityList.add(tarInten);
					}
					
					TIC += intensity;
					if( tarInten > baseInt ){
						baseInt = tarInten;
						baseMZ = tarMass;
					}
				}
				line = mgfbr.readLine();
//				System.out.println(line);
				//System.out.println(mass + " " +  intensity);
			}
			
			spectrumt.peakMassList.trimToSize();
			spectrumt.peakIntensityList.trimToSize();
			
			spectrumt.TIC = TIC;
			spectrumt.peaksCount = spectrumt.peakMassList.size();
			spectrumt.baseIntensity = baseInt;
			spectrumt.baseMZ = baseMZ;
			spectrumt.mzLimit = ((int)(MW/100)+1)*100;
//			spectrum.SME = SME;
			//System.out.println(spectrum.peaksCount);
			
			setPeakRankInLocalArea();
			
			//for(int i = 0; i< spectrum.localRank.length; i ++){
			//	System.out.println(spectrum.localRank[i]);
			//}
			
		} catch (IOException e) {
			
			e.printStackTrace();
		}
		return spectrumt;
	}

	private void setPeakRankInLocalArea() {
		int i = 0;
		int localSize = 100; 
		
		//int j = 0;
		spectrumt.localRank = new int[spectrumt.peaksCount];
		
		int localIndex = (int)(spectrumt.peakMassList.get(i)/localSize) + 1;
		while( i < spectrumt.peaksCount  ){
			int go = 0;
			
			pIndex = new int[spectrumt.peaksCount];
			inten = new double[spectrumt.peaksCount];
			while( localIndex > spectrumt.peakMassList.get(i)/localSize){
				pIndex[go] = i;
				inten[go] = spectrumt.peakIntensityList.get(i);
				i ++;
				go++;
				if(i >= spectrumt.peaksCount)
					break;
			}
			sortSpectrumByIntencity( 0, go-1);
			for(int k = 0; k < go; k++){
				spectrumt.localRank[pIndex[k]] = go - (k+1);
			}
			//for(int k =0;k < go; k++){
			//	System.out.println(j + " : " + spectrum.localRank[k] + "\t|||||\t" + spectrum.peakIntensityList.get(j));
			//	j++;
			//}
			//System.out.println("========================="+spectrum.peakMassList.get(i-1)+"=======================");
			
			localIndex++;
		}
		return;
	}
//	占쏙옙 index占쏙옙, 占쏙옙 intensity占쏙옙 swap占싼댐옙. 
	private void swap(int i, int j){
		//System.out.println(pIndex[i] +" : " + inten[i] + "\t\t" + pIndex[j] + " : " + inten[j]);
		int index = pIndex[i];
		pIndex[i] = pIndex[j];
		pIndex[j] = index;
		
		double temp = inten[i];
		inten[i] = inten[j];
		inten[j] = temp;
		//System.out.println(pIndex[i] +" : " + inten[i] + "\t\t" + pIndex[j] + " : " + inten[j]);
		//System.out.println();
		return;
	}
	 
	public static String strippedseq(String seq){
		String stripseq="";
		for(int i =0; i< seq.length() ; i++){
			if(Character.isAlphabetic(seq.charAt(i))){
				stripseq+=seq.charAt(i);
			}
			else if(seq.charAt(i)=='-'){
//				stripseq+=seq.charAt(i);
			}
		}
		
		
		return stripseq;
	}
	 
	 //	quick sort占쏙옙 intensity 占쏙옙占쏙옙占쏙옙 index占쏙옙 占쏙옙占쏙옙占싼댐옙.
	private void sortSpectrumByIntencity(int left, int right){
		double pivot;
		int i, j;
		//System.out.println("Right = " + right + " , left = " + left);
		if( right > left ){
			if(right - left < 2)
				return;
			
			pivot = inten[left+(right-left+1)/2];
			
			i = left-1;
			j = right+1;
			while(true){
				
				while( inten[++i] < pivot );
				while( inten[--j] > pivot );
				
				if( i >= j ) break;
				
				swap(i, j);
			}
			sortSpectrumByIntencity(left, j);
			sortSpectrumByIntencity(j+1, right);
			
		}
	}
}
