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

import java.io.*;

import java.util.*;

import dst.AnnoPeakT;
import dst.ModificationT;
import dst.PeptideT;
import dst.SearchSummaryT;
import dst.SpectrumT;

public class ExtractFeature {
	
	
	public static final double PROTON = 1.007276;
	public static final double H_2_O = 1.007825035 * 2 + 15.99491463;
	public static final double N_H_3 = 1.007825035 * 3 + 14.003074;
	public static final double phospho_neutral_loss=97.976896;

	public static double NTermOff = PROTON;
	public static double CTermOff = H_2_O + PROTON;
	public static double NTermFixedMod = 0;


	public static final char[] aminos = { 'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'K', 'M', 'F', 'P', 'S',
			'T', 'W', 'Y', 'V', 'U' };
	public static String DECOY = "";
	public static final HashMap<String, Integer> PeptideForm = new HashMap<String, Integer>();
	public static HashMap<String, Integer> numPSMperPep = new HashMap<String, Integer>();
	public static HashMap<String, ArrayList<String>> PeplistperProt = new HashMap<String, ArrayList<String>>();
	public static HashMap<String, String> ProtlistPerPep = new HashMap<String, String>();
	public static int minPepLen = 6;
	public static String header="";

	public static int SpecFile_id=-1;
	public static int Scan_id=-1;
	public static int Charge_id=-1;
	public static int PrecursorMZ_id=-1;
	public static int Peptide_id=-1;
	public static int Protein_id=-1;
	public static int flag_precursorMZ=0;
	public static ArrayList<Integer>add_feature_ids=new ArrayList<Integer>();
	public static ArrayList<String>add_feature_names=new ArrayList<String>();
	
	public static boolean phosphodata=false;
	
	private ReadResultFiles rere;

	private String[] features = { "TIC", "Intensity_of_highest_peak", "Intensity_of_matched_highest_peak",
			"sum_of_intensity_of_all_matched_ions\t" + "sum(intensity_of_y)\t" + "sum(intensity_of_y++)\t"
					+ "sum(intensity_of_y)nl\t" + "sum(intensity_of_b)\t" + "sum(intensity_of_b++)\t"
					+ "sum(intensity_of_b)nl"
//			+ "sum(intensity_of_y++)nl\t"
//			+ "\t"
//			+ "sum(intensity_of_b++)nl"
			, "#all_matched_ions\t" + "#matched_y\t" + "#matched_y++\t" + "#matched_y_nl\t" + "#matched_b\t"
					+ "#matched_b++\t" + "#matched_b_nl"

//			+ "#matched_y++_nl\t"

//			+ "\t"
//			+ "#matched_b++_nl"
			,
			"longest_consecutive_all_ions\t" + "longest_consecutive_y\t" + "longest_consecutive_y++\t"
					+ "longest_consecutive_y_nl\t" + "longest_consecutive_b\t" + "longest_consecutive_b++\t"
					+ "longest_consecutive_b_nl\t"

//			+ "longest_consecutive_y++_nl\t"

//			+ "longest_consecutive_b++_nl\t"
					+ "longest_consecutive_matched",

			"median(all_matched_fragment_ion_errors/matched_fragment)\t"
					+ "median(all_|matched_fragment_ion_errors|/matched_fragment)",
			"mean(all_matched_fragment_ion_errors/matched_fragment)\t"
					+ "mean(all_|matched_fragment_ion_errors|/matched_fragment)",
//			"isotope error corrected mass"
	};

	static int frac_col;
	static int title_col;

	@SuppressWarnings("static-access")
	public void readPSM(String result, String msms, int frac_col, int title_col, int pep_col, int cs_col, int MW_col,
			double frag_tol, String output, long pre_time) throws IOException {
		this.frac_col = frac_col;
		this.title_col = title_col;
		List<String> sortedFile = sortPSMFile(result);

		rere = new ReadResultFiles(pre_time);
		rere.spectrum.TOLERANCE = frag_tol;
		try {
			rere.ss = new SearchSummaryT();
			rere.ss.aaMass = ReadResultFiles.aminoacid_masss.clone();
			PrintWriter err = new PrintWriter(new BufferedWriter(new FileWriter(new File("error.txt"))));

			// Fixed Mod assign
//			rere.ss.aaMass['C'-'A'] += 57.107;

			PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(new File(output))));
			File fMGF = new File(msms);

			File[] fractions;
			if (fMGF.isDirectory()) {
				fractions = fMGF.listFiles();
			} else {
				fractions = new File[] { fMGF };
			}
			int mgfid = -1;
			String prevFrac = "";
			for (String line : sortedFile) {
				String[] tok = line.split("\\t");

				String mgf_frac = tok[frac_col];
				mgf_frac = mgf_frac.substring(0, mgf_frac.indexOf('.'));

				String title = tok[title_col]; // title

				double observedMW = Double.parseDouble(tok[MW_col]); // PrecursorM
				int charge = Integer.parseInt(tok[cs_col]); // charge
				String peptide = tok[pep_col]; // Peptide

				if (mgf_frac.compareTo(prevFrac) != 0) {
					System.out.print((System.currentTimeMillis() - pre_time) / 1000 + " seconds\t");
					System.out.println("Read specturm file : " + mgf_frac);
					ReadResultFiles.specID = 0;
					prevFrac = mgf_frac;
					mgfid++;
					ReadResultFiles.mgfbr = new BufferedReader(new FileReader(fractions[mgfid]));
				}

				boolean[] modified_chk = new boolean[peptide.length()];
				double[] modified_position = new double[peptide.length()];
				for (int k = 0; k < peptide.length(); k++) {
					modified_position[k] = 0;
					modified_chk[k] = false;
				}

				int site = -1;
				StringBuffer seq = new StringBuffer();
				StringBuffer delta = new StringBuffer();
				int i = 0;
				if (peptide.startsWith("+") || peptide.startsWith("-")) {
					StringBuffer ntermM = new StringBuffer();
					for (i = 0; i < peptide.length(); i++) {
						if (!Character.isAlphabetic(peptide.charAt(i))) {
							ntermM.append(peptide.charAt(i));
						} else
							break;
					}

					modified_position[0] += Double.parseDouble(ntermM.toString());
					modified_chk[0] = true;
				}

				for (; i < peptide.length(); i++) {
					if (Character.isAlphabetic(peptide.charAt(i))) {
						delta = new StringBuffer();
						seq.append(peptide.charAt(i));
						site++;
					} else {
						delta.append(peptide.charAt(i));
						int j = 0;
						for (j = i + 1; j < peptide.length(); j++) {
							if (Character.isAlphabetic(peptide.charAt(j)))
								break;
							delta.append(peptide.charAt(j));
						}
						modified_position[site] += Double.parseDouble(delta.toString());
						modified_chk[site] = true;
						i = j - 1;
					}
				}
				SpectrumT spectrumt = rere.getProcessedSpectrum(title, (observedMW - 1.007276) * charge, observedMW,
						frag_tol);

				// write output file

				PeptideT peptidet = new PeptideT();
				peptidet.sequence = seq.toString();
				peptidet.modified_chk = modified_chk;
				peptidet.modified_position = modified_position;
				peptidet.modset = new ArrayList<ModificationT>();

				for (i = 0; i < modified_position.length; i++) {
					if (modified_position[i] != 0) {
						ModificationT mod = new ModificationT();
						mod.diff = modified_position[i];
						mod.site = peptide.charAt(i) + "";
						mod.name = mod.site + modified_position[i];
						mod.position = "ANYWHERE";
						mod.varID = i;
					}
				}

				int pLen = peptidet.sequence.length();
				peptidet.theo_b_ions = new double[pLen];
				peptidet.theo_y_ions = new double[pLen];

				int k = 0;
				double tarMass = NTermOff;

				for (k = 0; k < pLen - 1; k++) {
					tarMass += rere.ss.getAAMass(peptidet.sequence.charAt(k));
					tarMass += peptidet.modified_position[k];
					peptidet.theo_b_ions[k] = tarMass;
				}
				tarMass = CTermOff;

				for (k = pLen - 1; k > 0; k--) {
					tarMass += rere.ss.getAAMass(peptidet.sequence.charAt(k));
					tarMass += peptidet.modified_position[k];
					peptidet.theo_y_ions[k] = tarMass;
				}

				Vector<AnnoPeakT> annos = getAnnotatedPeaks(spectrumt, peptidet, rere.ss);

				double prev_mass = -1;
				for (int j = 0; j < annos.size(); j++) {

					if (annos.get(j).type.compareTo("ION|Y1") == 0) {
					} else if (annos.get(j).type.compareTo("ION|B1") == 0) {
					} else if (annos.get(j).type.startsWith("ION|Y2") && prev_mass != annos.get(j).mass) {
					} else if (annos.get(j).type.startsWith("ION|B2") && prev_mass != annos.get(j).mass) {
					} else if (annos.get(j).type.startsWith("ION|Y3") && prev_mass != annos.get(j).mass) {
					} else if (annos.get(j).type.startsWith("ION|B3") && prev_mass != annos.get(j).mass) {
					} else if (annos.get(j).type.startsWith("ION|Y") && prev_mass != annos.get(j).mass) {
					} else if (annos.get(j).type.startsWith("ION|B") && prev_mass != annos.get(j).mass) {
					} else if (annos.get(j).type.startsWith("LOSS|Y") && prev_mass != annos.get(j).mass) {
						out.println();
					} else if (annos.get(j).type.startsWith("LOSS|B") && prev_mass != annos.get(j).mass) {
					}
					prev_mass = annos.get(j).mass;
				}
			}

			out.close();
			err.close();
		} catch (Exception e) {
			System.out.println(1234);
			e.printStackTrace();
			System.err.println(e.getMessage());
		}
	}
	public void readSearchResult(String result, String msms, double frag_tol, String output, long pre_time) throws IOException {
	
		System.out.println("Sort PSM search results by file fraction and spectrum index!");
	
		List<String> sortedFile = sortFile(result);
		
		System.out.println("End sorting PSM search results! "+sortedFile.size());
		
		int CNT_less2=0;
		int CNT_above2=0;
		
		int psm_cnt = 0;
		double sum_y_rate = 0;
		double sum_b_rate = 0;
		double sum_yall_rate = 0;
		double sum_ball_rate = 0;

		String PreTitle = "";

		rere = new ReadResultFiles(pre_time);
		rere.spectrum.TOLERANCE = frag_tol;
		try {                                                                                                                               
			rere.ss = new SearchSummaryT();
			rere.ss.aaMass = ReadResultFiles.aminoacid_masss.clone();
			
			PrintWriter err = new PrintWriter(new BufferedWriter(new FileWriter(new File("error.txt"))));

			PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(new File(output))));

			out.append(header);
			int Col_CNT=header.split("\\t").length;
			
			out.append("\tTitle\tPrecursorM\tPrecursorMZ\tCalMW\tDeltaMass\tAbsDeltaMass\tPeptideLength\tTIC");
			out.append("\tMaxIntALL\tMaxYionInt\tMaxBionInt\tSumYmatchInt\tSumBmatchInt\tFracYmatchInt\tFracBmatchInt");
			
			out.append("\tSeqCoverYion\tSeqCoverBion\tConsecutiveYion\tConsecutiveBion");
			out.append("\tMassErrMean\tMassErrSD\tNumofAnnoPeaks");//\tMaxMassErr\tMinMassErr");

			//for (int n = 0; n < aminos.length; n++) {
			//	out.append("\t" + aminos[n]);
			//}
			//out.append("\tLnNumofPSM\tLnMaxNumofPepsinBelongedProt\tTriptic\tNumofModif\tTorD\n");

			out.append("\tMissedCleavage\tTriptic\tTorD\tPeptide\n");
			File fMGF = new File(msms);

			if (fMGF.isDirectory()) {

				File[] fractions = fMGF.listFiles();
				int mgfid = -1;
				String prevFrac = "";
				for (String line : sortedFile) {

//					if (line.startsWith("#Spec"))
//						continue;
					psm_cnt++;
//					System.out.println("Process: " + line);
					
					String stre[]=line.split("\\t");
					String mgf_frac=stre[SpecFile_id];
					String scanNum=stre[Scan_id];
					String title=mgf_frac.split(".mgf")[0]+"."+scanNum+"."+scanNum+"."+stre[Charge_id];
					Double Precursor=Double.parseDouble(stre[PrecursorMZ_id]);
					int charge=Integer.parseInt(stre[Charge_id]);
					String peptide=checkpep(stre[Peptide_id]);
					
				
					String protein=stre[Protein_id]; 
						
					Double PrecursorMZ=0.0,ObserveM=0.0;
				//	System.out.println(peptide +" "+title);
					if(flag_precursorMZ==0)
					{
						ObserveM=Precursor;
						PrecursorMZ=ObserveM/charge + 1.007276;
					}
					else
					{
						ObserveM=(PrecursorMZ - 1.007276) * charge;
						PrecursorMZ=Precursor;
					}

					if (!title.equalsIgnoreCase(PreTitle)) {
						PreTitle = title;

						if (mgf_frac.compareTo(prevFrac) != 0) {
							System.out.print((System.currentTimeMillis() - pre_time) / 1000 + " seconds\t");
							System.out.println("Read specturm file : " + mgf_frac);
							ReadResultFiles.specID = 0;
							prevFrac = mgf_frac;
							mgfid++;
							System.out.println("FileName: " + mgf_frac);
							ReadResultFiles.mgfbr = new BufferedReader(new FileReader(msms+mgf_frac));
							
							System.out.print((System.currentTimeMillis() - pre_time) / 1000 + " seconds\t");
							System.out.println("End reading");
						}

						int peplen=(this.extractStripSeq(peptide)).length();
						boolean[] modified_chk = new boolean[peplen];
						double[] modified_position = new double[peplen];
						boolean[] phospho_chk=new boolean[peplen];
						for (int k = 0; k < peplen; k++) {
 							modified_position[k] = 0;
							modified_chk[k] = false;
							phospho_chk[k]=false;
						}

						int site = -1;
						StringBuffer seq = new StringBuffer();
						StringBuffer delta = new StringBuffer();
						int i = 0;
						if (peptide.startsWith("+") || peptide.startsWith("-")) {
							StringBuffer ntermM = new StringBuffer();
							for (i = 0; i < peptide.length(); i++) {
								if (!Character.isAlphabetic(peptide.charAt(i))) {
									ntermM.append(peptide.charAt(i));
								} else
									break;
							}

							modified_position[0] += Double.parseDouble(ntermM.toString());
							modified_chk[0] = true;
						}
						
						for (; i < peptide.length(); i++) {
							if (Character.isAlphabetic(peptide.charAt(i))) {
								delta = new StringBuffer();
								seq.append(peptide.charAt(i));
								site++;
							} else {
								delta.append(peptide.charAt(i));
								int j = 0;
								for (j = i + 1; j < peptide.length(); j++) {
									if (Character.isAlphabetic(peptide.charAt(j)))
										break;
									delta.append(peptide.charAt(j));
								}
								modified_position[site] += Double.parseDouble(delta.toString());
								modified_chk[site] = true;
								
								if(Math.round(Double.parseDouble(delta.toString())*1000)/1000.0==79.966)
									phospho_chk[site]=true;
								
								i = j - 1;
							}
						}
						SpectrumT spectrumt = rere.getProcessedSpectrum(title, ObserveM,
								PrecursorMZ, frag_tol);
						spectrumt.cs=charge;
						
						if(spectrumt.havescan=true)
						{
							PeptideT peptidet = new PeptideT();
							peptidet.sequence = seq.toString();
							peptidet.MissedCleavages = this.Cnt_MissedCleavages(peptidet.sequence);
							peptidet.modified_chk = modified_chk;
							peptidet.modified_position = modified_position;
							peptidet.phospho_chk=phospho_chk;
							peptidet.Nterm_Cterm_Tryptic = this.checkNterm_Cterm_Tryptic(protein, stre[Peptide_id]);
							peptidet.modset = new ArrayList<ModificationT>();

							for (i = 0; i < modified_position.length; i++) {
								if (modified_position[i] != 0) {
									ModificationT mod = new ModificationT();
									mod.diff = modified_position[i];
									mod.site = peptidet.sequence.charAt(i) + "";
									mod.name = mod.site + modified_position[i];
									mod.position = "ANYWHERE";
									mod.varID = i;
								}
							}

							int pLen = peptidet.sequence.length();

							peptidet.theo_b_ions = new double[pLen];
							peptidet.theo_y_ions = new double[pLen];

							int k = 0;
							
							double tarMass = NTermOff;
							

							for (k = 0; k < pLen - 1; k++) {
								tarMass += rere.ss.getAAMass(peptidet.sequence.charAt(k));
								tarMass += peptidet.modified_position[k];
								peptidet.theo_b_ions[k] = tarMass;
								// System.out.println("b-ion " + i + " : " + tarMass);
							}
							
							peptidet.mw=calM(peptide);
						
							tarMass = CTermOff;

							for (k = pLen - 1; k > 0; k--) {
								tarMass += rere.ss.getAAMass(peptidet.sequence.charAt(k));
								tarMass += peptidet.modified_position[k];
								peptidet.theo_y_ions[k] = tarMass;
							}
							peptidet.pLen = pLen;
							
							
							//Vector<AnnoPeakT> annos = getAnnotatedPeaks(spectrumt, peptidet, rere.ss);
							
							Vector<AnnoPeakT> annos = getAnnotatedPeaks_Phopho(spectrumt, peptidet, rere.ss);
							
						//	System.out.println("number of annos:"+annos.size() + " # of annos P:"+annosP.size());
							
							
							peptidet = this.getMatchedPeakInfo(annos, peptidet);
							peptidet.deltaM=(PrecursorMZ - 1.007276) * charge-peptidet.mw;

							HashMap<Character, Integer> tmpPepInfo = this.extractAmino(peptidet.sequence, aminos);

							out.append(line);
							
							int tmp_cnt=line.split("\\t").length;
							
							for(int z=0;z<Col_CNT-tmp_cnt;z++)
							{
								out.append("\t");
							}
							
							out.append("\t"+title+"\t"+ObserveM+"\t"+PrecursorMZ+"\t"+peptidet.mw+"\t"+peptidet.deltaM+"\t"+Math.abs(peptidet.deltaM)+"\t"+pLen+"\t"+Math.log(spectrumt.TIC)/Math.log(2));
							
							Double tmpMaxYInt=peptidet.MaxYInt;
							if(tmpMaxYInt>0)
								tmpMaxYInt=Math.log(tmpMaxYInt)/Math.log(2);
							else tmpMaxYInt=-1.0;
							
							Double tmpMaxBInt=peptidet.MaxBInt;
							if(tmpMaxBInt>0)
								tmpMaxBInt=Math.log(tmpMaxBInt)/Math.log(2);
							else tmpMaxBInt=-1.0;
							
							Double tmpSumYmathedInt=peptidet.SumYmathedInt;
							if(tmpSumYmathedInt>0)
								tmpSumYmathedInt=Math.log(tmpSumYmathedInt)/Math.log(2);
							else
								tmpSumYmathedInt=-1.0;
							
							Double tmpSumBmathchedInt=peptidet.SumBmathchedInt;
							if(tmpSumBmathchedInt!=0)
								tmpSumBmathchedInt=Math.log(tmpSumBmathchedInt)/Math.log(2);
							else
								tmpSumBmathchedInt=-1.0;
													
							out.append("\t"+Math.log(spectrumt.baseIntensity)/Math.log(2)+"\t"+tmpMaxYInt+"\t"+tmpMaxBInt);
							out.append("\t"+tmpSumYmathedInt+"\t"+tmpSumBmathchedInt);
							out.append("\t"+(peptidet.SumYmathedInt/(Double)spectrumt.TIC)+"\t"+(peptidet.SumBmathchedInt/(Double)spectrumt.TIC));
							out.append("\t"+peptidet.SeqCoverY+"\t"+peptidet.SeqCoverB);
							out.append("\t"+peptidet.ConsecutiveY+"\t"+peptidet.ConsecutiveB);
							out.append("\t"+peptidet.MeanMassEr+"\t"+peptidet.SDMassEr+"\t"+peptidet.cnt_mze);

//							for(int n=0;n<aminos.length;n++)
//							{
//								if(tmpPepInfo.containsKey(aminos[n]))
//									out.append("\t"+(tmpPepInfo.get(aminos[n])*1.0/pLen));
//								else
//									out.append("\t0");
//							}
//							
//							out.append("\t"+Math.log(numPSMperPep.get(peptide)));
//							
//							int max=-1; String tmpProt[]=protein.split(";");
//							for(int tmp_i=0;tmp_i<tmpProt.length;tmp_i++)
//							{
//								if(PeplistperProt.get(tmpProt[tmp_i]).size()> max)
//									max=PeplistperProt.get(tmpProt[tmp_i]).size();
//							}
//							out.append("\t"+Math.log(max));
							out.append("\t"+peptidet.MissedCleavages+"\t"+peptidet.Nterm_Cterm_Tryptic);

							out.append("\t"+this.checkProtein(protein));
							out.append("\t"+peptide);
							out.append("\n");
							
//							ArrayList<Double> tmpmz=new ArrayList<Double>();
//							
//							tmpmz=putInfo(tmpmz,peptidet.theo_b_ions);
//							tmpmz=putInfo(tmpmz,peptidet.theo_y_ions);
//							
//							Collections.sort(tmpmz);
//							
//							fw_anno.write("BEGIN IONS\nTITLE="+title+"\n");
//							fw_anno.write("CHARGE="+charge+"\n");
//							fw_anno.write("seq="+peptide+"\n");
//							fw_anno.write("TD="+this.checkProtein(Protein)+"\n");
//							fw_anno.write("THEO_mz=");
//							
//							for(int z=0;z<tmpmz.size();z++)
//							{	
//								fw_anno.write(tmpmz.get(z)+" ");
//							}
//							fw_anno.write("\n");
//							fw_anno.write("Basepeak_intensity="+spectrumt.baseIntensity+"\n");
//							//fw_anno.write("IONTYPE=cs1\n");
//							
//							
//							Double premass=0.0;
//							for(int z=0;z<annos.size();z++)
//							{
//								//if(!(annos.get(z).mass==premass) && annos.get(z).type.startsWith("ION") &&  (annos.get(z).priority==1 || annos.get(z).priority==2 || annos.get(z).priority==3))
//								{
//									fw_anno.write(annos.get(z).mass+" "+annos.get(z).inten+"\t"+annos.get(z).type+"_"+annos.get(z).priority+"\n");
//									premass=annos.get(z).mass;
//								}
//							}
//							fw_anno.write("\nEND IONS\n");
						}

					

					}
				}
			}
			
			out.close();
			err.close();
	//		fw_anno.close();
		} catch (Exception e) {
			System.out.println("No scans in the mgf file");
			e.printStackTrace();
			System.err.println(e.getMessage());
		}

	}

	private String checkpep(String peptide) {
		// TODO Auto-generated method stub
		String tmp=peptide;
		
		char[] AA=peptide.toCharArray();
		
		if((Character.isAlphabetic(AA[0]) || AA[0]=='-')&& (Character.isAlphabetic(AA[AA.length-1]) || AA[AA.length-1]=='-') && AA[1]=='.' && AA[AA.length-2]=='.')
			tmp=peptide.substring(2,AA.length-2);
	 
			
		return tmp;
	}

	private Double calM(String peptide) {
		Double mass=0.0;
		
		StringBuffer delta = new StringBuffer();
		for (int i=0; i < peptide.length(); i++) {
			if (Character.isAlphabetic(peptide.charAt(i))) {
				mass+=rere.ss.getAAMass(peptide.charAt(i));
				delta = new StringBuffer();
			} else {
				delta.append(peptide.charAt(i));
				int j = 0;
				for (j = i + 1; j < peptide.length(); j++) {
					if (Character.isAlphabetic(peptide.charAt(j)))
						break;
					delta.append(peptide.charAt(j));
				}
				mass += Double.parseDouble(delta.toString());
				i = j - 1;
			}
		}
		//System.out.println(peptide+"\t"+(mass+H_2_O));
		return mass+H_2_O;
	}

	private List<String> sortPSMFile(String result) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(new File(result)));
		String line;

		List<String> newLines = new ArrayList<String>();

		br.readLine(); // Header skip

		while ((line = br.readLine()) != null) {
			newLines.add(line);
		}
		br.close();

		Collections.sort(newLines, new Comparator<String>() {

			@Override
			public int compare(String o1, String o2) {
				String[] sp1 = o1.split("\\s+");
				String[] sp2 = o2.split("\\s+");

				String frac1 = sp1[frac_col].substring(0, sp1[frac_col].lastIndexOf('.'));
				String frac2 = sp2[frac_col].substring(0, sp2[frac_col].lastIndexOf('.'));

				if (frac1.compareTo(frac2) < 0) {
					return -1;
				} else if (frac1.compareTo(frac2) > 0) {
					return 1;
				}

				String ind1 = sp1[title_col];
				String ind2 = sp2[title_col];

				if (ind1.compareTo(ind2) < 0)
					return -1;
				else if (ind1.compareTo(ind2) > 0) {
					return 1;
				}
				return 0;
			}
		});

		return newLines;
	}
	private HashMap<String, String> PUT2(String protein, String modPEP,
			HashMap<String, String> protlistperPep2) {

	String prot[]=protein.split("\\;");
	
	if(protlistperPep2.containsKey(modPEP))
	{
		String tmp=protlistperPep2.get(modPEP);
		
		for(int i=0;i<prot.length;i++)
		{
			if(!tmp.contains(prot[i]))
				tmp+=";"+prot[i];
		}
		protlistperPep2.put(modPEP, tmp);
	}
	else
	{
		protlistperPep2.put(modPEP, protein);
	}
	
	
		return protlistperPep2;
	}

private HashMap<String, ArrayList<String>> PUT(String protein, String modPEP,
			HashMap<String, ArrayList<String>> peplistperProt) {
		
	
		String prot[]=protein.split(";");
		
		for(int i=0;i<prot.length;i++)
		{
			if(peplistperProt.containsKey(prot[i]))
			{
				ArrayList<String> tmp=peplistperProt.get(prot[i]);
				if(!tmp.contains(modPEP))
					tmp.add(modPEP);
					
				peplistperProt.put(prot[i],tmp);
			}
			else
			{
				ArrayList<String> tmp=new ArrayList<String>();
				
				tmp.add(modPEP);
					
				peplistperProt.put(prot[i],tmp);
					
			}
			
		}
		return peplistperProt;
	}
	private List<String> sortFile(String result) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(new File(result)));
		
		String line;
		

		List<String> newLines = new ArrayList<String>();

		 header= br.readLine(); // Header skip
	
		while ((line = br.readLine()) != null) {
		
			String str[]=line.split("\\t");
			
//			str[0]=str[0].split("\\/")[str[0].split("\\/").length-1];
//			
//			String tmppline=str[0];
//			for(int i=1;i<str.length;i++)
//				{tmppline+="\t"+str[i];}
//			
//			line=tmppline;
			newLines.add(line);
			
			String pep=str[this.Peptide_id];
			String protein=str[this.Protein_id];
			
			
			
			if(numPSMperPep.containsKey(pep))
			{
				numPSMperPep.put(pep,numPSMperPep.get(pep)+1);
			}
			else
			{
				numPSMperPep.put(pep,1);
			}
			
			PeplistperProt= PUT(protein,pep,PeplistperProt);
			ProtlistPerPep=PUT2(protein, pep, ProtlistPerPep);
			
		}
		br.close();

		Collections.sort(newLines, new Comparator<String>() {

			@Override
			public int compare(String o1, String o2) {
				String[] sp1 = o1.split("\\t");
				String[] sp2 = o2.split("\\t");

				String frac1 = sp1[SpecFile_id];
				String frac2 = sp2[SpecFile_id];

				if (frac1.compareTo(frac2) < 0) {
					return -1;
				} else if (frac1.compareTo(frac2) > 0) {
					return 1;
				}

				int ind1 = Integer.parseInt(sp1[Scan_id]);
				int ind2 = Integer.parseInt(sp2[Scan_id]);

				if (ind1 < ind2)
					return -1;
				else if (ind1 > ind2) {
					return 1;
				}

			return 0;
			}
		});

		return newLines;
	}
	
	public ExtractFeature(String result_file, String msms, String output, long pre_time, double fragTol, String requiredInfo, String decoyPrefix, Boolean phospho) {
		
			try {
				String str[]=requiredInfo.split("\\,");
				this.SpecFile_id=Integer.parseInt(str[0])-1;
				this.Scan_id=Integer.parseInt(str[1])-1;
				this.Charge_id=Integer.parseInt(str[2])-1;
				this.PrecursorMZ_id=Integer.parseInt(str[3])-1;
				this.Peptide_id=Integer.parseInt(str[4])-1;
				this.Protein_id=Integer.parseInt(str[5])-1;
				this.DECOY=decoyPrefix;
				this.phosphodata=phospho;
									
				readSearchResult(result_file, msms, fragTol, output, pre_time);
				
			} catch (IOException e) {
				e.printStackTrace();
			}	 

	}

	private ArrayList<Double> putInfo(ArrayList<Double> tmpmz, double[] theo_ions) {

		for (int i = 0; i < theo_ions.length; i++) {
			if (!tmpmz.contains(theo_ions[i]) && theo_ions[i] > 0)
				tmpmz.add(theo_ions[i]);
			if (!tmpmz.contains(theo_ions[i] - H_2_O) && theo_ions[i] - H_2_O > 0)
				tmpmz.add(theo_ions[i] - H_2_O);
			if (!tmpmz.contains(theo_ions[i] - N_H_3) && theo_ions[i] - N_H_3 > 0)
				tmpmz.add(theo_ions[i] - N_H_3);

		}
		return tmpmz;
	}

	private String checkProtein(String protein) {
		String TorD = "D";

		String str[] = protein.split(";");
		for (int i = 0; i < str.length; i++) {
			if (!str[i].startsWith(DECOY)) {
				TorD = "T";
				break;
			}
		}
		return TorD;
	}

	private PeptideT getMatchedPeakInfo(Vector<AnnoPeakT> annos, PeptideT peptidet) {

		peptidet.MaxBInt = peptidet.MaxYInt = -1.0;
		peptidet.SumBmathchedInt = peptidet.SumYmathedInt = 0.0;

		int yion[] = new int[peptidet.pLen - 1];
		int bion[] = new int[peptidet.pLen - 1];

		yion = this.init(yion);
		bion = this.init(bion);

		ArrayList<Double> mze = new ArrayList<Double>();

		for (int i = 0; i < annos.size(); i++) {
			AnnoPeakT annoPeak = annos.get(i);

			mze.add(annoPeak.abmz);

			if (annoPeak.anno.startsWith("y")) {
				peptidet.SumYmathedInt += annoPeak.inten;
				
				yion[Character.getNumericValue(annoPeak.anno.charAt(1)) - 1] = 1;

				if (peptidet.MaxYInt < annoPeak.inten)
					peptidet.MaxYInt = annoPeak.inten;

			} else if (annoPeak.anno.startsWith("b")) {
				// System.out.println(peptidet.pLen+"
				// "+Character.getNumericValue(annoPeak.anno.charAt(1))+" "+annoPeak.anno+"
				// "+annoPeak.type);
				peptidet.SumBmathchedInt += annoPeak.inten;
				
				bion[Character.getNumericValue(annoPeak.anno.charAt(1)) - 1] = 1;

				if (peptidet.MaxBInt < annoPeak.inten)
					;
				peptidet.MaxBInt = annoPeak.inten;
			}
		}

		peptidet.MeanMassEr = peptidet.sum_abmxe / peptidet.cnt_mze;
		peptidet.SDMassEr = sd(mze);

		peptidet.ConsecutiveB = this.ConsecutiveCheck(bion);
		peptidet.ConsecutiveY = this.ConsecutiveCheck(yion);

		peptidet.SeqCoverY = this.SeqCover(yion);
		peptidet.SeqCoverB = this.SeqCover(bion);

		return peptidet;
	}

	private int[] init(int[] ion) {

		for (int i = 0; i < ion.length; i++)
			ion[i] = 0;

		return ion;
	}

	private double ConsecutiveCheck(int[] ion) {
		int cnt = 0;
		int pre = 0;
		int max = 0;
		for (int i = 0; i < ion.length; i++) {
			if (pre == 0 && ion[i] == 1) {
				if (max < cnt)
					max = cnt;

				cnt = 1;
			} else if (pre == 1 && ion[i] == 1)
				cnt++;

			pre = ion[i];
		}

		if (cnt > max)
			max = cnt;

		return ((max * 1.0) / (ion.length * 1.0));
	}

	public static double mean(ArrayList<Double> table) {
		Double sum = 0.0;
		for (int i = 0; i < table.size(); i++) {
			sum += table.get(i);
		}
		return (sum / (table.size() * 1.0));
	}

	public static double sd(ArrayList<Double> table) {
		double temp = 0;
		double mean = mean(table);

//        if(table.size()!=1)
//        {
		for (int i = 0; i < table.size(); i++) {
			temp += Math.pow(table.get(i) - mean, 2);
		}

		return Math.sqrt(temp / (table.size() - 1));
//        }
//        else
//        	return 0.025;

	}

	private double SeqCover(int[] ion) {

		int cnt = 0;

		for (int i = 0; i < ion.length; i++) {
			if (ion[i] == 1)
				cnt++;
		}

		return ((cnt * 1.0) / (ion.length * 1.0));
	}

	private HashMap<Character, Integer> extractAmino(String sequence, char[] aminos) {

		HashMap<Character, Integer> tmp = new HashMap<Character, Integer>();

		char[] seq = sequence.toCharArray();

		for (int i = 0; i < aminos.length; i++) {
			int cnt = 0;
			for (int j = 0; j < seq.length; j++) {
				if (seq[j] == aminos[i])
					cnt++;
			}
			if (cnt != 0)
				tmp.put(aminos[i], cnt);
		}

		// TODO Auto-generated method stub
		return tmp;
	}

	private Integer Cnt_MissedCleavages(String sequence) {

		int cnt = 0;

		char tmps[] = sequence.toCharArray();

		for (int i = 0; i < tmps.length - 1; i++) {
			if (tmps[i] == 'K' || tmps[i] == 'R') {
				if (i <tmps.length-1) {
					if (tmps[i + 1] != 'P')
						cnt++;
				}
			}
		}
		return cnt;
	}

	private String extractStripSeq(String pep) {

		char[] tmp = pep.toCharArray();
		String seq = "";
		for (int i = 0; i < tmp.length; i++) {
			if (Character.isAlphabetic(tmp[i]))
				seq += tmp[i];
		}
		return seq;
	}

	private int cntModif(boolean[] modified_chk) {
		int count = 0;
		for (int i = 0; i < modified_chk.length; i++) {
			if (modified_chk[i] == true)
				count++;
		}
		return count;
	}

	private String checkNterm_Cterm_Tryptic(String protein, String peptide) {

		String info = "";

		// System.out.println(protein+" aaa "+peptide );

		if ((protein.contains("pre=R") || protein.contains("pre=K") || protein.contains("pre=-"))
				&& (peptide.charAt(peptide.length() - 1) == 'K' || peptide.charAt(peptide.length() - 1) == 'R')
				|| protein.contains("next=-"))
			info = "2";
		else if ((protein.contains("pre=R") || protein.contains("pre=K") || protein.contains("pre=-")))
			info = "1";
		else if ((peptide.charAt(peptide.length() - 1) == 'K' || peptide.charAt(peptide.length() - 1) == 'R'
				|| protein.contains("next=-")))
			info = "0";
		else
			info = "-1";
		
		
		String tmp=peptide;
				
				char[] AA=peptide.toCharArray();
				
				if((Character.isAlphabetic(AA[0]) || AA[0]=='-' ) && (Character.isAlphabetic(AA[AA.length-1]) || AA[AA.length-1]=='-' )&& AA[1]=='.' && AA[AA.length-2]=='.')
				{
					if((AA[2]=='R' || AA[2]=='K' || AA[0]=='-') && (AA[AA.length-3]=='R' ||AA[AA.length-3]=='K' || AA[AA.length-3]=='-'))
						info="2";
					else if((AA[2]=='R' || AA[2]=='K' || AA[0]=='-') || (AA[AA.length-3]=='R' ||AA[AA.length-3]=='K' || AA[AA.length-3]=='-'))
						info="1";
					else
						info="0";
					
				}
				
		return info;
	}

//	private List<String> sortFile(String result, int mod) throws Exception {
//		BufferedReader br = new BufferedReader(new FileReader(new File(result)));
//		String line;
//
//		int pep_id = 0;
//		HashMap<String, ArrayList<String>> pepinfo = new HashMap<String, ArrayList<String>>();
//		List<String> newLines = new ArrayList<String>();
//
//		if (mod == 2) {
//			while ((line = br.readLine()) != null) {
//				String str[] = line.split("\\t");
//				if (line.startsWith("#Spec")) {
//					HeaderInfo = line;
//					for (int i = 0; i < str.length; i++) {
//						if (str[i].equals("Peptide")) {
//							pep_id = i;
//							break;
//						}
//					}
//				} else {
//					newLines.add(line);
//					String pep = str[pep_id];
//					String strip_pep = this.extractStripSeq(pep);
//
//					if (pepinfo.containsKey(strip_pep)) {
//						if (!pepinfo.get(strip_pep).contains(pep))
//							pepinfo.get(strip_pep).add(pep);
//					} else {
//						ArrayList<String> al = new ArrayList<String>();
//						al.add(pep);
//						pepinfo.put(strip_pep, al);
//					}
//				}
//
//			}
//
//			br.close();
//
//		} else if (mod == 3) {
//			ReadXML RX = new ReadXML(result);
//
//			pepinfo = RX.pepinfo;
//			newLines = RX.lists;
//			numPSMperPep = RX.numPSMperPep;
//			PeplistperProt = RX.PeplistperProt;
//			ProtlistPerPep = RX.ProtlistperPep;
//
//			HeaderInfo = RX.Header;
//		}
//
//		Iterator it = pepinfo.keySet().iterator();
//
//		while (it.hasNext()) {
//			String key = (String) it.next();
//			PeptideForm.put(key, pepinfo.get(key).size());
//		}
//
//		Collections.sort(newLines, new Comparator<String>() {
//
//			@Override
//			public int compare(String o1, String o2) {
//				String[] sp1 = o1.split("\\s+");
//				String[] sp2 = o2.split("\\s+");
//
//				String frac1 = sp1[0].substring(0, sp1[0].indexOf('.'));
//				String frac2 = sp2[0].substring(0, sp2[0].indexOf('.'));
//
//				if (frac1.compareTo(frac2) < 0) {
//					return -1;
//				} else if (frac1.compareTo(frac2) > 0) {
//					return 1;
//				}
//
//				int ind1 = Integer.parseInt(sp1[1].split("=")[1]);
//				int ind2 = Integer.parseInt(sp2[1].split("=")[1]);
//
//				if (ind1 < ind2)
//					return -1;
//				else if (ind1 > ind2) {
//					return 1;
//				}
//
//				int snum1 = Integer.parseInt(sp1[2]);
//				int snum2 = Integer.parseInt(sp2[2]);
//
//				if (snum1 < snum2)
//					return -1;
//				else if (snum1 > snum2) {
//					return 1;
//				}
//
//				return 0;
//			}
//		});
//		pepinfo.clear();
//
//		return newLines;
//
//	}
//
//	private void readMODplus(String result, String msms, double frag_tol, String output, long pre_time)
//			throws IOException {
//		// TODO Auto-generated method stub
//		System.out.println("Sort MODplus results by file fraction and spectrum index!");
//		List<String> sortedFile = sortMODplusFile(result);
//		System.out.println("End sorting MODplus results!");
//		int psm_cnt = 0;
//		double sum_y_rate = 0;
//		double sum_b_rate = 0;
//		double sum_yall_rate = 0;
//		double sum_ball_rate = 0;
//		try {
//			PrintWriter err = new PrintWriter(new BufferedWriter(new FileWriter(new File("error.txt"))));
//
//			PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(new File(output))));
//
//			File fMGF = new File(msms);
//
//			if (fMGF.isDirectory()) {
//
//				File[] fractions = fMGF.listFiles();
//				int mgfid = -1;
//				String prevFrac = "";
//				for (String line : sortedFile) {
//
//					if (line.startsWith("Spec"))
//						continue;
//					psm_cnt++;
////					StringTokenizer tokens = new StringTokenizer(line, "\t");
//					String[] tok = line.split("\\t");
//					String mgf_frac = tok[0];
//					mgf_frac = mgf_frac.substring(0, mgf_frac.indexOf('.'));
//
//					int specID = Integer.parseInt(tok[1]);
//					double observedMW = Double.parseDouble(tok[2]);
//					int charge = Integer.parseInt(tok[3]);
////					tok[4];	//	CalculatedMW
////					tok[5];	//	DeltaM
////					tok[6];	//	Score
//					String prob = tok[7]; // Prob
//					String peptide = tok[8]; // Peptide
//					peptide = peptide.substring(2, peptide.length() - 2);
////					tok[9];	//	Protein
//					String Modification = tok[10]; // Modification
//					String scanNum = tok[11]; // scanNum
//					String title = tok[12]; // title
//					title = title.substring(title.indexOf("FN"));
//					if (mgf_frac.compareTo(prevFrac) != 0) {
//						System.out.print((System.currentTimeMillis() - pre_time) / 1000 + " seconds\t");
//						System.out.println("Read specturm file : " + mgf_frac);
//						ReadResultFiles.specID = 0;
//						prevFrac = mgf_frac;
//						mgfid++;
////						System.out.println("FileName: " + fractions[mgfid].getName());
//						ReadResultFiles.mgfbr = new BufferedReader(new FileReader(fractions[mgfid]));
////						scanDB = new ScanContainer(charge);
////						if( fractions[mgfid].getName().toLowerCase().endsWith(".mgf") ){
////							scanDB.parseFromMGF(fractions[mgfid].getAbsolutePath());
////						}
////						System.out.print((System.currentTimeMillis() - pre_time)/1000 + " seconds\t");
////						System.out.println("End reading");
//					}
//
//					boolean[] modified_chk = new boolean[peptide.length()];
//					double[] modified_position = new double[peptide.length()];
//					for (int k = 0; k < peptide.length(); k++) {
//						modified_position[k] = 0;
//						modified_chk[k] = false;
//					}
//
//					int site = -1;
//					StringBuffer seq = new StringBuffer();
//					StringBuffer delta = new StringBuffer();
//					int i = 0;
//					if (peptide.startsWith("+") || peptide.startsWith("-")) {
//						StringBuffer ntermM = new StringBuffer();
//						for (i = 0; i < peptide.length(); i++) {
//							if (!Character.isAlphabetic(peptide.charAt(i))) {
//								ntermM.append(peptide.charAt(i));
//							} else
//								break;
//						}
//
//						modified_position[0] += Double.parseDouble(ntermM.toString());
//						modified_chk[0] = true;
//					}
//
//					for (; i < peptide.length(); i++) {
//						if (Character.isAlphabetic(peptide.charAt(i))) {
//							delta = new StringBuffer();
//							seq.append(peptide.charAt(i));
//							site++;
//						} else {
//							delta.append(peptide.charAt(i));
//							int j = 0;
//							for (j = i + 1; j < peptide.length(); j++) {
//								if (Character.isAlphabetic(peptide.charAt(j)))
//									break;
//								delta.append(peptide.charAt(j));
//							}
//							modified_position[site] += Double.parseDouble(delta.toString());
//							modified_chk[site] = true;
//							i = j - 1;
//						}
//					}
////					System.out.print((System.currentTimeMillis() - pre_time)/1000 + " seconds\t");
////					System.out.println("Start reading MGF and matching peptide with spectrum : " + title);
//					SpectrumT spectrumt = rere.getProcessedSpectrum(title, observedMW, observedMW / charge + 1.007276,
//							frag_tol);
////					System.out.print((System.currentTimeMillis() - pre_time)/1000 + " seconds\t");
////					System.out.println("End reading peptide with spectrum : " + title);
////					System.out.println();
//
//					// write output file
//
//					PeptideT peptidet = new PeptideT();
//					peptidet.sequence = seq.toString();
//					peptidet.modified_chk = modified_chk;
//					peptidet.modified_position = modified_position;
//					peptidet.modset = new ArrayList<ModificationT>();
//
//					for (i = 0; i < modified_position.length; i++) {
//						if (modified_position[i] != 0) {
//							ModificationT mod = new ModificationT();
//							mod.diff = modified_position[i];
//							mod.site = peptide.charAt(i) + "";
//							mod.name = mod.site + modified_position[i];
//							mod.position = "ANYWHERE";
//							mod.varID = i;
//						}
//					}
//
//					int pLen = peptidet.sequence.length();
//					peptidet.theo_b_ions = new double[pLen];
//					peptidet.theo_y_ions = new double[pLen];
//
//					int k = 0;
//					double tarMass = NTermOff + NTermFixedMod;
//
//					for (k = 0; k < pLen - 1; k++) {
//						tarMass += rere.ss.getAAMass(peptidet.sequence.charAt(k));
//						tarMass += peptidet.modified_position[k];
//						peptidet.theo_b_ions[k] = tarMass;
//						// System.out.println("b-ion " + i + " : " + tarMass);
//					}
//					tarMass = CTermOff;
//
//					for (k = pLen - 1; k > 0; k--) {
//						tarMass += rere.ss.getAAMass(peptidet.sequence.charAt(k));
//						tarMass += peptidet.modified_position[k];
//						peptidet.theo_y_ions[k] = tarMass;
//					}
//
//					Vector<AnnoPeakT> annos = getAnnotatedPeaks(spectrumt, peptidet, rere.ss);
////					System.out.println(sum_mxe + "\t" + sum_abmxe + "\t" + cnt_mze);
////					System.out.println(min_mxe + "\t" + max_mxe + "\t" + cnt_mze);
////					System.out.println(min_abmxe + "\t" + max_abmxe + "\t" + cnt_mze);
////					out.append(">");
////					out.append(specID+"");
////					out.append("\t");
////					out.append(scanNum);
////					out.append("\t");
////					out.append(title);
////					out.append("\t");
////					out.append(peptide);
////					out.append("\t");
////					out.append(prob);
////					out.append("\n");
//
////					for jonghun
////					Collections.sort(annos, new Comparator<AnnoPeakT>() {
////
////						@Override
////						public int compare(AnnoPeakT anno1, AnnoPeakT anno2) {
////							if( anno1.type.startsWith("ION") && !anno2.type.startsWith("ION") ){
////								return -1;
////							}else if( !anno1.type.startsWith("ION") && anno2.type.startsWith("ION") ){
////								return 1;
////							}else if( anno1.type.startsWith("ION") && anno2.type.startsWith("ION") ){
////								if( (anno1.type.compareTo("ION|Y1")==0 || anno1.type.compareTo("ION|B1")==0) && (anno2.type.compareTo("ION|Y1")==0 || anno2.type.compareTo("ION|B1")==0) ){
////									int ion1 = Integer.parseInt(anno1.anno.substring(1));
////									int ion2 = Integer.parseInt(anno2.anno.substring(1));
////									char type1 =  anno1.anno.charAt(0);
////									char type2 =  anno2.anno.charAt(0);
////									
////									if( type1 < type2 ){
////										return -1;
////									}else if( type1 > type2 ){
////										return 1;
////									}else{
////										if( ion1 < ion2 ){
////											return -1;
////										}else if( ion1 > ion2 ){
////											return 1;
////										}else{
////											return 0;
////										}
////									}
////								}else if( (anno1.type.compareTo("ION|Y1")==0 || anno1.type.compareTo("ION|B1")==0) ){
////									return -1;
////								}else if( (anno2.type.compareTo("ION|Y1")==0 || anno2.type.compareTo("ION|B1")==0) ){
////									return 1;
////								}else {
////									char type1 =  anno1.anno.charAt(0);
////									char type2 =  anno2.anno.charAt(0);
////									
////									if( type1 < type2 ){
////										return -1;
////									}else if( type1 > type2 ){
////										return 1;
////									}else{
////										return anno1.anno.compareTo(anno2.anno);
////									}
////								}
////								
////							}else{
////								if( anno1.type.startsWith("LOSS") && !anno2.type.startsWith("LOSS") ){
////									return -1;
////								}else if( !anno1.type.startsWith("LOSS") && anno2.type.startsWith("LOSS") ){
////									return 1;
////								}else{
////									return anno1.anno.compareTo(anno2.anno);
////								}
////							}
////							
////							
////							
////						}
////					});
//
////					for( int j = 0 ; j < annos.size(); j++ ){
////						if( !(annos.get(j).type.compareTo("ION|Y1")==0 || annos.get(j).type.compareTo("ION|B1")==0) )
////							continue;
////						out.append(annos.get(j).anno);
////						out.append("\t");
////						out.append(annos.get(j).mass+"");
////						out.append("\t");
////						out.append(annos.get(j).inten+"");
////						out.append("\n");
////					}
//					double prev_mass = -1;
//					int sum_y1 = 0, cnt_y1 = 0, sum_b1 = 0, cnt_b1 = 0;
//					int sum_y2 = 0, cnt_y2 = 0, sum_b2 = 0, cnt_b2 = 0;
//					int sum_y3 = 0, cnt_y3 = 0, sum_b3 = 0, cnt_b3 = 0;
//					int sum_y4 = 0, cnt_y4 = 0, sum_b4 = 0, cnt_b4 = 0;
//					int sum_ynl = 0, cnt_ynl = 0, sum_bnl = 0, cnt_bnl = 0;
//					for (int j = 0; j < annos.size(); j++) {
//
//						if (annos.get(j).type.compareTo("ION|Y1") == 0) {
//							sum_y1 += annos.get(j).inten;
//							cnt_y1++;
//						} else if (annos.get(j).type.compareTo("ION|B1") == 0) {
//							sum_b1 += annos.get(j).inten;
//							cnt_b1++;
//						} else if (annos.get(j).type.startsWith("ION|Y2") && prev_mass != annos.get(j).mass) {
//							sum_y2 += annos.get(j).inten;
//							cnt_y2++;
//						} else if (annos.get(j).type.startsWith("ION|B2") && prev_mass != annos.get(j).mass) {
//							sum_b2 += annos.get(j).inten;
//							cnt_b2++;
//						} else if (annos.get(j).type.startsWith("ION|Y3") && prev_mass != annos.get(j).mass) {
//							sum_y3 += annos.get(j).inten;
//							cnt_y3++;
//						} else if (annos.get(j).type.startsWith("ION|B3") && prev_mass != annos.get(j).mass) {
//							sum_b3 += annos.get(j).inten;
//							cnt_b3++;
//						} else if (annos.get(j).type.startsWith("ION|Y") && prev_mass != annos.get(j).mass) {
//							sum_y4 += annos.get(j).inten;
//							cnt_y4++;
//						} else if (annos.get(j).type.startsWith("ION|B") && prev_mass != annos.get(j).mass) {
//							sum_b4 += annos.get(j).inten;
//							cnt_b4++;
//						} else if (annos.get(j).type.startsWith("LOSS|Y") && prev_mass != annos.get(j).mass) {
//							sum_ynl += annos.get(j).inten;
//							cnt_ynl++;
//						} else if (annos.get(j).type.startsWith("LOSS|B") && prev_mass != annos.get(j).mass) {
//							sum_bnl += annos.get(j).inten;
//							cnt_bnl++;
//						}
//						prev_mass = annos.get(j).mass;
//					}
////					out.append(((Double)spectrumt.TIC).toString());
////					out.append("\t");
////					out.append(((Double)spectrumt.baseIntensity).toString());
////					out.append("\t");
////					out.append(((Integer)cnt_y1).toString());    
////					out.append("\t");
////					out.append(((Integer)cnt_y2).toString());    
////					out.append("\t");
////					out.append(((Integer)cnt_y3).toString());    
////					out.append("\t");
////					out.append(((Integer)cnt_y4).toString());    
////					out.append("\t");
////					out.append(((Integer)cnt_ynl).toString());    
////					out.append("\t");
////					out.append(((Integer)sum_y1).toString());    
////					out.append("\t");
////					out.append(((Integer)sum_y2).toString());    
////					out.append("\t");
////					out.append(((Integer)sum_y3).toString());    
////					out.append("\t");
////					out.append(((Integer)sum_y4).toString());    
////					out.append("\t");
////					out.append(((Integer)sum_ynl).toString());    
////					out.append("\t");
////					out.append(((Integer)cnt_b1).toString());    
////					out.append("\t");
////					out.append(((Integer)cnt_b2).toString());    
////					out.append("\t");
////					out.append(((Integer)cnt_b3).toString());    
////					out.append("\t");
////					out.append(((Integer)cnt_b4).toString());    
////					out.append("\t");
////					out.append(((Integer)cnt_bnl).toString());    
////					out.append("\t");
////					out.append(((Integer)sum_b1).toString());    
////					out.append("\t");
////					out.append(((Integer)sum_b2).toString());    
////					out.append("\t");
////					out.append(((Integer)sum_b3).toString());    
////					out.append("\t");
////					out.append(((Integer)sum_b4).toString());    
////					out.append("\t");
////					out.append(((Integer)sum_bnl).toString());    
////					out.append("\t");
////					out.append(((Integer)pLen).toString());
////					out.append("\t");
//					out.append(((Double) ((double) cnt_y1 / (pLen - 1))).toString());
//					out.append("\t");
////					out.append(((Double)((double)(cnt_y1+cnt_y2+cnt_y3+cnt_y4+cnt_ynl)/(pLen-1))).toString());
////					out.append("\t");
//					out.append(((Double) ((double) cnt_b1 / (pLen - 1))).toString());
//					out.append("\t");
////					out.append(((Double)((double)(cnt_b1+cnt_b2+cnt_b3+cnt_b4+cnt_bnl)/(pLen-1))).toString());
////					out.append("\t");
//					out.append("\n");
//					sum_y_rate += ((double) cnt_y1 / (pLen - 1));
//					sum_b_rate += ((double) cnt_b1 / (pLen - 1));
//					sum_yall_rate += ((double) (cnt_y1 + cnt_y2 + cnt_y3 + cnt_y4 + cnt_ynl) / (pLen - 1));
//					sum_ball_rate += ((double) (cnt_b1 + cnt_b2 + cnt_b3 + cnt_b4 + cnt_bnl) / (pLen - 1));
//
////					if ((cnt_y1+cnt_y2+cnt_y3+cnt_y4+cnt_ynl+cnt_b1+cnt_b2+cnt_b3 + cnt_b4 + cnt_bnl) == 0) {
////						for( int j = 0 ; j < annos.size(); j++ ){
////							System.out.println(annos.get(j).anno + "\t" + annos.get(j).mass + "\t" + prob);
////						}
//// 					}
//				}
//			}
//			out.close();
//			err.close();
//			System.out.println((double) (sum_y_rate / psm_cnt));
//			System.out.println((double) (sum_b_rate / psm_cnt));
//			System.out.println((double) (sum_yall_rate / psm_cnt));
//			System.out.println((double) (sum_ball_rate / psm_cnt));
//		} catch (Exception e) {
//			System.out.println(1234);
//			e.printStackTrace();
//			System.err.println(e.getMessage());
//		}
//	}
//
//	private List<String> sortMODplusFile(String result) throws IOException {
//		BufferedReader br = new BufferedReader(new FileReader(new File(result)));
//		String line;
//
//		List<String> newLines = new ArrayList<String>();
//
//		while ((line = br.readLine()) != null) {
//			if (line.startsWith("Spec")) {
//				continue;
//			}
//			newLines.add(line);
//		}
//		br.close();
//
//		Collections.sort(newLines, new Comparator<String>() {
//
//			@Override
//			public int compare(String o1, String o2) {
//				String[] sp1 = o1.split("\\s+");
//				String[] sp2 = o2.split("\\s+");
//
//				String frac1 = sp1[0].substring(0, sp1[0].indexOf('.'));
//				String frac2 = sp2[0].substring(0, sp2[0].indexOf('.'));
//
//				if (frac1.compareTo(frac2) < 0) {
//					return -1;
//				} else if (frac1.compareTo(frac2) > 0) {
//					return 1;
//				}
//
//				int ind1 = Integer.parseInt(sp1[1]);
//				int ind2 = Integer.parseInt(sp2[1]);
//
//				if (ind1 < ind2)
//					return -1;
//				else if (ind1 > ind2) {
//					return 1;
//				}
//
//				int cs1 = Integer.parseInt(sp1[3]);
//				int cs2 = Integer.parseInt(sp2[3]);
//
//				if (cs1 < cs2)
//					return -1;
//				else if (cs1 > cs2) {
//					return 1;
//				}
//
//				return 0;
//			}
//		});
//
//		return newLines;
//	}

	public Vector<AnnoPeakT> getAnnotatedPeaks(SpectrumT spectrumT, PeptideT peptide, SearchSummaryT ss) {
		int i, matIndex = -1;

		Vector<AnnoPeakT> annoPeaks = new Vector<AnnoPeakT>();

		// 占싱뤄옙占쏙옙占쏙옙占쏙옙 占쏙옙占쏙옙 y-ion占쏙옙 占쏙옙占쏙옙占쏙옙占쏙옙 찾占승댐옙.
		
		for (i = 0; i < peptide.theo_y_ions.length; i++) {

			// y-ion占쏙옙 m/z占쏙옙占쏙옙 占쏙옙치(tolerance 占쏙옙 占싱삼옙 占쏙옙占싱놂옙占쏙옙 占십댐옙)占싹댐옙 peak占쏙옙 찾占쏙옙 占쏙옙占쏙옙占싼댐옙.
			matIndex = getMatchedPeak(spectrumT, peptide.theo_y_ions[i], 0.00);

			if (matIndex > -1) {
				AnnoPeakT annoPeak = new AnnoPeakT();
				annoPeak.type = "ION|Y1";
				annoPeak.mass = spectrumT.peakMassList.get(matIndex);
				annoPeak.inten = spectrumT.peakIntensityList.get(matIndex);
				annoPeak.priority = 1;
				annoPeak.anno = "y" + (peptide.sequence.length() - i);
				annoPeaks.add(annoPeak);

				double mze = peptide.theo_y_ions[i] - annoPeak.mass;
				annoPeak.abmz = Math.abs(mze);

				peptide.sum_mxe += mze;
				peptide.sum_abmxe += Math.abs(mze);
				peptide.cnt_mze++;
				if (mze > peptide.max_mxe)
					peptide.max_mxe = mze;
				if (Math.abs(mze) > peptide.max_abmxe)
					peptide.max_abmxe = Math.abs(mze);
				if (mze < peptide.min_mxe)
					peptide.min_mxe = mze;
				if (Math.abs(mze) < peptide.min_abmxe)
					peptide.min_abmxe = Math.abs(mze);

			}

			// charge state占쏙옙 2 占싱삼옙占쏙옙 占쏙옙占쏙옙占� m/z占쏙옙占쏙옙 占쏙옙占쏙옙臼占� match 占쏙옙占쏙옙占쏙옙 확占쏙옙占싹울옙 占쏙옙치占싹댐옙 占쏙옙占�
			// 占쏙옙占쏙옙占싼댐옙.
			for (int tcs = 2; tcs <= spectrumT.cs; tcs++) {
				double cutOff = (tcs == spectrumT.cs) ? 0.5 : 0.000;
				matIndex = getMatchedPeak(spectrumT, (peptide.theo_y_ions[i] + PROTON * (tcs - 1)) / tcs, cutOff);
				if (matIndex > -1) {
					AnnoPeakT annoPeak = new AnnoPeakT();
					annoPeak.type = "ION|Y" + tcs;
					annoPeak.mass = spectrumT.peakMassList.get(matIndex);
					annoPeak.inten = spectrumT.peakIntensityList.get(matIndex);
					annoPeak.priority = 3;
					if (tcs == 2)
						annoPeak.anno = "y" + (peptide.sequence.length() - i) + "++";
					// annotation 표占썩에 charge state占쏙옙 표占쏙옙占싼댐옙.
					else
						annoPeak.anno = "y" + (peptide.sequence.length() - i) + "(" + tcs + "+)";
					annoPeaks.add(annoPeak);

					double mze = (peptide.theo_y_ions[i] + PROTON * (tcs - 1)) / tcs - annoPeak.mass;
					peptide.sum_mxe += mze;
					peptide.sum_abmxe += Math.abs(mze);
					peptide.cnt_mze++;
					if (mze > peptide.max_mxe)
						peptide.max_mxe = mze;
					if (Math.abs(mze) > peptide.max_abmxe)
						peptide.max_abmxe = Math.abs(mze);
					if (mze < peptide.min_mxe)
						peptide.min_mxe = mze;
					if (Math.abs(mze) < peptide.min_abmxe)
						peptide.min_abmxe = Math.abs(mze);

				}
			}

			// neutral loss占쏙옙 占싹어났占쏙옙 占쏙옙占쏙옙占� m/z占쏙옙占쏙옙 占쏙옙占쏙옙臼占� match 占쏙옙占쏙옙占쏙옙 확占쏙옙占싼댐옙.
			AnnoPeakT LOSS, NH3;
			LOSS = new AnnoPeakT();
			NH3 = new AnnoPeakT();
			LOSS.inten = 0;
			NH3.inten = 0;
			matIndex = getMatchedPeak(spectrumT, peptide.theo_y_ions[i] - H_2_O, .1);
			// System.out.println(matIndex + " " + (peptide.theo_y_ions[i]-H_2_O));
			if (matIndex > -1) {
				LOSS.type = "LOSS|Y";
				LOSS.mass = spectrumT.peakMassList.get(matIndex);
				LOSS.inten = spectrumT.peakIntensityList.get(matIndex);
				LOSS.anno = "y" + (peptide.sequence.length() - i) + "-H2O";
			}
			matIndex = getMatchedPeak(spectrumT, peptide.theo_y_ions[i] - N_H_3, .1);
			if (matIndex > -1) {
				NH3.type = "LOSS|Y";
				NH3.mass = spectrumT.peakMassList.get(matIndex);
				NH3.inten = spectrumT.peakIntensityList.get(matIndex);
				NH3.anno = "y" + (peptide.sequence.length() - i) + "-NH3";
			}
			if (LOSS.inten < NH3.inten) {
				LOSS.type = NH3.type;
				LOSS.mass = NH3.mass;
				LOSS.inten = NH3.inten;
				LOSS.anno = NH3.anno;
			}
			if (LOSS.inten != 0) {
				AnnoPeakT annoPeak = new AnnoPeakT();
				annoPeak.type = LOSS.type;
				annoPeak.mass = LOSS.mass;
				annoPeak.inten = LOSS.inten;
				annoPeak.priority = 6;
				annoPeak.anno = LOSS.anno;
				annoPeaks.add(annoPeak);

				double mze = 0.0;
				if (annoPeak.anno.endsWith("-NH3"))
					mze = peptide.theo_y_ions[i] - N_H_3 - LOSS.mass;
				else
					mze = peptide.theo_y_ions[i] - H_2_O - LOSS.mass;
				peptide.sum_mxe += mze;
				peptide.sum_abmxe += Math.abs(mze);
				peptide.cnt_mze++;
				if (mze > peptide.max_mxe)
					peptide.max_mxe = mze;
				if (Math.abs(mze) > peptide.max_abmxe)
					peptide.max_abmxe = Math.abs(mze);
				if (mze < peptide.min_mxe)
					peptide.min_mxe = mze;
				if (Math.abs(mze) < peptide.min_abmxe)
					peptide.min_abmxe = Math.abs(mze);
			}

		}

		// System.out.println();

		// 占싱뤄옙占쏙옙占쏙옙占쏙옙 占쏙옙占쏙옙 b-ion占쏙옙 占쏙옙占쏙옙占쏙옙占쏙옙 찾占승댐옙.
		for (i = 0; i < peptide.theo_b_ions.length; i++) {

			// b-ion占쏙옙 m/z占쏙옙占쏙옙 占쏙옙치(tolerance 占쏙옙 占싱삼옙 占쏙옙占싱놂옙占쏙옙 占십댐옙)占싹댐옙 peak占쏙옙 찾占쏙옙 占쏙옙占쏙옙占싼댐옙.
			matIndex = getMatchedPeak(spectrumT, peptide.theo_b_ions[i], 0.00);
			if (matIndex > -1) {
				AnnoPeakT annoPeak = new AnnoPeakT();
				annoPeak.type = "ION|B1";
				annoPeak.mass = spectrumT.peakMassList.get(matIndex);
				annoPeak.inten = spectrumT.peakIntensityList.get(matIndex);
				annoPeak.priority = 2;
				annoPeak.anno = "b" + (i + 1);
				annoPeaks.add(annoPeak);

				double mze = peptide.theo_b_ions[i] - annoPeak.mass;
				peptide.sum_mxe += mze;
				peptide.sum_abmxe += Math.abs(mze);
				peptide.cnt_mze++;
				if (mze > peptide.max_mxe)
					peptide.max_mxe = mze;
				if (Math.abs(mze) > peptide.max_abmxe)
					peptide.max_abmxe = Math.abs(mze);
				if (mze < peptide.min_mxe)
					peptide.min_mxe = mze;
				if (Math.abs(mze) < peptide.min_abmxe)
					peptide.min_abmxe = Math.abs(mze);
			}

			// charge state占쏙옙 2 占싱삼옙占쏙옙 占쏙옙占쏙옙占� m/z占쏙옙占쏙옙 占쏙옙占쏙옙臼占� match 占쏙옙占쏙옙占쏙옙 확占쏙옙占싹울옙 占쏙옙치占싹댐옙 占쏙옙占�
			// 占쏙옙占쏙옙占싼댐옙.
			for (int tcs = 2; tcs <= spectrumT.cs; tcs++) {
				double cutOff = (tcs == spectrumT.cs) ? .5 : .000;

				matIndex = getMatchedPeak(spectrumT, (peptide.theo_b_ions[i] + PROTON * (tcs - 1)) / tcs, cutOff);
				if (matIndex > -1) {
					AnnoPeakT annoPeak = new AnnoPeakT();
					annoPeak.type = "ION|B" + tcs;
					annoPeak.mass = spectrumT.peakMassList.get(matIndex);
					annoPeak.inten = spectrumT.peakIntensityList.get(matIndex);
					annoPeak.priority = 3;
					if (tcs == 2)
						annoPeak.anno = "b" + (i + 1) + "++";
					else
						annoPeak.anno = "b" + (i + 1) + "(" + tcs + "+)";
					annoPeaks.add(annoPeak);

					double mze = (peptide.theo_b_ions[i] + PROTON * (tcs - 1)) / tcs - annoPeak.mass;
					peptide.sum_mxe += mze;
					peptide.sum_abmxe += Math.abs(mze);
					peptide.cnt_mze++;
					if (mze > peptide.max_mxe)
						peptide.max_mxe = mze;
					if (Math.abs(mze) > peptide.max_abmxe)
						peptide.max_abmxe = Math.abs(mze);
					if (mze < peptide.min_mxe)
						peptide.min_mxe = mze;
					if (Math.abs(mze) < peptide.min_abmxe)
						peptide.min_abmxe = Math.abs(mze);
				}
			}
			// neutral loss占쏙옙 占싹어났占쏙옙 占쏙옙占쏙옙占� m/z占쏙옙占쏙옙 占쏙옙占쏙옙臼占� match 占쏙옙占쏙옙占쏙옙 확占쏙옙占싼댐옙.
			AnnoPeakT LOSS, NH3;
			LOSS = new AnnoPeakT();
			NH3 = new AnnoPeakT();
			LOSS.inten = 0;
			NH3.inten = 0;
			matIndex = getMatchedPeak(spectrumT, peptide.theo_b_ions[i] - H_2_O, .1);
			if (matIndex > -1) {
				LOSS.type = "LOSS|B";
				LOSS.mass = spectrumT.peakMassList.get(matIndex);
				LOSS.inten = spectrumT.peakIntensityList.get(matIndex);
				LOSS.anno = "b" + (i + 1) + "-H2O";
			}
			matIndex = getMatchedPeak(spectrumT, peptide.theo_b_ions[i] - N_H_3, .1);
			if (matIndex > -1) {
				NH3.type = "LOSS|B";
				NH3.mass = spectrumT.peakMassList.get(matIndex);
				NH3.inten = spectrumT.peakIntensityList.get(matIndex);
				NH3.anno = "b" + (i + 1) + "-NH3";
			}
			if (LOSS.inten < NH3.inten) {
				LOSS.type = NH3.type;
				LOSS.mass = NH3.mass;
				LOSS.inten = NH3.inten;
				LOSS.anno = NH3.anno;
			}
			if (LOSS.inten != 0) {
				AnnoPeakT annoPeak = new AnnoPeakT();
				annoPeak.type = LOSS.type;
				annoPeak.mass = LOSS.mass;
				annoPeak.inten = LOSS.inten;
				annoPeak.priority = 6;
				annoPeak.anno = LOSS.anno;
				annoPeaks.add(annoPeak);

				double mze = 0.0;
				if (annoPeak.anno.endsWith("-NH3"))
					mze = peptide.theo_b_ions[i] - N_H_3 - LOSS.mass;
				else
					mze = peptide.theo_b_ions[i] - H_2_O - LOSS.mass;
				peptide.sum_mxe += mze;
				peptide.sum_abmxe += Math.abs(mze);
				peptide.cnt_mze++;
				if (mze > peptide.max_mxe)
					peptide.max_mxe = mze;
				if (Math.abs(mze) > peptide.max_abmxe)
					peptide.max_abmxe = Math.abs(mze);
				if (mze < peptide.min_mxe)
					peptide.min_mxe = mze;
				if (Math.abs(mze) < peptide.min_abmxe)
					peptide.min_abmxe = Math.abs(mze);
			}
		}

		// a-ion占쏙옙 占쏙옙占쏙옙占쏙옙 확占쏙옙占싼댐옙.
		for (i = 1; i < 4 && i < peptide.theo_b_ions.length; i++) {
			matIndex = getMatchedPeak(spectrumT, peptide.theo_b_ions[i] - 27.9949, .1);
			if (matIndex > -1) {
				AnnoPeakT annoPeak = new AnnoPeakT();
				annoPeak.type = "ION|A";
				annoPeak.mass = spectrumT.peakMassList.get(matIndex);
				annoPeak.inten = spectrumT.peakIntensityList.get(matIndex);
				annoPeak.priority = 4;
				annoPeak.anno = "a" + (i + 1);
				annoPeaks.add(annoPeak);

				double mze = peptide.theo_b_ions[i] - 27.9949 - annoPeak.mass;
				peptide.sum_mxe += mze;
				peptide.sum_abmxe += Math.abs(mze);
				peptide.cnt_mze++;
				if (mze > peptide.max_mxe)
					peptide.max_mxe = mze;
				if (Math.abs(mze) > peptide.max_abmxe)
					peptide.max_abmxe = Math.abs(mze);
				if (mze < peptide.min_mxe)
					peptide.min_mxe = mze;
				if (Math.abs(mze) < peptide.min_abmxe)
					peptide.min_abmxe = Math.abs(mze);
			}
		}

		// Immonium ion占쏙옙 占쏙옙占쏙옙占쏙옙 확占쏙옙
		HashSet<Character> imm = new HashSet<Character>();
		for (i = 0; i < peptide.sequence.length(); i++) {
			if (imm.contains(peptide.sequence.charAt(i)))
				continue;
			matIndex = getMatchedPeak(spectrumT,
					ss.getAAMass(peptide.sequence.charAt(i)) + peptide.modified_position[i] - 27, .1);
			if (matIndex > -1) {
				imm.add(peptide.sequence.charAt(i));
				AnnoPeakT annoPeak = new AnnoPeakT();
				annoPeak.type = "IMM|B";
				annoPeak.mass = spectrumT.peakMassList.get(matIndex);
				annoPeak.inten = spectrumT.peakIntensityList.get(matIndex);
				annoPeak.priority = 6;
				if (peptide.modified_position[i] == 0)
					annoPeak.anno = peptide.sequence.charAt(i) + "*";
				else if (peptide.modified_position[i] > 0)
					annoPeak.anno = peptide.sequence.charAt(i) + "+" + Math.round(peptide.modified_position[i]) + "*";
				else
					annoPeak.anno = "" + peptide.sequence.charAt(i) + Math.round(peptide.modified_position[i]) + "*";
				annoPeaks.add(annoPeak);

				double mze = ss.getAAMass(peptide.sequence.charAt(i)) + peptide.modified_position[i] - 27
						- annoPeak.mass;
				peptide.sum_mxe += mze;
				peptide.sum_abmxe += Math.abs(mze);
				peptide.cnt_mze++;
				if (mze > peptide.max_mxe)
					peptide.max_mxe = mze;
				if (Math.abs(mze) > peptide.max_abmxe)
					peptide.max_abmxe = Math.abs(mze);
				if (mze < peptide.min_mxe)
					peptide.min_mxe = mze;
				if (Math.abs(mze) < peptide.min_abmxe)
					peptide.min_abmxe = Math.abs(mze);
			}
		}
		if (annoPeaks.size() > 0)
			annoPeaks = qSort(annoPeaks);
		return annoPeaks;
	}
	public Vector<AnnoPeakT> getAnnotatedPeaks_Phopho(SpectrumT spectrumT, PeptideT peptide, SearchSummaryT ss) {
		int i, matIndex = -1;

		Vector<AnnoPeakT> annoPeaks = new Vector<AnnoPeakT>();

		boolean phospho_flag=false;
		
		// 占싱뤄옙占쏙옙占쏙옙占쏙옙 占쏙옙占쏙옙 y-ion占쏙옙 占쏙옙占쏙옙占쏙옙占쏙옙 찾占승댐옙.
		for (i = 0; i < peptide.theo_y_ions.length; i++) {

			// y-ion占쏙옙 m/z占쏙옙占쏙옙 占쏙옙치(tolerance 占쏙옙 占싱삼옙 占쏙옙占싱놂옙占쏙옙 占십댐옙)占싹댐옙 peak占쏙옙 찾占쏙옙 占쏙옙占쏙옙占싼댐옙.
			matIndex = getMatchedPeak(spectrumT, peptide.theo_y_ions[i], 0.00);

			if (matIndex > -1) {
				AnnoPeakT annoPeak = new AnnoPeakT();
				annoPeak.type = "ION|Y1";
				annoPeak.mass = spectrumT.peakMassList.get(matIndex);
				annoPeak.inten = spectrumT.peakIntensityList.get(matIndex);
				annoPeak.priority = 1;
				annoPeak.anno = "y" + (peptide.sequence.length() - i);
				annoPeaks.add(annoPeak);

				double mze = peptide.theo_y_ions[i] - annoPeak.mass;
				annoPeak.abmz = Math.abs(mze);

				peptide.sum_mxe += mze;
				peptide.sum_abmxe += Math.abs(mze);
				peptide.cnt_mze++;
				if (mze > peptide.max_mxe)
					peptide.max_mxe = mze;
				if (Math.abs(mze) > peptide.max_abmxe)
					peptide.max_abmxe = Math.abs(mze);
				if (mze < peptide.min_mxe)
					peptide.min_mxe = mze;
				if (Math.abs(mze) < peptide.min_abmxe)
					peptide.min_abmxe = Math.abs(mze);

			}

			// charge state占쏙옙 2 占싱삼옙占쏙옙 占쏙옙占쏙옙占� m/z占쏙옙占쏙옙 占쏙옙占쏙옙臼占� match 占쏙옙占쏙옙占쏙옙 확占쏙옙占싹울옙 占쏙옙치占싹댐옙 占쏙옙占�
			// 占쏙옙占쏙옙占싼댐옙.
			for (int tcs = 2; tcs <= spectrumT.cs; tcs++) {
				double cutOff = (tcs == spectrumT.cs) ? 0.5 : 0.000;
				
				matIndex = getMatchedPeak(spectrumT, (peptide.theo_y_ions[i] + PROTON * (tcs - 1)) / tcs, cutOff);
				
				if (matIndex > -1) {
					AnnoPeakT annoPeak = new AnnoPeakT();
					annoPeak.type = "ION|Y" + tcs;
					annoPeak.mass = spectrumT.peakMassList.get(matIndex);
					annoPeak.inten = spectrumT.peakIntensityList.get(matIndex);
					annoPeak.priority = 3;
					if (tcs == 2)
						annoPeak.anno = "y" + (peptide.sequence.length() - i) + "++";
					// annotation 표占썩에 charge state占쏙옙 표占쏙옙占싼댐옙.
					else
						annoPeak.anno = "y" + (peptide.sequence.length() - i) + "(" + tcs + "+)";
					annoPeaks.add(annoPeak);

					double mze = (peptide.theo_y_ions[i] + PROTON * (tcs - 1)) / tcs - annoPeak.mass;
					peptide.sum_mxe += mze;
					peptide.sum_abmxe += Math.abs(mze);
					peptide.cnt_mze++;
					if (mze > peptide.max_mxe)
						peptide.max_mxe = mze;
					if (Math.abs(mze) > peptide.max_abmxe)
						peptide.max_abmxe = Math.abs(mze);
					if (mze < peptide.min_mxe)
						peptide.min_mxe = mze;
					if (Math.abs(mze) < peptide.min_abmxe)
						peptide.min_abmxe = Math.abs(mze);

				}
				
			}
			///////////////////////////////// nh3. h20
			
			for (int tcs = 1; tcs <= spectrumT.cs-1; tcs++) {
				double cutOff = (tcs == spectrumT.cs) ? 0.5 : 0.000;
				
				AnnoPeakT LOSS, NH3;
				LOSS = new AnnoPeakT();
				NH3 = new AnnoPeakT();
				LOSS.inten = 0;
				NH3.inten = 0;
				
				if(tcs==1)
					cutOff=.1;
				
				matIndex = getMatchedPeak(spectrumT, ((peptide.theo_y_ions[i]-N_H_3) + PROTON * (tcs - 1)) / tcs, cutOff);
				
				if (matIndex > -1) {
					NH3.type = "LOSS|Y";
					NH3.mass = spectrumT.peakMassList.get(matIndex);
					NH3.inten = spectrumT.peakIntensityList.get(matIndex);
					if (tcs == 2)
						NH3.anno = "y" + (peptide.sequence.length() - i) + "-NH3++";
					// annotation 표占썩에 charge state占쏙옙 표占쏙옙占싼댐옙.
					else
						NH3.anno = "y" + (peptide.sequence.length() - i) + "-NH3(" + tcs + "+)";
					
				}
				
				matIndex = getMatchedPeak(spectrumT, ((peptide.theo_y_ions[i]-H_2_O) + PROTON * (tcs - 1)) / tcs, cutOff);
				
				// System.out.println(matIndex + " " + (peptide.theo_y_ions[i]-H_2_O));
				if (matIndex > -1) {
					LOSS.type = "LOSS|Y";
					LOSS.mass = spectrumT.peakMassList.get(matIndex);
					LOSS.inten = spectrumT.peakIntensityList.get(matIndex);
					if (tcs == 2)
						LOSS.anno = "y" + (peptide.sequence.length() - i) + "-H20++";
					// annotation 표占썩에 charge state占쏙옙 표占쏙옙占싼댐옙.
					else
						LOSS.anno = "y" + (peptide.sequence.length() - i) + "-H2O(" + tcs + "+)";
			      }
				if (LOSS.inten < NH3.inten) {
					LOSS.type = NH3.type;
					LOSS.mass = NH3.mass;
					LOSS.inten = NH3.inten;
					LOSS.anno = NH3.anno;
				}
				if (LOSS.inten != 0) {
					AnnoPeakT annoPeak = new AnnoPeakT();
					annoPeak.type = LOSS.type;
					annoPeak.mass = LOSS.mass;
					annoPeak.inten = LOSS.inten;
					annoPeak.priority = 6;
					annoPeak.anno = LOSS.anno;
					annoPeaks.add(annoPeak);

					double mze = 0.0;
					if (annoPeak.anno.endsWith("-NH3"))
						mze = peptide.theo_y_ions[i] - N_H_3 - LOSS.mass;
					else if(annoPeak.anno.endsWith("-H2O"))
						mze = peptide.theo_y_ions[i] - H_2_O - LOSS.mass;
					
					
					peptide.sum_mxe += mze;
					peptide.sum_abmxe += Math.abs(mze);
					peptide.cnt_mze++;
					if (mze > peptide.max_mxe)
						peptide.max_mxe = mze;
					if (Math.abs(mze) > peptide.max_abmxe)
						peptide.max_abmxe = Math.abs(mze);
					if (mze < peptide.min_mxe)
						peptide.min_mxe = mze;
					if (Math.abs(mze) < peptide.min_abmxe)
						peptide.min_abmxe = Math.abs(mze);
				}
				
			}
			
				///////////////////////////////// pnh3. h20
			

			
			// neutral loss占쏙옙 占싹어났占쏙옙 占쏙옙占쏙옙占� m/z占쏙옙占쏙옙 占쏙옙占쏙옙臼占� match 占쏙옙占쏙옙占쏙옙 확占쏙옙占싼댐옙.
			
			if(this.phosphodata==true)
			{
		
				if(phospho_ioncheck(i, peptide.phospho_chk)==true)
				{
					
					for (int tcs = 1; tcs <= spectrumT.cs-1; tcs++) {
						double cutOff = (tcs == spectrumT.cs) ? 0.5 : 0.000;
						
						AnnoPeakT pLOSS;
						pLOSS=new AnnoPeakT();
						pLOSS.inten=0;
						
						
						if(tcs==1)
							cutOff=.1;
						
						matIndex = getMatchedPeak(spectrumT, ((peptide.theo_y_ions[i]-phospho_neutral_loss) + PROTON * (tcs - 1)) / tcs, cutOff);
						
						if (matIndex > -1) {
							pLOSS.type = "pLOSS|Y";
							pLOSS.mass = spectrumT.peakMassList.get(matIndex);
							pLOSS.inten = spectrumT.peakIntensityList.get(matIndex);
							if (tcs == 2)
								pLOSS.anno = "y" + (peptide.sequence.length() - i) + "-pH20++";
							// annotation 표占썩에 charge state占쏙옙 표占쏙옙占싼댐옙.
							else
								pLOSS.anno = "y" + (peptide.sequence.length() - i) + "-pH2O(" + tcs + "+)";
							
						}
						
						if (pLOSS.inten != 0) {
							AnnoPeakT annoPeak = new AnnoPeakT();
							annoPeak.type = pLOSS.type;
							annoPeak.mass = pLOSS.mass;
							annoPeak.inten = pLOSS.inten;
							annoPeak.priority = 6;
							annoPeak.anno = pLOSS.anno;
							annoPeaks.add(annoPeak);
		
							double mze = peptide.theo_y_ions[i]-phospho_neutral_loss -pLOSS.mass;
													
							peptide.sum_mxe += mze;
							peptide.sum_abmxe += Math.abs(mze);
							peptide.cnt_mze++;
							if (mze > peptide.max_mxe)
								peptide.max_mxe = mze;
							if (Math.abs(mze) > peptide.max_abmxe)
								peptide.max_abmxe = Math.abs(mze);
							if (mze < peptide.min_mxe)
								peptide.min_mxe = mze;
							if (Math.abs(mze) < peptide.min_abmxe)
								peptide.min_abmxe = Math.abs(mze);		
					}
						
					}
			
				}
				
			}

		}
		
		phospho_flag=false;

		// System.out.println();

		// 占싱뤄옙占쏙옙占쏙옙占쏙옙 占쏙옙占쏙옙 b-ion占쏙옙 占쏙옙占쏙옙占쏙옙占쏙옙 찾占승댐옙.
		for (i = 0; i < peptide.theo_b_ions.length; i++) {

			// b-ion占쏙옙 m/z占쏙옙占쏙옙 占쏙옙치(tolerance 占쏙옙 占싱삼옙 占쏙옙占싱놂옙占쏙옙 占십댐옙)占싹댐옙 peak占쏙옙 찾占쏙옙 占쏙옙占쏙옙占싼댐옙.
			matIndex = getMatchedPeak(spectrumT, peptide.theo_b_ions[i], 0.00);
			if (matIndex > -1) {
				AnnoPeakT annoPeak = new AnnoPeakT();
				annoPeak.type = "ION|B1";
				annoPeak.mass = spectrumT.peakMassList.get(matIndex);
				annoPeak.inten = spectrumT.peakIntensityList.get(matIndex);
				annoPeak.priority = 2;
				annoPeak.anno = "b" + (i + 1);
				annoPeaks.add(annoPeak);

				double mze = peptide.theo_b_ions[i] - annoPeak.mass;
				peptide.sum_mxe += mze;
				peptide.sum_abmxe += Math.abs(mze);
				peptide.cnt_mze++;
				if (mze > peptide.max_mxe)
					peptide.max_mxe = mze;
				if (Math.abs(mze) > peptide.max_abmxe)
					peptide.max_abmxe = Math.abs(mze);
				if (mze < peptide.min_mxe)
					peptide.min_mxe = mze;
				if (Math.abs(mze) < peptide.min_abmxe)
					peptide.min_abmxe = Math.abs(mze);
			}

			// charge state占쏙옙 2 占싱삼옙占쏙옙 占쏙옙占쏙옙占� m/z占쏙옙占쏙옙 占쏙옙占쏙옙臼占� match 占쏙옙占쏙옙占쏙옙 확占쏙옙占싹울옙 占쏙옙치占싹댐옙 占쏙옙占�
			// 占쏙옙占쏙옙占싼댐옙.
			for (int tcs = 2; tcs <= spectrumT.cs; tcs++) {
				double cutOff = (tcs == spectrumT.cs) ? .5 : .000;

				matIndex = getMatchedPeak(spectrumT, (peptide.theo_b_ions[i] + PROTON * (tcs - 1)) / tcs, cutOff);
				if (matIndex > -1) {
					AnnoPeakT annoPeak = new AnnoPeakT();
					annoPeak.type = "ION|B" + tcs;
					annoPeak.mass = spectrumT.peakMassList.get(matIndex);
					annoPeak.inten = spectrumT.peakIntensityList.get(matIndex);
					annoPeak.priority = 3;
					if (tcs == 2)
						annoPeak.anno = "b" + (i + 1) + "++";
					else
						annoPeak.anno = "b" + (i + 1) + "(" + tcs + "+)";
					annoPeaks.add(annoPeak);

					double mze = (peptide.theo_b_ions[i] + PROTON * (tcs - 1)) / tcs - annoPeak.mass;
					peptide.sum_mxe += mze;
					peptide.sum_abmxe += Math.abs(mze);
					peptide.cnt_mze++;
					if (mze > peptide.max_mxe)
						peptide.max_mxe = mze;
					if (Math.abs(mze) > peptide.max_abmxe)
						peptide.max_abmxe = Math.abs(mze);
					if (mze < peptide.min_mxe)
						peptide.min_mxe = mze;
					if (Math.abs(mze) < peptide.min_abmxe)
						peptide.min_abmxe = Math.abs(mze);
				}
			}
			// neutral loss占쏙옙 占싹어났占쏙옙 占쏙옙占쏙옙占� m/z占쏙옙占쏙옙 占쏙옙占쏙옙臼占� match 占쏙옙占쏙옙占쏙옙 확占쏙옙占싼댐옙.
			
			//nh3, h2o
			for (int tcs = 1; tcs <= spectrumT.cs-1; tcs++) {
				double cutOff = (tcs == spectrumT.cs) ? .5 : .000;

				if(tcs==1)
					cutOff=.1;
				
				AnnoPeakT LOSS, NH3;
				LOSS = new AnnoPeakT();
				NH3 = new AnnoPeakT();
				LOSS.inten = 0;
				NH3.inten = 0;
				
				matIndex = getMatchedPeak(spectrumT, ((peptide.theo_b_ions[i]-N_H_3) + PROTON * (tcs - 1)) / tcs, cutOff);
			
				if (matIndex > -1) {
					NH3.type = "LOSS|B";
					NH3.mass = spectrumT.peakMassList.get(matIndex);
					NH3.inten = spectrumT.peakIntensityList.get(matIndex);
					if (tcs == 2)
						NH3.anno = "b" + (i + 1) + "-NH3++";
					else
						NH3.anno = "b" + (i + 1) + "-NH3(" + tcs + "+)";
				}
				
				matIndex = getMatchedPeak(spectrumT, ((peptide.theo_b_ions[i]-H_2_O) + PROTON * (tcs - 1)) / tcs, cutOff);
				
				if (matIndex > -1) {
						LOSS.type = "LOSS|B";
						LOSS.mass = spectrumT.peakMassList.get(matIndex);
						LOSS.inten = spectrumT.peakIntensityList.get(matIndex);
						if (tcs == 2)
							LOSS.anno = "b" + (i + 1) + "-H2O++";
						else
							LOSS.anno = "b" + (i + 1) + "-H2O(" + tcs + "+)";
					}
				
				
				if (LOSS.inten < NH3.inten) {
					LOSS.type = NH3.type;
					LOSS.mass = NH3.mass;
					LOSS.inten = NH3.inten;
					LOSS.anno = NH3.anno;
				}
				if (LOSS.inten != 0) {
					AnnoPeakT annoPeak = new AnnoPeakT();
					annoPeak.type = LOSS.type;
					annoPeak.mass = LOSS.mass;
					annoPeak.inten = LOSS.inten;
					annoPeak.priority = 6;
					annoPeak.anno = LOSS.anno;
					annoPeaks.add(annoPeak);

					double mze = 0.0;
					if (annoPeak.anno.endsWith("-NH3"))
						mze = peptide.theo_b_ions[i] - N_H_3 - LOSS.mass;
					else if(annoPeak.anno.endsWith("-H2O"))
						mze = peptide.theo_b_ions[i] - H_2_O - LOSS.mass;
					
					
					
					peptide.sum_mxe += mze;
					peptide.sum_abmxe += Math.abs(mze);
					peptide.cnt_mze++;
					if (mze > peptide.max_mxe)
						peptide.max_mxe = mze;
					if (Math.abs(mze) > peptide.max_abmxe)
						peptide.max_abmxe = Math.abs(mze);
					if (mze < peptide.min_mxe)
						peptide.min_mxe = mze;
					if (Math.abs(mze) < peptide.min_abmxe)
						peptide.min_abmxe = Math.abs(mze);
				}
			}
			
			//ph20
			if(phosphodata=true)
			{
				if(peptide.phospho_chk[i]==true || phospho_flag==true)
				{
					phospho_flag=true;
					
					for (int tcs = 1; tcs <= spectrumT.cs-1; tcs++) {
						double cutOff = (tcs == spectrumT.cs) ? .5 : .000;
		
						AnnoPeakT pLOSS;
						pLOSS=new AnnoPeakT();
						pLOSS.inten=0;
						
						if(tcs==1)
							cutOff=.1;
						
						matIndex = getMatchedPeak(spectrumT, ((peptide.theo_b_ions[i]-phospho_neutral_loss) + PROTON * (tcs - 1)) / tcs, cutOff);
						
						if (matIndex > -1) {
							pLOSS.type = "pLOSS|B";
							pLOSS.mass = spectrumT.peakMassList.get(matIndex);
							pLOSS.inten = spectrumT.peakIntensityList.get(matIndex);
							if (tcs == 2)
								pLOSS.anno = "b" + (i + 1) + "-pH2O++";
							else
								pLOSS.anno = "b" + (i + 1) + "-pH2O(" + tcs + "+)";
						}
						if (pLOSS.inten != 0) {
							AnnoPeakT annoPeak = new AnnoPeakT();
							annoPeak.type = pLOSS.type;
							annoPeak.mass = pLOSS.mass;
							annoPeak.inten = pLOSS.inten;
							annoPeak.priority = 6;
							annoPeak.anno = pLOSS.anno;
							annoPeaks.add(annoPeak);
		
							double mze = peptide.theo_b_ions[i] - phospho_neutral_loss - pLOSS.mass;
							
							
							peptide.sum_mxe += mze;
							peptide.sum_abmxe += Math.abs(mze);
							peptide.cnt_mze++;
							if (mze > peptide.max_mxe)
								peptide.max_mxe = mze;
							if (Math.abs(mze) > peptide.max_abmxe)
								peptide.max_abmxe = Math.abs(mze);
							if (mze < peptide.min_mxe)
								peptide.min_mxe = mze;
							if (Math.abs(mze) < peptide.min_abmxe)
								peptide.min_abmxe = Math.abs(mze);
						
					}
					}
				}
				}
			
		}

		// a-ion占쏙옙 占쏙옙占쏙옙占쏙옙 확占쏙옙占싼댐옙.
		for (i = 1; i < 4 && i < peptide.theo_b_ions.length; i++) {
			matIndex = getMatchedPeak(spectrumT, peptide.theo_b_ions[i] - 27.9949, .1);
			if (matIndex > -1) {
				AnnoPeakT annoPeak = new AnnoPeakT();
				annoPeak.type = "ION|A";
				annoPeak.mass = spectrumT.peakMassList.get(matIndex);
				annoPeak.inten = spectrumT.peakIntensityList.get(matIndex);
				annoPeak.priority = 4;
				annoPeak.anno = "a" + (i + 1);
				annoPeaks.add(annoPeak);

				double mze = peptide.theo_b_ions[i] - 27.9949 - annoPeak.mass;
				peptide.sum_mxe += mze;
				peptide.sum_abmxe += Math.abs(mze);
				peptide.cnt_mze++;
				if (mze > peptide.max_mxe)
					peptide.max_mxe = mze;
				if (Math.abs(mze) > peptide.max_abmxe)
					peptide.max_abmxe = Math.abs(mze);
				if (mze < peptide.min_mxe)
					peptide.min_mxe = mze;
				if (Math.abs(mze) < peptide.min_abmxe)
					peptide.min_abmxe = Math.abs(mze);
			}
		}
		// Immonium ion占쏙옙 占쏙옙占쏙옙占쏙옙 확占쏙옙
		HashSet<Character> imm = new HashSet<Character>();
		for (i = 0; i < peptide.sequence.length(); i++) {
			if (imm.contains(peptide.sequence.charAt(i)))
				continue;
			matIndex = getMatchedPeak(spectrumT,
					ss.getAAMass(peptide.sequence.charAt(i)) + peptide.modified_position[i] - 27, .1);
			if (matIndex > -1) {
				imm.add(peptide.sequence.charAt(i));
				AnnoPeakT annoPeak = new AnnoPeakT();
				annoPeak.type = "IMM|B";
				annoPeak.mass = spectrumT.peakMassList.get(matIndex);
				annoPeak.inten = spectrumT.peakIntensityList.get(matIndex);
				annoPeak.priority = 6;
				if (peptide.modified_position[i] == 0)
					annoPeak.anno = peptide.sequence.charAt(i) + "*";
				else if (peptide.modified_position[i] > 0)
					annoPeak.anno = peptide.sequence.charAt(i) + "+" + Math.round(peptide.modified_position[i]) + "*";
				else
					annoPeak.anno = "" + peptide.sequence.charAt(i) + Math.round(peptide.modified_position[i]) + "*";
				annoPeaks.add(annoPeak);

				double mze = ss.getAAMass(peptide.sequence.charAt(i)) + peptide.modified_position[i] - 27
						- annoPeak.mass;
				peptide.sum_mxe += mze;
				peptide.sum_abmxe += Math.abs(mze);
				peptide.cnt_mze++;
				if (mze > peptide.max_mxe)
					peptide.max_mxe = mze;
				if (Math.abs(mze) > peptide.max_abmxe)
					peptide.max_abmxe = Math.abs(mze);
				if (mze < peptide.min_mxe)
					peptide.min_mxe = mze;
				if (Math.abs(mze) < peptide.min_abmxe)
					peptide.min_abmxe = Math.abs(mze);
			}
		}
		if (annoPeaks.size() > 0)
			annoPeaks = qSort(annoPeaks);
		return annoPeaks;
	}

private boolean phospho_ioncheck(int index, boolean[] phospho_chk) {
		boolean flag=false;
	    for(int i=index;i<phospho_chk.length;i++)
	    {
	    	if(phospho_chk[i]==true)
	    		flag=true;
	    }
		return true;
	}
private boolean check_phosphoion(double[] modified_position, int pos, String iontype, boolean phospho_flag) {
		
	
	
	if(iontype.equalsIgnoreCase("y"))
	{
		if(phospho_flag==false)
		{
			if(Math.round(modified_position[modified_position.length-pos+1]*1000)/1000.0==79.966)
			{
				
			}
		}
			
			
	}
	
	
	
		return phospho_flag;
	}
//	Vector 占쏙옙占쏙옙 占쏙옙 AnnoPeak占쌘료구占쏙옙占쏙옙 占쏙옙치占쏙옙 swap占싼댐옙.
	private Vector<AnnoPeakT> swap(Vector<AnnoPeakT> vector, int i, int j) {
		vector.set(j, vector.set(i, vector.get(j)));
		return vector;
	}

	// Annotation占쏙옙 占썼열占쏙옙 mass占쏙옙 priority占쏙옙占쏙옙占쏙옙 sorting占싼댐옙.
	public Vector<AnnoPeakT> qSort(Vector<AnnoPeakT> annoPeaks) {
		if (annoPeaks.size() == 0)
			return new Vector<AnnoPeakT>();
		else if (annoPeaks.size() == 1)
			return annoPeaks;
		else if (annoPeaks.size() == 2) {
			if (compareAnnoPeaks(annoPeaks.get(0), annoPeaks.get(1)) > 0)
				annoPeaks = swap(annoPeaks, 0, 1);
			return annoPeaks;
		}
		AnnoPeakT pivot = annoPeaks.lastElement();
		Vector<AnnoPeakT> lefts = new Vector<AnnoPeakT>();
		Vector<AnnoPeakT> rights = new Vector<AnnoPeakT>();

		for (int i = 0; i < annoPeaks.size() - 1; i++) {
			if (compareAnnoPeaks(annoPeaks.get(i), pivot) >= 0)
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

	// 특占쏙옙 m/z占쏙옙占쏙옙 占쏙옙占쏙옙 peak占쏙옙 spectrum占쏙옙占쏙옙 찾占쏙옙 index占쏙옙 占싼곤옙占쌔댐옙.
	public int getMatchedPeak(SpectrumT spectrumT, double d, double cutOff) {
		int x = rere.spectrum.getMatchedPeak(spectrumT, d, cutOff);
		return x;
	}

	private int compareAnnoPeaks(AnnoPeakT a, AnnoPeakT b) {
		double x = a.mass, y = b.mass;
		int xp = a.priority, yp = b.priority;

		if (x > y)
			return 1;
		else if (x == y) {
			if (xp > yp)
				return 1;
			else if (xp == yp)
				return 0;
			else
				return -1;
		} else
			return -1;
	}

	private int getLongestConsecutive(int[] array, int threshold) {
		int max = 0;
		int prev = -1;
		int len = 0;
		for (int i = 0; i < array.length; i++) {
			if (array[i] >= threshold && prev + 1 == i) {
				len++;
				if (max < len)
					max = len;
				prev = i;
			} else if (array[i] >= threshold) {
				if (max < len)
					max = len;
				len = 1;
			}
		}
		return max;
	}
	
}
