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
import java.util.*;

//import org.apache.commons.jexl2.parser.ParserTokenManager;

import java.io.*;

public class ReadXML {

	ArrayList<String> lists=new ArrayList<String>();
	HashMap<String,ArrayList<String>> pepinfo=new HashMap<String,ArrayList<String>>();
	HashMap<String,Integer> numPSMperPep=new HashMap<String,Integer>();
	HashMap<String,ArrayList<String>> PeplistperProt=new HashMap<String,ArrayList<String>>();
	HashMap<String,String> ProtlistperPep=new HashMap<String,String>();
	String files;
	String Header="";
	
	public ReadXML(String file) throws Exception
	{
		files=file;
		
		XML();
		
	}

	private void XML() throws Exception {
		
		FileReader fr=new FileReader(files);
		
		BufferedReader br=new BufferedReader(fr);
		
		
		System.out.println("read pepXML");
	    String line="";
	    
	    int flag=0;
	    int hflag=0;
	    int hitflag=0;
	    
	    String info="";
	    
	    
	    String hyperscore = "";
	    String nextscore="";
	    
	    
	    String xcorr = "";
	    String deltacn="";
	    String deltacnstar="";
	    String spscore="";
	    String sprank="";
	    String expect="";
	    String modPEP="";
	    String protein="";
	    String OriModPEP="";
	    String peptide="";
	    
	    
	    while((line=br.readLine())!=null)
	    {
	    	if(flag==0)
	    	{
	    		if(line.contains("<spectrum_query"))
	    		{
	    			flag=1;
	    			String Spectrum=this.extract(line,"spectrum");
	    			
	    			String SpectrumNativeID=this.extract(line,"spectrumNativeID");
	    			String Scan=this.extract(line,"start_scan");
	    			String index=this.extract(line,"index");
	    			String Precursor_neutral_mass=this.extract(line,"precursor_neutral_mass");
	    			String Charge=this.extract(line,"assumed_charge");
	    			
	    			//SpectrumNativeID=Spectrum.split("\\.")[0]+"."+Scan+"."+Scan+"."+Charge;
	    			
	    			if(hflag==0)
	    			{
	    				Header="#SpectrumFile\tFraction\tScan\tSpectrum\tSpectrumNativeID\tIndex\tPrecursor_neutral_mass\tCharge";
	    			}
	    			String tmp=Spectrum.split("\\.")[0].split("_")[Spectrum.split("\\.")[0].split("_").length-1];
	    			if(SpectrumNativeID.equals(""))
	    			{
	    				SpectrumNativeID=Spectrum.split("\\.")[0]+"."+Scan+"."+Scan+"."+Charge;
	    			}
	    			info=SpectrumNativeID.split("\\.")[0]+".mgf\tFraction="+tmp+"\t"+Scan+"\t"+Spectrum+"\t"+SpectrumNativeID+"\t"+index+"\t"+Precursor_neutral_mass+"\t"+Charge;
	    			
	    		}
	    		
	    	}
	    	else
	    	{
	    		if(line.contains("hit_rank=\"1\""))
	    		{
	    			hitflag=1;
	    			String pep=this.extract(line, "peptide"); 
	    			
	    			//if(pep.equals("TGAPCR"))
	    			
	    			peptide=pep;
	    			
	    			String pre=this.extract(line, "peptide_prev_aa");
	    			String next=this.extract(line,"peptide_next_aa");
	    			protein=this.extract(line,"protein")+"(pre="+pre+",next="+next+")";
	    			//double numProt=Double.parseDouble(this.extract(line,"num_tot_proteins"));
	    			String calc_neutral_pep_mass=this.extract(line, "calc_neutral_pep_mass");
	    			String MassDiff=this.extract(line, "massdiff");
//	    			String num_tol_term=this.extract(line,"num_tol_term");
	    			String num_missed_cleavage=this.extract(line,"num_missed_cleavages");
	    			String num_matched_ions=this.extract(line,"num_matched_ions");
	    			String tot_num_ions=this.extract(line,"tot_num_ions");
	    			String num_matched_peptides=this.extract(line,"num_matched_peptides");
	    			
	    			
	    			if(!pepinfo.containsKey(pep))
	    			{
	    				ArrayList<String> al=new ArrayList<String>();
	    				al.add(pep);
	    				
	    				pepinfo.put(pep,al);	
	    			}
	 
	    			if(hflag==0)
	    			{
	    				Header+="\tPeptide\tCal_neutral_pep_mass\tMassDiff\tAbsolutMassDiss\tNum_missed_cleavage\tIonFrac\tLnNumPepInDB";
	    			}
	    			
	    			info+="\t"+pep+"\t"+calc_neutral_pep_mass+"\t"+MassDiff+"\t"+Math.abs(Double.parseDouble(MassDiff))+"\t"+num_missed_cleavage+"\t"+(Double.parseDouble(num_matched_ions)/Double.parseDouble(tot_num_ions))+"\t"+Math.log(Double.parseDouble(num_matched_peptides));
	    			modPEP=pep;
	    			//OriModPEP=pep;
	    			peptide=pep;
	    			
	    		}
	    		else if(hitflag==1)
	    		{
	    			if(line.contains("<search_score"))
	    			{
	    				String ScoreType=this.extract(line,"name");
	    				String value=this.extract(line,"value");
	    				
	    				if(ScoreType.equals("xcorr"))
	    					xcorr=value;
	    				else if(ScoreType.equals("deltacn"))
	    					deltacn=value;
	    				
	    				else if(ScoreType.equals("hyperscore"))
	    					hyperscore=value;
	    				else if(ScoreType.equals("nextscore"))
	    					nextscore=value;
	    				
	    				else if(ScoreType.equals("deltacnstar"))
	    					deltacnstar=value;
	    				else if(ScoreType.equals("spscore"))
	    					spscore=value;
	    				else if(ScoreType.equals("sprank"))
	    					sprank=value;
	    				else if(ScoreType.equals("expect"))
	    					expect=value;
	    			}
	    			else if(line.contains("<alternative_protein"))
	    			{
	    				protein+=";"+this.extract(line,"protein");
	    			}
	    			else if(line.contains("</search_hit>"))
	    			{
	    				if(hflag==0)
	    					Header+="\tModPep\tXcorr\tDeltacn\tDeltacnstar\tSpscore\tLnSprank\tLnEvalue\tEvalue";
	    					    				
	    				info+="\t"+modPEP+"\t"+xcorr+"\t"+deltacn+"\t"+deltacnstar+"\t"+spscore+"\t"+Math.log(Double.parseDouble(sprank))+"\t"+Math.log(Double.parseDouble(expect))+"\t"+expect;
	    				
	    				hitflag=0;
	    				hflag=1;
	    				flag=0;
	    				
	    				this.lists.add(info);
	    				info="";
	    				
	    				if(numPSMperPep.containsKey(modPEP))
		    			{
		    				numPSMperPep.put(modPEP,numPSMperPep.get(modPEP)+1);
		    			}
		    			else
		    			{
		    				numPSMperPep.put(modPEP,1);
		    			}
		    			
	    				PeplistperProt= PUT(protein,modPEP,PeplistperProt);
	    				ProtlistperPep=PUT2(protein, modPEP, ProtlistperPep);
	    				
	    				if(!pepinfo.get(peptide).contains(modPEP))
	    					pepinfo.get(peptide).add(modPEP);	
	    				
	    			}
	    			
	    			else if( line.contains("<mod_aminoacid_mass"))//line.contains("<modification_info") ||
	    			{
	    					String seq=modPEP;
	    					String mass="";
	    					String tmpSeq="";
	    					if(line.contains("static"))
	    						mass=this.extract(line, "static");
	    					else
	    						mass=this.extract(line, "variable");
	    					
	    					int pos=Integer.parseInt(this.extract(line, "position"));
	    						
	    					if(mass.equals("42.0106"))
	    					{
	    						tmpSeq="+42.010600"+seq;
	    					}
	    					else
	    					{	char tmp[]=seq.toCharArray();
	    					
	    						int idpos=0;
	    						//modPEP="";
	    					    for(int i=0;i<seq.length();i++)
	    					    {
	    					    	if(tmp[i]<='Z' && tmp[i]>='A')
	    					    	{
	    					    		idpos++;
	    					    		
	    					    	}
	    					    	
	    					    	if(idpos==pos)
	    					    		tmpSeq+=tmp[i]+"+"+mass;
	    					    	else
	    					    		tmpSeq+=tmp[i];
	    					    }
	     					}
	    				
	    					modPEP=tmpSeq;
	    				//System.out.println(modPEP);
	    				

	    			}
	    			
	    			
	    		}
	    		if(line.contains("</search_result>"))
	    		{
	    			
	    			flag=0;hitflag=0;
	    		}
	    	}
	    }
	    
	    br.close();
			
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

public String extract(String string,String tag) {
		
		String info="";
		String str[]=string.split(" ");
		for(int i=0;i<str.length;i++)
		{
			if(str[i].contains(tag+"="))
			{
				//System.out.println(str[i]);
				info=str[i].split("\\=")[1];

				if(!info.equals("\"\""))
				{
					info=info.substring(1,info.length());
					info=info.split("\"")[0];
				}
				break;
			}
		}

	return info;
	}	
}
