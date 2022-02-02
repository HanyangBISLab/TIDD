/*   TIDD: A Java(TM) based peptide feature extraction software.
 *
 *   Written by H. Li <hllee@hanyang.ac.kr>
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


public class TIDD {
	
	static public long pre_time = System.currentTimeMillis();
	protected static void printUsage(){	
		System.out.println("--------------------------------------------------------------------------");
		System.out.println("Usage    : java -jar TIDD.jar <options> <attributes>");
		System.out.println("Options  :");
		System.out.println("  -i <filename> : PSM search result file [Required]");
		System.out.println("            : tab-deliminated tsv file ");
		System.out.println("  -dir <filename> : MS/MS data file used for search [Required]");
		System.out.println("            : Available formats [mgf] ");
		System.out.println("  -o <output_filename>: Output file path for extracted features [Optional]");
		System.out.println("            : If not specified, it's automatically generated");
		System.out.println("            ");
		System.out.println("  -id <list of index with seperator \",\" >: list of column index in the tab-deleminated PSM result file [required]");
		System.out.println("            : SpectrumFile column index,ScanNum column index,charge column index, precursorMZ column index,peptide column index,Protein column index ");
		System.out.println("            ");
		System.out.println("  -dec <String>: decoy prefix [required]");
		System.out.println("            ");
		System.out.println("  -ft <double>: fragment tolerance (Unit:Da) [required]");
		System.out.println("            ");
		System.out.println("  -phospho <boolean>: true for phospho data [optional]");
		System.out.println("                    : default values is false");
		System.out.println("            ");
		System.out.println();
		System.out.println();
		System.out.println("Example1 : java -jar TIDD.jar -i psm_test.tsv -d C:/Test/ -id 1,2,3,4,5,11 -dec XXX_ -ft 0.025");
		System.out.println();
		
	}
	public TIDD(String[] mode) throws Exception{
		
		String result_file = "";
		String output_file = "";
		String param = "";
		String msms = "";
		String requiredInfo="";
		String addInfo="";
		String decoyPrefix="";
		Boolean phospho=false;
		
		boolean[] range = null;
		double fragTol = 0.02;
	
		for( int i = 0; i < mode.length; i++ ){
			if( mode[i].compareTo("-i") == 0 ){
				i++;
				result_file = mode[i];
				
			}else if( mode[i].compareTo("-id") == 0 ){
				i++;
				requiredInfo = mode[i];
			}else if( mode[i].compareTo("-o") == 0 ){
				i++;
				output_file = mode[i];
			
			}else if( mode[i].compareTo("-dir") == 0 ){
				i++;
				msms = mode[i];
			}else if( mode[i].compareTo("-ft") == 0 ){
				i++;
				fragTol = Double.parseDouble(mode[i]);
			}else if( mode[i].compareTo("-h") == 0 ){
				printUsage();
				System.exit(0);
			}else if( mode[i].compareTo("-dec") == 0 ){
				i++;
				decoyPrefix=mode[i];
			}else if( mode[i].compareTo("-phospho") == 0 ){
				i++;
				if(mode[i].equalsIgnoreCase("true"))
					phospho=true;
			
			}else{
				printUsage();
				System.exit(1);
			}
		}
		if( result_file.isEmpty() || msms.isEmpty() || requiredInfo.split("\\,").length!=6 ){
			printUsage();		
			System.exit(1);
		}
		if( output_file.isEmpty() ){
			output_file = result_file.substring(0, result_file.lastIndexOf('.')) + "_Feature.tsv";
		}
		//-i //home/hllee/Documents/TIDD/modif/phosphp/hela/Hela_phospho_comet/Hela_format.tsv  -dir /home/hllee/Documents/TIDD/modif/phosphp/hela/Hela_newmgf/ -id 1,3,29,4,27,28  -dec DECOY_ -ft 0.025 -phospho true
		//-i /home/hllee/Documents/TIDD/modif/hek/drive-download-20211227T034102Z-001/hek293_modplus.tsv  -dir /home/hllee/Documents/TIDD/modif/hek/drive-download-20211227T034102Z-001/ -id 1,3,5,4,12,15  -dec XXX_ -ft 0.025
		//-i /home/hllee/Documents/TIDD/modif/phosphp/AML_phospho_comet/AML_phospho_format.tsv  -dir /home/hllee/Documents/TIDD/modif/phosphp/mgf/ -id 1,3,29,4,27,28  -dec DECOY_ -ft 0.025 -phospho true
		new ExtractFeature(result_file, msms, output_file, pre_time, fragTol, requiredInfo,decoyPrefix, phospho);
	}

	public static void main(String[] args) throws Exception {
		

		if( args.length == 0 )
			printUsage();
		else
		   new TIDD(args);
	}
}
