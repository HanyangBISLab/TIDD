/*   Xcorr: A Java(TM) based PSM similarity score-Xcorr calculation software.
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

package xcorr;

import java.util.ArrayList;
import java.io.*;

import java.util.*;
/*
 * You need to assign those parameters below 
 * 1 strip peptide sequence 
 * 2 peak (m/z and intensity) 
 * 3 modification (modification location and 'exact' modification mass) 
 * 4 obMW 
 * 5 charge
 */


/**
 * output : Xcorr column added txt file
 * 
 * @author Jo
 *
 */
public class Run {

  private static int TITLE_COL = -1;
  private static int PEPTIDE_COL = -1;
  private static int PRECURSOR_COL = -1;
  private static int CHARGE_COL = -1;
  private static int MGFNAME_COL=-1;
  private static String HEADER="";

  private int resultCnt = 0;
  private int processedScanCnt = 0;
  
  private static boolean containsHeaderInResultFile = true;
  /*
   * In tide search, for low resolution, the default is 1.0005079. for high-resolution data, 0.02 is
   * recommended.
   */
  private static double binWidth = 1.0005079;
  private final int weightScaling = 20;

  public static void main(String[] args) throws IOException, InterruptedException {

    if (args.length != 10) {
      System.out.println("Invalid arguments error, please check the arguments\n");
      printUsage();
      System.exit(-1);
    }

    String spectrumFile = args[0];
    String resultFile = args[1];
    binWidth = Double.parseDouble(args[2]);
    Constants.binOffset = Double.parseDouble(args[3]);

    TITLE_COL = Integer.parseInt(args[4]);
    PEPTIDE_COL = Integer.parseInt(args[5]);
    PRECURSOR_COL = Integer.parseInt(args[6]);
    CHARGE_COL = Integer.parseInt(args[7]);
    MGFNAME_COL=Integer.parseInt(args[8]);
    
    if (Boolean.parseBoolean(args[9]) == false) {
    	containsHeaderInResultFile = false;
    }
    
    new Run(resultFile, spectrumFile);
  }

  private static void printUsage() {
    System.out.println(
        "Usage: java -jar xcorr.jar spectrumDir resultFile binWidth binOffset titleCol peptideCol precursorMZCol chargeCol mgffilenameCol header");
    System.out.println("-----Arguments description------");
    System.out.println("spectrumFile should be a mgf directory");
    System.out.println("resultFile should be a tsv file");
    System.out.println("Column number is a zero-based number");
    System.out.println(
        "binWidth is default is 1.0005079 for low resolution. 0.02 is recommended for high resolution");
    System.out
        .println("binOffset defulat is 0.4 for low resolution. 0.0 is fine for high resolution");
    System.out.println("header value can be either true (result file has a header) or false (no header in the result file)");
  }

  public Run(String resultFile, String spectrumFile) throws IOException {

    System.out.println(resultFile + " is being processed");

    BufferedReader brResultFile = new BufferedReader(new FileReader(resultFile));
    
    BufferedWriter bwResultOutput = new BufferedWriter(new FileWriter(
        resultFile.substring(0, resultFile.lastIndexOf(".")) + "_xcorrColAdded.tsv"));

    // Pass header
    if (containsHeaderInResultFile) {
    	HEADER=brResultFile.readLine();
    	bwResultOutput.write(HEADER+"\t"+"CalXcorr\n");
    }

    String resultLine, specLine;
    String[] splitResultLine;

    // Create a map
    // key: title
    // value: a result line (tsv file)
    int resultCount = 0;
    
    ArrayList<String> mgffiles=new ArrayList<String>();
    
    Map<String, String> specTitleMap = new HashMap<String, String>();
    while ((resultLine = brResultFile.readLine()) != null) {
	  resultCount ++;
      specTitleMap.put(resultLine.split("\t")[TITLE_COL], resultLine);
      if(!mgffiles.contains(resultLine.split("\t")[MGFNAME_COL]))
    	  mgffiles.add(resultLine.split("\t")[MGFNAME_COL]);
    }
    System.out.println("Total PSM count:" + resultCount);
    System.out.println("Total number of MGF files:"+mgffiles.size());
   
    
    for(int m=0;m<mgffiles.size();m++)
    {
    	BufferedReader brSpecFile = new BufferedReader(new FileReader(spectrumFile+mgffiles.get(m)));
    	  while ((specLine = brSpecFile.readLine()) != null) {

    	      if (specLine.startsWith("TITLE")) {
    	        String specTitle = specLine.split("=")[1];

    	        if (specTitleMap.containsKey(specTitle)) {

    	          ArrayList<Peak> peakList = new ArrayList<>();

    	          while ((specLine = brSpecFile.readLine()).contains("=")) {
    	            // Read line until peak list
    	          } ;

    	          Peak rawPeak = new Peak(Double.parseDouble(specLine.split("\\s+")[0]),
    	              Double.parseDouble(specLine.split("\\s+")[1]));
    	          peakList.add(rawPeak);
    	          while (!(specLine = brSpecFile.readLine()).startsWith("END")) {
    	            rawPeak = new Peak(Double.parseDouble(specLine.split("\\s+")[0]),
    	                Double.parseDouble(specLine.split("\\s+")[1]));
    	            peakList.add(rawPeak);
    	          }

    	          if (peakList.size() == 0) {
    	            System.out.println("Couldn't find peak list. Title : " + specTitle);
    	            System.exit(0);
    	          }

    	          if (specTitleMap.containsKey(specTitle)) {
    	            processedScanCnt++;
    	  		    if ((processedScanCnt % 100) == 0) {
    	  	          System.out.println(processedScanCnt + "/" + resultCount);
    	  		    }

    	            resultLine = specTitleMap.get(specTitle);

    	            splitResultLine = resultLine.split("\t");

    	            String title = splitResultLine[TITLE_COL];
    	            String peptide = splitResultLine[PEPTIDE_COL];
    	            String precursor = splitResultLine[PRECURSOR_COL];
    	            int charge = Integer.parseInt(splitResultLine[CHARGE_COL]);

    	            String strippedPeptideSequence = getStripSeq_IL_Replaced(peptide);
    	            ArrayList<Modification> modList = getModiList(peptide);

    	            double obMW = (Double.parseDouble(precursor) - Constants.Proton) * charge;

    	            int threoPeak[] = XCorr.generateThreoPeak(strippedPeptideSequence, binWidth, modList);
    	            double peptideMass = XCorr.getPeptideMass(strippedPeptideSequence, modList);
    	            int evidence[] = XCorr.CreateEvidenceVector(peakList, peptideMass, obMW, binWidth, charge);

    	            int refactoredXcorrValue = 0;
    	            for (int i = 0; i < threoPeak.length; i++) {
    	              if (threoPeak[i] >= evidence.length) {
    	                break;
    	              }
    	              refactoredXcorrValue += evidence[threoPeak[i]];
    	            }
    	            bwResultOutput
    	                .write(resultLine + "\t" + (double) refactoredXcorrValue / weightScaling + "\n");
    	          }
    	        }
    	      }

    	    }
    	  brSpecFile.close();
    }
  
    System.out.println(processedScanCnt + "/" + resultCount);
    
    System.out.println("Processed scan count : " + processedScanCnt);
    System.out.println(resultFile.substring(0, resultFile.lastIndexOf(".")) + "_xcorrColAdded.tsv"
        + " is generated!");
    System.out.println();


    
    brResultFile.close();
    bwResultOutput.close();
  }


  public void normalizePeakList(ArrayList<Peak> peakList) {
    double basePeakInten = 0;

    for (Peak p : peakList) {
      if (p.getIntensity() > basePeakInten) {
        basePeakInten = p.getIntensity();
      }
    }

    for (Peak p : peakList) {
      p.setNormIntensity(p.getIntensity() / basePeakInten);
    }

  }

  public ArrayList<Modification> getModiList(String peptide) {
    ArrayList<Modification> modiList = new ArrayList<>();

    int index = 0;
    int AAIndex = 1;
    String modiMass = "";
    String modiMass2 = "";
    
    // Nterm modification
    int nTermModCount = 0;
    boolean hasTwoNtermModi = false;
    
    if (peptide.startsWith("+")) {
      if (!Character.isAlphabetic(peptide.charAt(index))) {
        while (!Character.isAlphabetic(peptide.charAt(index))) {
          
          if (peptide.charAt(index) == '+' || peptide.charAt(index) == '-') {
            nTermModCount += 1;

            if (nTermModCount > 1) {
              hasTwoNtermModi = true;
            }
          }
          
          if (hasTwoNtermModi == false) {
            modiMass += peptide.charAt((index++));
          } else {
            modiMass2 += peptide.charAt((index++));
          }
        }
        
        Double modMass = 0.0;
        if (hasTwoNtermModi == false) {
        	modMass = Double.parseDouble(modiMass);
        }
        else {
        	modMass = Double.parseDouble(modiMass) + Double.parseDouble(modiMass2);
        }
        Modification modi = new Modification("", "", 0, 0, "", modMass);
        modiList.add(modi);
      }
    }

    // The other residues
    modiMass = "";
    modiMass2 = "";
    boolean hasTwoModiInTheResidue = false;
    int residueModCount = 0;
    
    for (; index < peptide.length(); index++) {
      modiMass = "";
      residueModCount = 0;
      hasTwoModiInTheResidue = false;
      
      if (!Character.isAlphabetic(peptide.charAt(index))) {
        while (index < peptide.length() && !Character.isAlphabetic(peptide.charAt(index))) {
          
          if (peptide.charAt(index) == '+' || peptide.charAt(index) == '-') {
            residueModCount += 1;

            if (residueModCount > 1) {
              hasTwoModiInTheResidue = true;
            }
          }
          
          if (hasTwoModiInTheResidue == false) {
            modiMass += peptide.charAt((index++));
          } else {
            modiMass2 += peptide.charAt((index++));
          }
        }
        
        Double modMass = 0.0;
        if (hasTwoModiInTheResidue == false) {
          modMass = Double.parseDouble(modiMass);
        }
        else {
          modMass = Double.parseDouble(modiMass) + Double.parseDouble(modiMass2);
        }

        Modification modi = new Modification("", "", 0, AAIndex - 1, "", modMass);
        modiList.add(modi);
      }
      AAIndex++;
    }
    
    return modiList;
  }

  public static String getStripSeq_IL_Replaced(String peptide) {
    String stripSeq = "";

    for (int i = 0; i < peptide.trim().length(); i++) {
      if (Character.isLetter(peptide.charAt(i))) {
        stripSeq += peptide.charAt(i);
      }
    }
    return stripSeq;
  }

}
