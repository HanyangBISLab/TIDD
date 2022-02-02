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
import java.util.List;

public class XCorr {
 
    public static int[] CreateEvidenceVector(ArrayList<Peak> peak, double peptideMass, double obMW, double binWidth, int charge) {
        
        int MAX_XCORR_OFFSET = 75;
        int maxPrecurMass = mass2bin(peptideMass, 1, binWidth);
        double pepMassMonoMean = (maxPrecurMass - 0.5 + Constants.binOffset) * binWidth;
        maxPrecurMass = mass2bin(obMW + MAX_XCORR_OFFSET + 30, 1, binWidth) + 50;

        double evidenceIntScale = 500.0; //original value is 500
        int nRegion = 10; 
        double maxIntensPerRegion = 50.0;
        double precursorMZExclude = 15.0;
        double massNH3Mono = 17.02655;     // mass of NH3 (monoisotopic)
        double massCOMono =  27.9949;      // mass of CO (monoisotopic)
        double massH2OMono = 18.010564684; // mass of water (monoisotopic)
        double massHMono = 1.0078246;      // mass of hydrogen (monoisotopic)
        double BYHeight = 50.0;
        double NH3LossHeight = 10.0;    
        double COLossHeight = 10.0;    // for creating a ions on the fly from b ions
        double H2OLossHeight = 10.0;

        int ma;
        int pc;
        int ionBin;
        double bIonMass;
        double yIonMass;
        double ionMZMultiCharge;
        double ionMassNH3Loss;
        double ionMassCOLoss;
        double ionMassH2OLoss;

        double[] evidence = new double[maxPrecurMass];
        double[] intensArrayObs = new double[maxPrecurMass];
        int[] intensRegion = new int[maxPrecurMass];
        int[] evidenceInt = new int[maxPrecurMass];

        for (ma = 0; ma < maxPrecurMass; ma++) {
            evidence[ma] = 0.0;
            evidenceInt[ma] = 0;
            intensArrayObs[ma] = 0.0;
            intensRegion[ma] = -1;
        }

        double precurMz = obMW/charge + Constants.Proton;;
        int nIon = peak.size();
        int precurCharge = charge;
        double experimentalMassCutoff = precurMz * precurCharge + 50.0;
        
        double maxIonMass = 0.0;
        double maxIonIntens = 0.0;
        for (int ion = 0; ion < nIon; ion++) {
            double ionMass = peak.get(ion).getMz();
            double ionIntens = peak.get(ion).getIntensity();
            if (ionMass >= experimentalMassCutoff) {
                continue;
            }
            if (maxIonMass < ionMass) {
                maxIonMass = ionMass;
            }
            if (maxIonIntens < ionIntens) {
                maxIonIntens = ionIntens;
            }
        }
        int regionSelector = (int)Math.floor( (int)Math.floor(((double)maxIonMass / binWidth) + 1.0 - Constants.binOffset ) / (double)nRegion);
        for (int ion = 0; ion < nIon; ion++) {
            double ionMass = peak.get(ion).getMz();
            double ionIntens = peak.get(ion).getIntensity();
            if (ionMass >= experimentalMassCutoff) {
                continue;
            }
            if (ionMass > precurMz - precursorMZExclude && ionMass < precurMz + precursorMZExclude) {
                continue;
            }
            ionBin = (int)Math.floor((ionMass / binWidth) + 1.0 - Constants.binOffset); //REVIEW same as mass2bin method
            int region = (int)Math.floor((double)(ionBin) / (double)regionSelector);
            if (region >= nRegion) {
                region = nRegion - 1;
            }
            intensRegion[ionBin] = region;
            if (intensArrayObs[ionBin] < ionIntens) {
                intensArrayObs[ionBin] = ionIntens;
            }
        }

        maxIonIntens = Math.sqrt(maxIonIntens);
        for (ma = 0; ma < maxPrecurMass; ma++) {
            intensArrayObs[ma] = Math.sqrt(intensArrayObs[ma]);
            if (intensArrayObs[ma] <= 0.05 * maxIonIntens) {
                intensArrayObs[ma] = 0.0;
            }
        }

        double[] maxRegion = new double[nRegion];
        for (int re = 0; re < nRegion; re++) {
            maxRegion[re] = 0.0;
        }
        for (ma = 0; ma < maxPrecurMass; ma++) {
            int reg = intensRegion[ma];
            if (reg >= 0 && maxRegion[reg] < intensArrayObs[ma]) {
                maxRegion[reg] = intensArrayObs[ma];
            }
        }
        for (ma = 0; ma < maxPrecurMass; ma++) {
            int reg = intensRegion[ma];
            if (reg >= 0 && maxRegion[reg] > 0.0) {
                intensArrayObs[ma] *= (maxIntensPerRegion / maxRegion[reg]);
            }
        }
        //delete [] maxRegion;

        // ***** Adapted from tide/spectrum_preprocess2.cc.
        // TODO replace, if possible, with call to 
        // static void SubtractBackground(double* observed, int end).
        // Note numerous small changes from Tide code.
        double multiplier = 1.0 / (MAX_XCORR_OFFSET * 2.0 + 1.0);
        double total = 0.0;
        double[] partial_sums = new double[maxPrecurMass];
        for (int i = 0; i < maxPrecurMass; ++i) {
            partial_sums[ i ] = ( total += intensArrayObs[i]); //cumulative sum
        }
        for (int i = 0; i < maxPrecurMass; ++i) {
            int right_index = Math.min(maxPrecurMass - 1, i + MAX_XCORR_OFFSET);
            int left_index = Math.max(0, i - MAX_XCORR_OFFSET - 1);
            intensArrayObs[i] -= multiplier * (partial_sums[right_index] - partial_sums[left_index]); //REVIEW just diminish the value a little
        }
        //delete [] partial_sums;
        // *****
        int binFirst = mass2bin(30, 1, binWidth);
        int binLast = mass2bin(pepMassMonoMean - 47, 1, binWidth);

        for (ma = binFirst; ma <= binLast; ma++) {
            // b ion
            bIonMass = (ma - 0.5 + Constants.binOffset) * binWidth;
            ionBin = (int)Math.floor(bIonMass / binWidth + 1.0 - Constants.binOffset);
            evidence[ma] = evidence[ma] + intensArrayObs[ionBin] * BYHeight;
            for (pc = 3; pc <= precurCharge; pc++) {
                ionBin = mass2bin(bIonMass, pc-1, binWidth);
                evidence[ma] = evidence[ma] + intensArrayObs[ionBin] * BYHeight;
            }
            // y ion
            yIonMass = pepMassMonoMean + 2 * massHMono - bIonMass;
            ionBin = (int)Math.floor(yIonMass / binWidth + 1.0 - Constants.binOffset);
            evidence[ma] = evidence[ma] + intensArrayObs[ionBin] * BYHeight;
            for (pc = 3; pc <= precurCharge; pc++) {
                ionBin = mass2bin(yIonMass, pc-1, binWidth);
                evidence[ma] = evidence[ma] + intensArrayObs[ionBin] * BYHeight;
            }
            // NH3 loss from b ion
            ionMassNH3Loss = bIonMass - massNH3Mono;
            ionBin = (int)Math.floor(ionMassNH3Loss / binWidth + 1.0 - Constants.binOffset);
            evidence[ma] = evidence[ma] + intensArrayObs[ionBin] * NH3LossHeight;
            for (pc = 3; pc <= precurCharge; pc++) {
                ionBin = mass2bin(ionMassNH3Loss, pc-1, binWidth);
                evidence[ma] = evidence[ma] + intensArrayObs[ionBin] * NH3LossHeight;
            }
            // NH3 loss from y ion
            ionMassNH3Loss = yIonMass - massNH3Mono;
            ionBin = (int)Math.floor(ionMassNH3Loss / binWidth + 1.0 - Constants.binOffset);
            evidence[ma] = evidence[ma] + intensArrayObs[ionBin] * NH3LossHeight;
            for (pc = 3; pc <= precurCharge; pc++) {
                ionBin = mass2bin(ionMassNH3Loss, pc-1, binWidth);
                evidence[ma] = evidence[ma] + intensArrayObs[ionBin] * NH3LossHeight;
            }
            // CO and H2O loss from b ion
            ionMassCOLoss = bIonMass - massCOMono;
            ionMassH2OLoss = bIonMass - massH2OMono;
            ionBin = (int)Math.floor(ionMassCOLoss / binWidth + 1.0 - Constants.binOffset);
            evidence[ma] = evidence[ma] + intensArrayObs[ionBin] * COLossHeight;
            ionBin = (int)Math.floor(ionMassH2OLoss / binWidth + 1.0 - Constants.binOffset);
            evidence[ma] = evidence[ma] + intensArrayObs[ionBin] * H2OLossHeight;
            for (pc = 3; pc <= precurCharge; pc++) {
                ionBin = mass2bin(ionMassCOLoss, pc-1, binWidth);
                evidence[ma] = evidence[ma] + intensArrayObs[ionBin] * COLossHeight;
                ionBin = mass2bin(ionMassH2OLoss, pc-1, binWidth);
                evidence[ma] = evidence[ma] + intensArrayObs[ionBin] * H2OLossHeight;
            }
            // H2O loss from y ion
            ionMassH2OLoss = yIonMass - massH2OMono;
            ionBin = (int)Math.floor(ionMassH2OLoss / binWidth + 1.0 - Constants.binOffset);
            evidence[ma] = evidence[ma] + intensArrayObs[ionBin] * H2OLossHeight;
            for (pc = 3; pc <= precurCharge; pc++) {
                ionBin = mass2bin(ionMassH2OLoss, pc-1, binWidth);
                evidence[ma] = evidence[ma] + intensArrayObs[ionBin] * H2OLossHeight;
            }
        }

        // discretize evidence array
        for (ma = 0; ma < maxPrecurMass; ma++) {
            evidenceInt[ma] = (int)Math.floor(evidence[ma] / evidenceIntScale + 0.5);
        }

        return evidenceInt;
    }
    
    public static int[] generateThreoPeak(String seq, double binWidth) {
        
        int arraysize = seq.length() - 1;
        int threopeak[] = new int[arraysize];
        
        double mass = 0;
        for(int i=0; i<arraysize; i++) {
            mass += AminoAcid.getMass(seq.charAt(i));
            int ma = mass2bin(mass + Constants.Proton, 1, binWidth);
            threopeak[i] = ma;
        }
        
        return threopeak;
    }

    static int mass2bin(double mass, double charge, double binWidth){
        return (int)Math.floor( ((mass + (charge-1)*Constants.Proton)/(charge*binWidth)) + 1.0 - Constants.binOffset ); //REVIEW what is binOffset? the value is always zero.
    }

    public static int[] generateThreoPeak(String seq, double binWidth, List<Modification> modificationList) {
        
 //   	System.out.println(seq);
        int fragmentationSiteSize = seq.length() - 1;
        int threopeak[] = new int[fragmentationSiteSize];
        
        double mass = 0;
        int index = 0;
        
        if (modificationList.size() != 0) {
          //if it has n-term modification, the n-term mass should be added.
          if (modificationList.get(0).location == 0) { //meaning n-term modification.
            mass += modificationList.get(index).Mass;
            index++;
          }
        }
        
        for (int aminoIndex = 0; aminoIndex < fragmentationSiteSize; aminoIndex++) {
            char aa = seq.charAt(aminoIndex);
            mass += AminoAcid.getMass(aa);
            if(index < modificationList.size() && modificationList.get(index).location - 1 == aminoIndex) { // if it's not a N-term modification, location - 1 >= 0
                mass += modificationList.get(index).Mass;
                index++;
            }
            int ma = mass2bin(mass + Constants.Proton, 1, binWidth);
            threopeak[aminoIndex] = ma;
        }
        
        return threopeak;
    }

    public static double getPeptideMass(String strippedPeptideSequence,
        ArrayList<Modification> modList) {
      
      double pMass = 0;
      //Add all modification mass
      for (Modification modification : modList) {
        pMass += modification.Mass;
      }
      //Add all residue mass
      for (char residue : strippedPeptideSequence.toCharArray()) {
        pMass += AminoAcid.getMass(residue);
      }
      pMass += Constants.H2O;
      
      return pMass;
    }
}