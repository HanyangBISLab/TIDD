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

import java.text.DecimalFormat;

public class Constants {
	
  /*
   * Used for xcorr calculation, 0.4 is the default value with low resolution MS2
   * with high resolution, 0.0 should be fine.
   * https://groups.google.com/forum/#!topic/crux-users/xv7kx75zp1s
   */
    public static double binOffset = 0.4;
  
	public static int	 			FragMethod;

	public static String 			SPECTRUM_LOCAL_PATH;
	public static String 			SPECTRUM_FILE_NAME;
	public static String 			INSTRUMENTS_NAME = "QTOF";
	public static int	 			INSTRUMENTS_TYPE = 0; //TOF(0), LOW_TRAP(1), HIGH_TRAP(2)
	public static SpectraFileType 	SPECTRA_FILE_TYPE = SpectraFileType.MGF;	
	public enum SpectraFileType {
		PKL,		// read spectrums in SPECTRUM_FILE_NAME
		DTA,		// read all dta file from SPECTRUM_FILE_NAME(compressed file)
		MGF,		// read spectrums in SPECTRUM_FILE_NAME
		MS2,
		MZXML,		// read spectrums in SPECTRUM_FILE_NAME
		ZIPDTA, 
	}
	
	public enum MSMSType {
		QTOF,
		IONTRAP
	}
	
	public static String		isobaricTag = "";
	public static double[]		reporterMassOfIsobaricTag = null;
	
	public static int			targetDecoy=0;
	public static int			runMODmap = 0;
	
	public static int 			multiStagesSearch = 0;
	public static String 		firstSearchProgram = "";
	
	
	
	public static final double	UNIT_MASS = 1.;
	
	public static final double	Electron = 0.000549;
	public static final double	Hydrogen = 1.007825035;
	public static final double	Oxygen = 15.99491463;
	public static final double	Nitrogen = 14.003074;
	public static final double	Proton = Hydrogen-Electron;
	public static final double	HO = Hydrogen + Oxygen;	
	public static final double	H2O = Hydrogen*2 + Oxygen;	
	public static final double	NH3 = Hydrogen*3 + Nitrogen;		
	public static final double	IsotopeSpace = 1.00235;
	
	public static double		NTERM_FIX_MOD = 0;
	public static double		CTERM_FIX_MOD = 0;
	public static final double	B_ION_OFFSET = Proton;
	public static final double	Y_ION_OFFSET = H2O + Proton;
	public static final double	A_ION_OFFSET = Oxygen + 12.;
	
	public static double		minPeptideMass = 300.;
	public static double		maxPeptideMass = 8000.;//
	
	public static int			numberOfTrypticTermini = 2;
	public static int			missCleavages = 2;
	public static int			trypticSearch = 0;
	
	public static int			NoOfC13 = 0;
	public static int			rangeForIsotopeIncrement = 0;
	
	public static double		alkylatedToCys = 0;
	public static String		alkylationMethod;
	
	public static double		precursorAccuracy = 0.5;
	public static double		precursorTolerance = 5;
	public static double		precursorPPMTolerance = 5;
	public static double		PPMTolerance = 0;
	public static double		fragmentTolerance = 0.6;
	public static double		fragmentPPMTolerance = 0.6;
	public static double		gapTolerance = 0.6;
	public static double		tagChainPruningRate = 0.4;
	public static int			maxTagChainPerPept = 100;	
	public static double		minNormIntensity = 0.00;
	public static double		PTMTolerance = 0.1;
	public static double		maxPTMTolerance = 0.5;
	public static double		minPTMTolerance = 0.5;
	public static double		MaxfragmentTolerance = 0.025;
	
	
	public static boolean multiBlind= true;
	public static boolean ZeroMode=false;
	
	public static double		minModifiedMass = -precursorTolerance;
	public static double		maxModifiedMass = precursorTolerance;
	public static boolean		isInModifiedRange( double v ){
		if( minModifiedMass-gapTolerance < v && v < maxModifiedMass+gapTolerance ) return true;
		else if( Math.abs(v) <= gapTolerance ) return true;
		else return false;
	}
	
	public static int			MSResolution 	= 0; // if 1, high (FT, OrbiTrap)
	public static int			MSMSResolution 	= 0; // if 1, high (FT, OrbiTrap)
	public static int			MSMSAccurateMass = 0;
	
	//for De novo sequencing
	public static double		massToleranceForDenovo = 0.3;
	public static int 			MAX_TAG_SIZE = 100;
	public static double		selectionWindowSize   = 70;
	public static int			minNumOfPeaksInWindow = 4;
	public static int			maxNumOfPeaksInWindow = 8;
	public static int			minTagLength = 3;
	public static int			minTagLengthPeptideShouldContain = 3;
	public static boolean		Leu_indistinguishable_Ile = true;
	public static boolean		Lys_indistinguishable_Qln = true;

	public static String		PTM_FILE_NAME = "PTMDB.xml";
	
	// for Peptide DB
	public static final int		proteinIDModeSeqLength	= 3;
	public static final String 	SOURCE_PROTEIN_FILE_NAME = "sourceProtein.mprot";	
	// for PTM DB
	public static final int 	maxPTMSizePerGap		= 10;
	
	public static final String	UNIMOD_FILE_NAME = "unimod.xml";
	
	// for mother mass correction for LTQ/LCQ
	public static final double	MINIMUM_PRECURSOR_MASS_ERROR = -1.5;
	public static final double	MAXIMIM_PRECURSOR_MASS_ERROR = 1.5;

	// if true, write unidrawing only tag chains whose all gaps are annotated
	public static final boolean	writeAnnotatedTagChainOnly = false;
	
	public static final int		MINIMUM_SHARED_PEAK_COUNT = 2; 
	
	// for offset
	public static final int newLineCharSize = new String("\r\n").getBytes().length;
	
	public static double		nonModifiedDelta = massToleranceForDenovo;
	public static int[] 		maxPTMOccurrence = {1, 1, 2, 2, 3, 3, 3, 2};// = new int[7];
	
	public static int getMaxPTMOccurrence( int seqLength ){
		if( seqLength > 10 || seqLength < 2 ) return 1;
		else return 2;//*/
	}
	
	public static void	adjustParametersForInstrument(int type){
		if( type == 0 ) {
			massToleranceForDenovo = ( MSMSResolution == 0 )? 0.2 : 0.2;	
			minNumOfPeaksInWindow = 3;
			rNorm[0]= 6;
		}
		else {
			massToleranceForDenovo = ( MSMSResolution == 0 )? 0.3 : 0.01;	
			minNumOfPeaksInWindow = 4;
			rNorm[0]= 6;
		}
		if( massToleranceForDenovo > fragmentTolerance/2 )
			massToleranceForDenovo = fragmentTolerance/2;
		
	}
	
	public static boolean	fEqual(double v1, double v2){
		if( Math.abs(v1-v2) <= fragmentTolerance ) return true;
		else return false;
	}
	public static boolean	pEqual(double v1, double v2){
		if( Math.abs(v1-v2) <= precursorTolerance ) return true;
		else return false;
	}	
	
	public static final double[] rNorm= {6,
		2.928968, 1.928968, 1.428968, 1.095635, 0.845635,
		0.645635, 0.478968, 0.336111, 0.211111, 0.100000};
	
	public static double[] coEfft= {0.3159, -34.6288, 1.3209, -8.7609, 0., - 5.0206};		
	public static double getMODScore( double a, double b, double c, double d, double e){
		return coEfft[0]*a + coEfft[1]*b + coEfft[2]*c + coEfft[3]*d + coEfft[4]*e + coEfft[5]; 		
	}
	
	public static String	getString(double value){
		return new DecimalFormat("#.###").format(value).toString();
	}
	
	public static double	MASS_CAL_STD_THRESHOLD = 0.1;
	public static double	PTM_ADD_PENALTY = 0.2;
	

	public static double	getNotExplainedPenaltyWeight(){
		return 0.15;
	}
	public static final	double	ANALYSIS_VERSION = 0.8;
	
	public static boolean isWithinTolerance(double calc, double obsv, double tol){

		if( NoOfC13 == 0 ){
			if( Math.abs(calc-obsv) > tol ) return false;
		}
		else{
			double tempError = obsv - calc;		
			int isoerr = round( tempError / IsotopeSpace );		
			if( isoerr < 0 || NoOfC13 < isoerr ) return false;		
			if(	Math.abs( tempError - isoerr*IsotopeSpace ) > precursorAccuracy ) return false;
		}
		return true;
	}
	
	public static double PPMtoDalton(double mass, double ppm){	
		return mass/1000000*ppm;
	}

	public static int round(double a){
		if( a > 0 ) return (int)(a + 0.5);
		else return (int)(a - 0.5);
	}
}
