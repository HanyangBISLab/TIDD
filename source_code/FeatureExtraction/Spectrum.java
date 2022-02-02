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
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.StringTokenizer;
import java.util.Vector;

import dst.SpectrumT;

//	spectrum�� �� �����͸� �ٷ�� �Լ����� ������ �ִ�.
public class Spectrum {
	public final int MINRANK = 11;	//	peak�� annotation rank ����
	public double TOLERANCE = 0.6;	//	 ������
	private SpectrumT spectrum;		//	spectrum ��ü
	private int[] pIndex;
	private double[] inten;
	
	public int specIndex = -1;
	public int scanNum = -1;
	private static BufferedReader file_pointer;
	
	//	���� �о���� spectrum�� ��ȯ�Ѵ�.
	public SpectrumT getProcessedSpectrum(File mgf, int offset, double MW, double pmz, double frag_tol){
		String line = "";
		
		double baseInt = .0, baseMZ = .0, TIC = .0, SME = .0;
		double tarMass = .0, tarInten = .0, tolerance = frag_tol/2;
		//System.out.println("OK");
		
		//	mgf file ���� ��ġ�� ã�ư� mass�� intensity���� �о���δ�.
		try {
			spectrum = new SpectrumT();
			//RandomAccessFile rf = new RandomAccessFile(mgf, "r");
			BufferedReader br = null;
			try{
				br = new BufferedReader(new FileReader(mgf));
			}catch(Exception e){
				System.out.println("Cannot find MSMS file");
				System.out.println(mgf.getName());
				return null;
			}
			
			//long os = offset;
			
//			spectrum.maxME = -99999999.0;
//			spectrum.minME = 99999999.0;
			
			spectrum.peakMassList = new Vector<Double>();
			spectrum.peakIntensityList = new Vector<Double>();
			
			spectrum.originalPeakMassList = new Vector<Double>();
			spectrum.originalPeakIntensityList = new Vector<Double>();
			
			//rf.seek(os);
			br.skip(offset);
			while(!(line = br.readLine()).contains("END")){
				double mass = .0, intensity = .0;
				line.trim();
				if( line.isEmpty() )
					continue;
				
				//	mass�� intensity�� �����ϰ� �ִ� �����Ϳ��� ���� �޾ƿ� (�����Ϳ� ���� �ڿ� �߰��� ���� �ٴ� ��쵵 ������ �ٷ��� ����)
				StringTokenizer tokenizer = new StringTokenizer(line);
				mass = Double.parseDouble(tokenizer.nextToken());
//				System.out.println(line);
				intensity = Double.parseDouble(tokenizer.nextToken());
				
				/*
				 * 
				 * 2014-01-15���� �ڵ�
				 * //	m/z ���� intensity�� ������ split ���ڸ� �̿��Ͽ� �� �и�
				String[] temp = line.split("\\t");
				if( temp.length != 2 )	temp = line.split(" ");
				if( temp.length != 2 )	break;
				mass = Double.parseDouble(temp[0]);
				intensity = Double.parseDouble(temp[1]);
				*/
				
				/*
				 * 	2014-05-29�ڵ� �߰�
				 * 
				 *  ���� �����ִ� spectrum�� ����ϱ� ���� ���� ������
				 */
				spectrum.originalPeakMassList.add(mass);
				spectrum.originalPeakIntensityList.add(intensity);
				
				if( mass > MW )	break;
				
				if( mass < 1 || Math.abs(mass-pmz) < 3 ) continue;
				
				
				if( intensity != 0 ){
					if( (mass-tarMass) < tolerance ){
						double sum = tarInten + intensity;
						tarMass = tarMass * (tarInten / sum) + mass * (intensity / sum);
						tarInten += intensity;
						spectrum.peakMassList.add(tarMass);
						spectrum.peakIntensityList.add(tarInten);
						SME += mass-tarMass;
//						if( spectrum.minME > mass-tarMass ){
//							spectrum.minME = mass-tarMass;
//						}if( spectrum.maxME < mass-tarMass ){
//							spectrum.maxME = mass-tarMass;
//						}
					}else{
						tarMass = mass;
						tarInten = intensity;
						spectrum.peakMassList.add(tarMass);
						spectrum.peakIntensityList.add(tarInten);
					}
					
					TIC += intensity;
					if( tarInten > baseInt ){
						baseInt = tarInten;
						baseMZ = tarMass;
					}
				}
				
				//System.out.println(mass + " " +  intensity);
			}
			spectrum.peakMassList.trimToSize();
			spectrum.peakIntensityList.trimToSize();
			
			
			
//			spectrum.peakMassList = spectrum.originalPeakMassList;
//			spectrum.peakIntensityList = spectrum.originalPeakIntensityList;
			
			spectrum.TIC = TIC;
			spectrum.peaksCount = spectrum.peakMassList.size();
			spectrum.baseIntensity = baseInt;
			spectrum.baseMZ = baseMZ;
			spectrum.mzLimit = ((int)(MW/100)+1)*100;
//			spectrum.SME = SME;
			//System.out.println(spectrum.peaksCount);
			
			setPeakRankInLocalArea();
			
			//for(int i = 0; i< spectrum.localRank.length; i ++){
			//	System.out.println(spectrum.localRank[i]);
			//}
			
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return spectrum;
	}
	
	public SpectrumT getProcessedSpectrum(BufferedReader br, double MW, double pmz, double frag_tol){
		String line = "";
		
		double baseInt = .0, baseMZ = .0, TIC = .0, SME = .0;
		double tarMass = .0, tarInten = .0, tolerance = frag_tol/2;
		//System.out.println("OK");
		
		//	mgf file ���� ��ġ�� ã�ư� mass�� intensity���� �о���δ�.
		try {
			spectrum = new SpectrumT();
			//RandomAccessFile rf = new RandomAccessFile(mgf, "r");
						
			//long os = offset;
			
//			spectrum.maxME = -99999999.0;
//			spectrum.minME = 99999999.0;
			
			spectrum.peakMassList = new Vector<Double>();
			spectrum.peakIntensityList = new Vector<Double>();
			
			spectrum.originalPeakMassList = new Vector<Double>();
			spectrum.originalPeakIntensityList = new Vector<Double>();
			
			while(!(line = br.readLine()).contains("END")){
				double mass = .0, intensity = .0;
				line.trim();
				if( line.isEmpty() )
					continue;
				
				//	mass�� intensity�� �����ϰ� �ִ� �����Ϳ��� ���� �޾ƿ� (�����Ϳ� ���� �ڿ� �߰��� ���� �ٴ� ��쵵 ������ �ٷ��� ����)
				StringTokenizer tokenizer = new StringTokenizer(line);
				mass = Double.parseDouble(tokenizer.nextToken());
//				System.out.println(line);
				intensity = Double.parseDouble(tokenizer.nextToken());
				
				/*
				 * 
				 * 2014-01-15���� �ڵ�
				 * //	m/z ���� intensity�� ������ split ���ڸ� �̿��Ͽ� �� �и�
				String[] temp = line.split("\\t");
				if( temp.length != 2 )	temp = line.split(" ");
				if( temp.length != 2 )	break;
				mass = Double.parseDouble(temp[0]);
				intensity = Double.parseDouble(temp[1]);
				*/
				
				/*
				 * 	2014-05-29�ڵ� �߰�
				 * 
				 *  ���� �����ִ� spectrum�� ����ϱ� ���� ���� ������
				 */
				spectrum.originalPeakMassList.add(mass);
				spectrum.originalPeakIntensityList.add(intensity);
				
				if( mass > MW )	break;
				
				if( mass < 1 || Math.abs(mass-pmz) < 3 ) continue;
				
				
				if( intensity != 0 ){
					if( (mass-tarMass) < tolerance ){
						double sum = tarInten + intensity;
						tarMass = tarMass * (tarInten / sum) + mass * (intensity / sum);
						tarInten += intensity;
						spectrum.peakMassList.add(tarMass);
						spectrum.peakIntensityList.add(tarInten);
						SME += mass-tarMass;
//						if( spectrum.minME > mass-tarMass ){
//							spectrum.minME = mass-tarMass;
//						}if( spectrum.maxME < mass-tarMass ){
//							spectrum.maxME = mass-tarMass;
//						}
					}else{
						tarMass = mass;
						tarInten = intensity;
						spectrum.peakMassList.add(tarMass);
						spectrum.peakIntensityList.add(tarInten);
					}
					
					TIC += intensity;
					if( tarInten > baseInt ){
						baseInt = tarInten;
						baseMZ = tarMass;
					}
				}
				
				//System.out.println(mass + " " +  intensity);
			}
			spectrum.peakMassList.trimToSize();
			spectrum.peakIntensityList.trimToSize();
			
			spectrum.TIC = TIC;
			spectrum.peaksCount = spectrum.peakMassList.size();
			spectrum.baseIntensity = baseInt;
			spectrum.baseMZ = baseMZ;
			spectrum.mzLimit = ((int)(MW/100)+1)*100;
//			spectrum.SME = SME;
			//System.out.println(spectrum.peaksCount);
			
			setPeakRankInLocalArea();
			
			//for(int i = 0; i< spectrum.localRank.length; i ++){
			//	System.out.println(spectrum.localRank[i]);
			//}
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		return spectrum;
	}
	
	//	target mass���� ���� peak�� index�� ã�� �Լ�
	public int getMatchedPeak(SpectrumT spectrum, double tarPeak, double minimum){
		double maxPeak = -1;//, maxDelta = TOLERANCE+1;
		int matIndex = -1;
		
		//	target m/z���� ���������� �� ���� �����ʿ� �����ϴ� ���� ���� ���� ã�´�.
		int index = biSearchPeak(spectrum.peakMassList, tarPeak - TOLERANCE, spectrum.peaksCount);
		//System.out.println(spectrum.peakMassList.length + "  " + index + "  " + tarPeak);
		if( spectrum.peakMassList.get(index) < tarPeak - TOLERANCE ){
			return matIndex;
			//return index;
		}
		//int index = 0;
		
		
		//	target m/z���� ��� ������ ���� ������ ���� ���� ������ ���� abundant�� peak�� ã�� �� index�� ��ȯ�Ѵ�.
		while( spectrum.peakMassList.get(index) <= tarPeak + TOLERANCE ){
			if( spectrum.peakIntensityList.get(index) > maxPeak ){
				//	���� ������ ���� ���� �̻��� intensity�� ���� peak�鿡 ���ؼ��� ã�´�.
				if( spectrum.localRank[index] < MINRANK && spectrum.peakIntensityList.get(index) > spectrum.baseIntensity * minimum ){
					maxPeak = spectrum.peakIntensityList.get(index);
					matIndex = index;
				}
			}
			
			index ++;
			if( index == spectrum.peaksCount )
				break;
		}
		return matIndex;
	}
	
	//	binary search������� target�� ã�´�.
	public int biSearchPeak(Vector<Double> pairs, double left, int count) {
		int index;
		if( left <= pairs.get(0) )
			index = 0;
		else if( left > pairs.get(count-1) ){
			index = count-1;
		}
		else{
			int M, L =0, R = count-1;
			while( R - L > 1 ){
				M = (L + R) / 2;
				
				if( left <= pairs.get(M) )
					R = M;
				else
					L = M;
			}
			index = R;
		}
		return index;
	}
	
	//	������ ���� ���� peak�鿡 ranking�� �Ű� ��ȯ�Ѵ�.
	private void setPeakRankInLocalArea() {
		int i = 0;
		int localSize = 100; 
		
		//int j = 0;
		spectrum.localRank = new int[spectrum.peaksCount];
		
		int localIndex = (int)(spectrum.peakMassList.get(i)/localSize) + 1;
		while( i < spectrum.peaksCount  ){
			int go = 0;
			
			pIndex = new int[spectrum.peaksCount];
			inten = new double[spectrum.peaksCount];
			while( localIndex > spectrum.peakMassList.get(i)/localSize){
				pIndex[go] = i;
				inten[go] = spectrum.peakIntensityList.get(i);
				i ++;
				go++;
				if(i >= spectrum.peaksCount)
					break;
			}
			sortSpectrumByIntencity( 0, go-1);
			for(int k = 0; k < go; k++){
				spectrum.localRank[pIndex[k]] = go - (k+1);
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
	
	
	
	//	�� index��, �� intensity�� swap�Ѵ�. 
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
	
	 //	quick sort�� intensity ������ index�� �����Ѵ�.
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

	public void init(){
		specIndex = -1;
		scanNum = -1;
	}
	
	public BufferedReader getOffset(int index, int scan, String msms) {
		// TODO Auto-generated method stub
		try {
			
//			if( specIndex == -1 ){
//				specIndex = 0;
//				file_pointer = new BufferedReader(new FileReader(msms));
//			}
			
//			if( specIndex == index ){
//				return file_pointer;
//			}
			
			if( scanNum == -1 ){
				file_pointer = new BufferedReader(new FileReader(msms));
			}
			
			if( scanNum == scan ){
				return file_pointer;
			}
			String line;
			while( (line = file_pointer.readLine()) != null ){
				if( line.startsWith("BEGIN") ){
					specIndex++;
//					System.out.println(specIndex);
				}
				
				if( line.startsWith("TITLE") ){
					scanNum = Integer.parseInt(line.split("\\.")[1]);
				}
				
				if( scanNum == scan ){
					while( !line.startsWith("CHARGE") ) line = file_pointer.readLine();
					break;
				}
				
//				if( specIndex == index ){
//					while( !line.startsWith("CHARGE") ) line = file_pointer.readLine();
//					break;
//				}
				
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return file_pointer;
	}
	
}
