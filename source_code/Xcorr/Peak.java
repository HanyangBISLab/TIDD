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

import java.util.Comparator;

public class Peak implements Comparable<Peak> {
	double mz;
	double intensity;
	double normIntensity=0;
	int charge = 1;
	int index = -1;
	
	public String info = "";
	
	public Peak(int index, double mass, double intensity, int charge)
	{
		this.index = index;
		this.mz = mass;
		this.intensity = intensity;
		this.charge = charge;
	}
	
	public Peak(int index, double mass, double intensity) 
	{
		this(index, mass, intensity, 1);
	}

	public Peak(double mz, double intensity) {
		this.mz = mz;
		this.intensity = intensity;
	}

	public	int				getIndex()			{ return index; }
	public	double			getMz() 			{ return mz; }
	public	double			getIntensity() 		{ return intensity; }
	public	double			getNormIntensity()	{ return normIntensity; }
	public	int				getCharge() 		{ return charge; }
	public	double	getComplementMass(double motherMass) 
	{ 
		return motherMass - mz + 2*Constants.Proton; 
	}
	
	public	void set(double m, double i){ 
		mz = m;
		intensity = i;
	}
	public	void shiftMass(double m){ 
		mz += m;
	}
	// normalized intensity, and charge may be decided after Peak is made
	public void setNormIntensity(double normIntensity)
	{
		this.normIntensity = normIntensity;
	}
	public void setCharge(int charge)
	{
		this.charge = charge;
	}
	public Peak getShiftedPeak(double shiftMass)
	{
		Peak p = this.clone();
		p.mz += shiftMass;
		
		return p;
	}
	public int compareTo(Peak p)	// default comparator : mass
	{
		if( mz > p.mz )
			return 1;
		else if( mz < p.mz )
			return -1;		
		else
			return 0;
	}
	public static double getMassDifference(Peak p1, Peak p2)
	{
		return Math.abs(p1.mz - p2.mz);
	}
	
	public String toString() {
		return new String(Constants.getString(mz));		
	}
	
	public Peak	clone()
	{
		Peak p = new Peak(index, mz, intensity, charge);
		p.setNormIntensity(normIntensity);
		
		return p;
	}
	
	private double probability=0;
	public	double			getProbability()	{ return probability; }
	public	void			setProbability(double pa)	
	{ 
		probability= pa; 
	}
	
}

class MassComparator implements Comparator<Peak> {
	public int compare(Peak p1, Peak p2)
	{
		if(p1.mz > p2.mz) return 1;
		else if(p1.mz == p2.mz){
			if(p1.intensity > p2.intensity) return 1;
			else if(p1.intensity == p2.intensity) return 0;
			else return -1;
		}
		else return -1;	
	}
	public boolean equals(Peak p1, Peak p2)
	{
		return p1.mz == p2.mz && p1.intensity == p2.intensity;
	}
}

class IntensityComparator implements Comparator<Peak> {
	public int compare(Peak p1, Peak p2)
	{
		if( p1.intensity > p2.intensity ) return 1;
		else if(p1.intensity == p2.intensity){
			if(p1.mz > p2.mz)
				return 1;
			else if(p1.mz == p2.mz)
				return 0;
			else
				return -1;
		}
		else
			return -1;	//*/
	}
	public boolean equals(Peak p1, Peak p2)
	{
		return p1.mz == p2.mz && p1.intensity == p2.intensity;
	}
}

