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

public class Modification {
	// common
	public int Mod_Id = 0;
	public String Site = "";
	public String Position = "";
	public double Mass = 0.0;
	public int Spectrum_Count = 0;
	public int Peptide_Count = 0;
	public int location = 0; 	// index in file / location in pept / not used in table
	
	// Only in Known Mod_List_Table. 
	// In Unknown_Mod_List table, they have 'null' value.
	public String Mod_Name = "";
	public String Classification = "";
	public int static_num = 0;
	public String Reference = "";
	public String Motif = "";
	public String Fragmentation_pattern = "";
	
	// Only if parameter 'page' is -2,
	public int Total_Spectrum_Count = 0;
	
	public Modification() {
		// TODO Auto-generated constructor stub
	}
	
	public Modification(String[] str) {
		location = Byte.valueOf(str[0]).byteValue();
		Mod_Name = str[1];
		Mass = Double.valueOf(str[2]).doubleValue();
		Position = str[3];
		Site = str[4];
		Classification = str[5];			
	}

	@Override
	public String toString() {
		return "Modification [Mod_Id=" + Mod_Id + ", Site=" + Site
				+ ", Position=" + Position + ", Mass=" + Mass
				+ ", Spectrum_Count=" + Spectrum_Count + ", Peptide_Count="
				+ Peptide_Count + ", location=" + location + ", Mod_Name="
				+ Mod_Name + ", Classification=" + Classification
				+ ", static_num=" + static_num + ", Reference=" + Reference
				+ ", Motif=" + Motif + ", Fragmentation_pattern="
				+ Fragmentation_pattern + ", Total_Spectrum_Count="
				+ Total_Spectrum_Count + "]";
	}

	public Modification(String Name, String site, int fixed, int location, String position, double mass) {
		this.Mod_Name = Name;
		this.Site = site;
		this.static_num = fixed;
		this.location = location;
		this.Position = position;
		this.Mass = mass;
	}
	
//	Modification m = new Modification();
//	System.out.println(m);	// msl.models.Modification@1c95949c


//	@Override
//	public String toString() {
//		return "Modification [Mod_Id=" + Mod_Id + ", Site=" + Site
//				+ ", Position=" + Position + ", Mass=" + Mass
//				+ ", Spectrum_Count=" + Spectrum_Count + ", Peptide_Count="
//				+ Peptide_Count + ", location=" + location + ", Mod_Name="
//				+ Mod_Name + ", Classification=" + Classification
//				+ ", static_num=" + static_num + ", Reference=" + Reference
//				+ ", Motif=" + Motif + ", Fragmentation_pattern="
//				+ Fragmentation_pattern + ", Total_Spectrum_Count="
//				+ Total_Spectrum_Count + "]";
//	}

	
}
