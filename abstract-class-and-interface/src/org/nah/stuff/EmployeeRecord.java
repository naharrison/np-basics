package org.nah.stuff;

public class EmployeeRecord {
	
	public String name;
	public int yearsExperience;
	public String type;
	
	public EmployeeRecord() {
		this.name = "ted";
		this.yearsExperience = 2;
		this.type = "FT";
	}

	public EmployeeRecord(String name, int years, String type) {
		this.name = name;
		this.yearsExperience = years;
		this.type = type;
	}

}