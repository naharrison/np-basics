package org.nah.stuff;

public interface EmployeeFactory {
	
	public Employee makeEmployee (EmployeeRecord r) throws InvalidEmployeeType;

}