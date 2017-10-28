package org.nah.stuff;

public class Driver {
	public static void main(String[] args) throws InvalidEmployeeType {
		
		EmployeeFactoryImpl employeeFactory = new EmployeeFactoryImpl();

		Employee bigBob = employeeFactory.makeEmployee(new EmployeeRecord("bob", 4, "FT"));
		
		if(bigBob instanceof HourlyEmployee) System.out.println("hourly");
		else if(bigBob instanceof CommissionedEmployee) System.out.println("comm.");
		
		System.out.println("done");
		
	}
}