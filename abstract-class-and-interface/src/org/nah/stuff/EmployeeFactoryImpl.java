package org.nah.stuff;

//public class EmployeeFactoryImpl {
// the above seems to also work just fine
// so this could be called EmployeeFactory
// and I could get rid of the other EmployeeFactory...

public class EmployeeFactoryImpl implements EmployeeFactory {

	public Employee makeEmployee(EmployeeRecord r) throws InvalidEmployeeType {

		switch(r.type) {
			case "FT":
				return new CommissionedEmployee(r);
			case "PT":
				return new HourlyEmployee(r);
			default:
				throw new InvalidEmployeeType(r.type);
		}
	}

}