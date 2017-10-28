package org.nah.stuff;

public abstract class Employee {
	
	public abstract boolean isPayday();
	public abstract Money calculatePay();
	public abstract void deliverPay(Money pay);

}